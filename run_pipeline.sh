#!/usr/bin/env bash
#-----------------------------------------------------------------------------
# rnaseq.sh — bulk RNA-seq differential expression, week 1.
#
#   bash rnaseq.sh <samplesheet.csv> <outdir> [last-stage]
#
#   bash rnaseq.sh dev.csv out validate      # stage 0 only
#   bash rnaseq.sh dev.csv out align         # stages 0..3
#   bash rnaseq.sh dev.csv out               # all ten
#
# WHY THIS FILE LOOKS LIKE THIS
#
# Every construct in here was covered in the week-1 lecture. Nothing else is
# used -- no parallelism, no resume, no temp-then-rename, no JSON, no traps,
# no associative arrays. Those all arrive in week 2, on the cluster, where
# there is a reason to want them. The one exception is the two lines after
# `set -euo pipefail`, which stage 9 needs; the Canvas page "Your pipeline's
# manifest.json" explains them.
#
# It is also deliberately repetitive: every stage reads the samplesheet with
# its own `while IFS=, read` loop. That repetition is the thing week 2 removes,
# and removing it is easier to motivate once you have felt it.
#
# The rule it does follow everywhere: after each tool runs, LOOK AT WHAT CAME
# OUT. A tool that exits 0 having produced nothing is the failure this course
# is about.
#
# WHAT TO CHECK WHILE YOU DEVELOP, AND WHAT NOT TO
#
# On a development slice -- a few thousand reads per sample -- the only question
# this pipeline can answer is "does the code run". Every assertion below is of
# that kind: did a file appear, is it non-empty, is the record count whole, is
# the alignment RATE plausible. Rates survive a slice because they are ratios.
#
# Gene counts and differential expression do not: there is no depth. So do not
# read the biology off a development run. That happens once, on the real data.
#-----------------------------------------------------------------------------
REF=${REF:-/courses/BINF6610.202710/data/refs/grch38-1000g/GRCh38_full_analysis_set_plus_decoy_hla.fa}
REGION=${REGION:-chr20:1-10000000}
THREADS=${THREADS:-4}

set -euo pipefail

# The folder this file is in, whichever folder you run it from, and when this
# run started. Stage 9 finds lib/write_manifest.sh through HERE, and the
# manifest records RUN_STARTED.
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export RUN_STARTED=$(date -u +%Y-%m-%dT%H:%M:%SZ)

# The samplesheet has SIX columns, in this order:
#   sample_id, condition, replicate, library_type, r1_fastq, r2_fastq
# Every stage reads it with `read -r id cond rep lt r1 r2`, and `read` puts
# anything left over into the LAST variable -- so a sheet with extra columns
# would silently put them all in $r2. If you add columns, read them by name.
SHEET=${1:?usage: rnaseq.sh <samplesheet.csv> <outdir> [last-stage]}
OUT=${2:?usage: rnaseq.sh <samplesheet.csv> <outdir> [last-stage]}
LAST=${3:-publish}

QC="${OUT}/qc_raw"; TRIM="${OUT}/trim"; ALN="${OUT}/align"
LOG="${OUT}/logs"; RES="${OUT}/results"; POST="${OUT}/postprocess"; GVCF="${OUT}/gvcf"; MERGE="${OUT}/merge"
mkdir -p "$QC" "$TRIM" "$ALN" "$POST" "$LOG" "$GVCF" "$MERGE" "$RES"

#--- two helpers, and they are the only ones ---------------------------------
# Messages go to stderr so that a stage's stdout stays free for data.
log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die()  { printf 'error: %s\n' "$*" >&2; exit 65; }

# The ten stages, in the order they run. Each name has a matching function
# below: `validate` -> stage_validate, `align` -> stage_align, and so on.
STAGES=(validate qc_raw trim align postprocess quantify merge analyze qc_report publish)

# Catch a typo in the third argument before running anything.
known=0
for stage in "${STAGES[@]}"; do
    [[ "$stage" == "$LAST" ]] && known=1
done
(( known )) || die "unknown stage: ${LAST}"

#=============================================================================
# 0 · validate — check everything before computing anything
#=============================================================================
stage_validate() {
    local id cond rep lt r1 r2 problems=0 n1 n2 r1_ok r2_ok local dict="${REF%.*}.dict"

    while IFS=, read -r id cond rep lt r1 r2; do
        [[ -n "$id" ]] || { log "a row has no sample_id"; problems=$(( problems + 1 )); continue; }

        # the files exist and are not empty
        [[ -s "$r1" ]] || { log "$id: R1 missing or empty: $r1"; problems=$(( problems + 1 )); }
        if [[ "$lt" == paired ]]; then
            [[ -s "$r2" ]] || { log "$id: declared paired but R2 is missing"; problems=$(( problems + 1 )); }
        fi

        # the r2 column and the declared layout agree
        if [[ -z "$r2" && "$lt" == paired ]]; then
            log "$id: r2_fastq is empty but library_type says paired"
            problems=$(( problems + 1 ))
        fi

        # the gzip streams are whole -- both of them. A stream that fails here is
        # not read below: under pipefail, `gzip -dc` on it would end the script
        # instead of letting it report every problem.
        r1_ok=0 r2_ok=0
        if [[ -s "$r1" ]]; then
            if gzip -t "$r1" 2>/dev/null; then r1_ok=1
            else log "$id: R1 is not a valid gzip file"; problems=$(( problems + 1 )); fi
        fi
        if [[ "$lt" == paired && -s "$r2" ]]; then
            if gzip -t "$r2" 2>/dev/null; then r2_ok=1
            else log "$id: R2 is not a valid gzip file"; problems=$(( problems + 1 )); fi
        fi

        # the records are whole, and the mates agree
        if (( r1_ok )); then
            n1=$(gzip -dc "$r1" | wc -l)
            (( n1 % 4 == 0 )) || { log "$id: R1 has $n1 lines, not a whole number of records"
                                   problems=$(( problems + 1 )); }
            if (( r2_ok )); then
                n2=$(gzip -dc "$r2" | wc -l)
                (( n1 == n2 )) || { log "$id: R1 has $(( n1 / 4 )) reads, R2 has $(( n2 / 4 ))"
                                    problems=$(( problems + 1 )); }
            fi
        fi
    done < <(tail -n +2 "$SHEET")

    # duplicate sample ids. sort | uniq -d prints only the repeats.
    local dupes
    dupes=$(awk -F, 'NR>1 { print $1 }' "$SHEET" | sort | uniq -d)
    [[ -z "$dupes" ]] || { log "duplicate sample_id: $dupes"; problems=$(( problems + 1 )); }

    # the reference is where the config says it is
[[ -s "${REF}.fai" ]] || {
    log "reference FASTA index missing: ${REF}.fai"
    problems=$(( problems + 1 ))
}

# check if the reference dictionary exists
[[ -s "$dict" ]] || {
    log "reference dictionary missing: $dict"
    problems=$(( problems + 1 ))
}

# check the bwa index exists
for ext in amb ann bwt pac sa; do
    [[ -s "${REF}.${ext}" ]] || {
        log "BWA index missing: ${REF}.${ext}"
        problems=$(( problems + 1 ))
    }
done

    (( problems == 0 )) || die "validation failed with ${problems} problem(s)"
    log "validation passed"
}

#=============================================================================
# 1 · qc_raw — FastQC on the reads as they arrived
#=============================================================================
stage_qc_raw() {
    local id cond rep lt r1 r2 base
    while IFS=, read -r id cond rep lt r1 r2; do
        # stdout as well as stderr: fastqc prints "application/gzip" per file on
        # STDOUT, which otherwise lands in the terminal and looks like output.
        fastqc -q -o "$QC" "$r1" > "${LOG}/${id}.fastqc.log" 2>&1
        [[ "$lt" != paired ]] || fastqc -q -o "$QC" "$r2" >> "${LOG}/${id}.fastqc.log" 2>&1

        # fastqc can exit 0 and write nothing. Ask the disk.
        base=$(basename "$r1" .fastq.gz)
        [[ -s "${QC}/${base}_fastqc.zip" ]] || die "$id: fastqc produced no report"
        log "$id: qc done"
    done < <(tail -n +2 "$SHEET")
}

#=============================================================================
# 2 · trim — adapters and low-quality tails
#=============================================================================
stage_trim() {
    local id cond rep lt r1 r2 n
    while IFS=, read -r id cond rep lt r1 r2; do
        if [[ "$lt" == paired ]]; then
            fastp -i "$r1" -I "$r2" \
                  -o "${TRIM}/${id}_R1.fastq.gz" -O "${TRIM}/${id}_R2.fastq.gz" \
                  -j "${LOG}/${id}.fastp.json" -h "${LOG}/${id}.fastp.html" \
                  2> "${LOG}/${id}.fastp.log"
        else
            fastp -i "$r1" -o "${TRIM}/${id}_R1.fastq.gz" \
                  -j "${LOG}/${id}.fastp.json" -h "${LOG}/${id}.fastp.html" \
                  2> "${LOG}/${id}.fastp.log"
        fi

        # trimming can only remove reads. Zero left means something is wrong.
        n=$(gzip -dc "${TRIM}/${id}_R1.fastq.gz" | wc -l)
        (( n > 0 )) || die "$id: nothing survived trimming"
        log "$id: trimmed to $(( n / 4 )) reads"
    done < <(tail -n +2 "$SHEET")
}

#=============================================================================
# 3 · align — BWA straight into a sorted MEM
#=============================================================================
stage_align() {
    local id cond rep lt r1 r2 rate
    while IFS=, read -r id cond rep lt r1 r2; do
        trimmed_r1="${TRIM}/${id}_R1.fastq.gz"
        trimmed_r2="${TRIM}/${id}_R2.fastq.gz"
        bam="${ALN}/${id}.bam"

        log "$id: aligning reads with BWA"

        if [[ "$lt" == paired ]]; then

            bwa mem \
                -t "$THREADS" \
                -R "@RG\tID:${id}\tSM:${id}" \
                "$REF" \
                "$trimmed_r1" \
                "$trimmed_r2" |
                samtools view -b -o "$bam" -

        else

            bwa mem \
                -t "$THREADS" \
                -R "@RG\tID:${id}\tSM:${id}" \
                "$REF" \
                "$trimmed_r1" |
                samtools view -b -o "$bam" -

        fi

        [[ -s "$bam" ]] ||
            die "$id: alignment produced no BAM"

        log "$id: alignment complete"
    done < <(tail -n +2 "$SHEET")
}

#=============================================================================
# 4 · postprocess — flagstat, so the numbers exist before anyone asks
#=============================================================================
stage_postprocess() {
    local id cond rep lt r1 r2
    while IFS=, read -r id cond rep lt r1 r2; do
        bam="${ALN}/${id}.bam"
        sorted_bam="${POST}/${id}.sorted.bam"
        dedup_bam="${POST}/${id}.dedup.bam"
        metrics="${POST}/${id}.duplicate_metrics.txt"

        # Make sure Stage 3 produced the expected BAM.
        [[ -s "$bam" ]] ||
            die "$id: aligned BAM missing: $bam"

        # Sort by genomic coordinate.
        log "$id: sorting BAM"

        samtools sort \
            -@ "$THREADS" \
            -o "$sorted_bam" \
            "$bam"

        [[ -s "$sorted_bam" ]] ||
            die "$id: samtools sort produced no BAM"

        # Index the sorted BAM.
        log "$id: indexing sorted BAM"

        samtools index "$sorted_bam"

        [[ -s "${sorted_bam}.bai" ]] ||
            die "$id: failed to index sorted BAM"

        # Mark PCR/optical duplicates.
        log "$id: marking duplicates"

        gatk MarkDuplicates \
            -I "$sorted_bam" \
            -O "$dedup_bam" \
            -M "$metrics"

        [[ -s "$dedup_bam" ]] ||
            die "$id: MarkDuplicates produced no BAM"

        [[ -s "$metrics" ]] ||
            die "$id: MarkDuplicates produced no metrics file"

        # Index the final BAM used for variant calling.
        log "$id: indexing deduplicated BAM"

        samtools index "$dedup_bam"

        [[ -s "${dedup_bam}.bai" ]] ||
            die "$id: failed to index deduplicated BAM"

        log "$id: postprocessing complete"
    done < <(tail -n +2 "$SHEET")
}

#=============================================================================
# 5 · quantify — featureCounts, one sample at a time
#=============================================================================
stage_quantify() {
    local id cond rep lt r1 r2 flag genes
    while IFS=, read -r id cond rep lt r1 r2; do
        bam="${POST}/${id}.dedup.bam"
        gvcf="${GVCF}/${id}.g.vcf.gz"

        # Make sure Stage 4 produced the required BAM and index.
        [[ -s "$bam" ]] ||
            die "$id: deduplicated BAM missing: $bam"

        [[ -s "${bam}.bai" ]] ||
            die "$id: BAM index missing: ${bam}.bai"

        log "$id: calling variants with GATK HaplotypeCaller"

        gatk HaplotypeCaller \
            -R "$REF" \
            -I "$bam" \
            -O "$gvcf" \
            -L "$REGION" \
            -ERC GVCF

        [[ -s "$gvcf" ]] ||
            die "$id: HaplotypeCaller produced no GVCF"

        log "$id: GVCF complete"

    done < <(tail -n +2 "$SHEET")
}

#=============================================================================
# 6 · merge — one matrix for the whole cohort
#
# THIS IS THE BARRIER. It cannot start until every sample has finished stage 5,
# and the only reason it works here is that the loop above ran to completion.
# Week 2 makes that a declaration instead of a line position.
#=============================================================================
stage_merge() {
    local id cond rep lt r1 r2
    local gvcf
    local combined="${MERGE}/cohort.g.vcf.gz"
    local cohort="${MERGE}/cohort.vcf.gz"
    local -a gvcf_args=()

    # Collect all per-sample GVCFs.
    while IFS=, read -r id cond rep lt r1 r2; do

        gvcf="${GVCF}/${id}.g.vcf.gz"

        [[ -s "$gvcf" ]] ||
            die "$id: GVCF missing: $gvcf"

        [[ -s "${gvcf}.tbi" ]] ||
            die "$id: GVCF index missing: ${gvcf}.tbi"

        gvcf_args+=("-V" "$gvcf")

    done < <(tail -n +2 "$SHEET")

    # Combine the per-sample GVCFs.
    log "combining per-sample GVCFs"

    gatk CombineGVCFs \
        -R "$REF" \
        "${gvcf_args[@]}" \
        -O "$combined"

    [[ -s "$combined" ]] ||
        die "CombineGVCFs produced no cohort GVCF"

    # Joint genotype all samples.
    log "joint genotyping cohort"

    gatk GenotypeGVCFs \
        -R "$REF" \
        -V "$combined" \
        -O "$cohort"

    [[ -s "$cohort" ]] ||
        die "GenotypeGVCFs produced no cohort VCF"

    log "joint genotyping complete"
}

#=============================================================================
# 7 · analyze — DESeq2, in R, called with four arguments
#=============================================================================
stage_analyze() {
    local input="${MERGE}/cohort.vcf.gz"
    local output="${RES}/cohort.filtered.vcf.gz"

    [[ -s "$input" ]] ||
        die "cohort VCF missing: $input"

    log "applying hard filters to cohort VCF"

       gatk VariantFiltration \
        -R "$REF" \
        -V "$input" \
        -O "$output" \

    [[ -s "$output" ]] ||
        die "VariantFiltration produced no filtered VCF"

    log "variant filtering complete"
}

#=============================================================================
# 8 · qc_report — MultiQC over every log this run produced
#=============================================================================
stage_qc_report() {
    multiqc -q -o "$RES" "$QC" "$LOG" 2> "${LOG}/multiqc.log"
    [[ -s "${RES}/multiqc_report.html" ]] || die "multiqc produced no report"
    log "QC report written"
}

#=============================================================================
# 9 · publish — say what produced this, in a file anyone can read
#=============================================================================
stage_publish() {
    # The course's script writes results/manifest.json: the commit this code was
    # on, the samples, the reference and annotation, the machine, every published
    # file with its checksum, and each tool's version. RNA-seq is not called in a
    # region, so the fourth argument is left out.
    PIPELINE_NAME=rnaseq-de ANNOTATION="${GTF}" \
        bash "${HERE}/lib/write_manifest.sh" "${RES}" "${SHEET}" "${INDEX}"

    log "results in ${RES}:"
    ls -1 "$RES" >&2
}

#=============================================================================
# the driver — ten stages, in order, one after another
#=============================================================================
# Run the stages in order and stop after the one named on the command line.
# `"stage_${stage}"` calls the function whose name is built from the stage name.
n=0
for stage in "${STAGES[@]}"; do
    log "===== stage ${n} : ${stage} ====="
    "stage_${stage}"
    [[ "$stage" == "$LAST" ]] && break
    n=$(( n + 1 ))
done
log "done"
