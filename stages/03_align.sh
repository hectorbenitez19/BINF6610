# shellcheck shell=bash
set -euo pipefail   # already set by the driver that sources this; stated here so the file says so itself
# 3 · align — HISAT2 straight into a sorted BAM.
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
