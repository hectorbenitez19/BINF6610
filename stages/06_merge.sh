# shellcheck shell=bash
set -euo pipefail   # already set by the driver that sources this; stated here so the file says so itself
# 6 · merge — one matrix for the whole cohort.
#
# THIS IS THE BARRIER. It cannot start until every sample has finished stage 5.
# Last week that worked because the loop above it had run to completion. On the
# cluster the twelve samples are twelve separate jobs, so the barrier is a
# declaration -- `--dependency=afterok` -- and this stage only ever runs in the
# cohort job.
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
