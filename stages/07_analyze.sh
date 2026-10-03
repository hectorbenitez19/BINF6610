# shellcheck shell=bash
set -euo pipefail   # already set by the driver that sources this; stated here so the file says so itself
# 7 · analyze — DESeq2, in R, called with four arguments.
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