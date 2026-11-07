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
        --filter-name "QD20" \
        --filter-expression "QD < 2.0"

        

    [[ -s "$output" ]] ||
        die "VariantFiltration produced no filtered VCF"

    {
        printf 'CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\n'
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\n' \
            "$output"
    } > "${RES}/variants.tsv"

    log "variant filtering complete"

    log "caught by each filter: $(awk 'NR>1 && $6!="PASS" {print $6}' "${RES}/variants.tsv" | sort | uniq -c | awk '{printf "%s %s  ", $2, $1}')"
}