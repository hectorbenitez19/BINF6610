# shellcheck shell=bash
set -euo pipefail   # already set by the driver that sources this; stated here so the file says so itself
# 5 · quantify — featureCounts, one sample at a time.
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