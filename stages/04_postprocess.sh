# shellcheck shell=bash
set -euo pipefail   # already set by the driver that sources this; stated here so the file says so itself
# 4 · postprocess — flagstat, so the numbers exist before anyone asks.
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

