# shellcheck shell=bash
set -euo pipefail   # already set by the driver that sources this; stated here so the file says so itself
# 8 · qc_report — MultiQC over every log this run produced.
stage_qc_report() {
    multiqc -q -o "$RES" "$QC" "$LOG" 2> "${LOG}/multiqc.log"
    [[ -s "${RES}/multiqc_report.html" ]] || die "multiqc produced no report"
    log "QC report written"
}
