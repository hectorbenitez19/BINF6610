# shellcheck shell=bash
set -euo pipefail   # already set by the driver that sources this; stated here so the file says so itself
# lib/common.sh — what every stage needs. Sourced, never executed.
#
# This file exists because there are now TWO drivers: run_pipeline.sh walks the
# whole cohort, run_sample.sh walks one sample. Both call the same ten stages,
# so the stages moved into files and the things they share moved in here.

PIPE_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die() { log "ERROR: $*"; exit 65; }

# Configuration: every value overridable from the environment, every one with a
# default, so no path is written into the code.
REF_DIR=${REF_DIR:-/courses/BINF6610.202710/data/refs/grch38}
INDEX=${INDEX:-${REF_DIR}/hisat2/genome}
GTF=${GTF:-${REF_DIR}/annotation.gtf}
FASTQ_ROOT=${FASTQ_ROOT:-/courses/BINF6610.202710/data}
CONDITION_REF=${CONDITION_REF:-Healthy}
THREADS=${THREADS:-4}

# samtools sort runs in the same pipe as hisat2 and out of the same allocation,
# so the two have to divide the cores between them rather than each taking all.
SORT_THREADS=${SORT_THREADS:-2}
ALIGN_THREADS=$(( THREADS > SORT_THREADS ? THREADS - SORT_THREADS : 1 ))


# Output layout. Derived from OUT, which the driver sets.
setup_dirs() {
QC="${OUT}/qc_raw"; TRIM="${OUT}/trim"; ALN="${OUT}/align"
LOG="${OUT}/logs"; RES="${OUT}/results"; POST="${OUT}/postprocess"; GVCF="${OUT}/gvcf"; MERGE="${OUT}/merge"
mkdir -p "$QC" "$TRIM" "$ALN" "$POST" "$LOG" "$GVCF" "$MERGE" "$RES"
}

# rows <sheet> [sample_id] -> six fields per row, whatever shape the sheet is
#
# Two jobs, and both are here so that no stage has to do either.
#
# 1. Pick columns BY NAME. Week 1's sheet had exactly six columns and every
#    stage read them positionally. The cohort sheet on the cluster has twelve,
#    and `read` puts everything left over into the LAST variable -- so reading
#    it positionally puts six extra fields inside $r2 and the run dies on a
#    filename that does not exist. Reading by name works for six or twelve.
#
# 2. Make the fastq paths absolute. The sheet stores them relative to the data
#    directory, which works when you run from there and does not when Slurm
#    runs your script from somewhere else.
#
# Given a sample id it emits that row and nothing else, so a stage body does
# not need to know which driver called it.
rows() {
    local sheet=$1 only=${2:-}
    awk -F, -v want="$only" -v root="$FASTQ_ROOT" '
        BEGIN { n = split("sample_id condition replicate library_type r1_fastq r2_fastq", need, " ") }
        NR == 1 {
            for (i = 1; i <= NF; i++) col[$i] = i
            for (i = 1; i <= n; i++)
                if (!(need[i] in col)) {
                    print "samplesheet has no column named " need[i] > "/dev/stderr"
                    exit 65
                }
            next
        }
        want != "" && $col["sample_id"] != want { next }
        {
            r1 = $col["r1_fastq"]; r2 = $col["r2_fastq"]
            if (r1 != "" && r1 !~ /^\//) r1 = root "/" r1
            if (r2 != "" && r2 !~ /^\//) r2 = root "/" r2
            print $col["sample_id"] "," $col["condition"] "," $col["replicate"] \
                  "," $col["library_type"] "," r1 "," r2
        }
    ' "$sheet"
}
