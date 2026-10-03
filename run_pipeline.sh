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
SHEET=${1:?usage: run_pipeline.sh <samplesheet.csv> <outdir> [last-stage]}
OUT=${2:?usage: run_pipeline.sh <samplesheet.csv> <outdir> [last-stage]}
LAST=${3:-publish}

QC="${OUT}/qc_raw"
TRIM="${OUT}/trim"
ALN="${OUT}/align"
POST="${OUT}/postprocess"
GVCF="${OUT}/gvcf"
MERGE="${OUT}/merge"
LOG="${OUT}/logs"
RES="${OUT}/results"

mkdir -p "$QC" "$TRIM" "$ALN" "$POST" "$LOG" "$GVCF" "$MERGE" "$RES"

source "${HERE}/lib/common.sh"

STAGES=(validate qc_raw trim align postprocess quantify merge analyze qc_report publish)

known=0

for stage_file in "${HERE}"/stages/*.sh; do
    source "$stage_file"
done

# Catch a typo in the third argument before running anything.
known=0
for stage in "${STAGES[@]}"; do
    [[ "$stage" == "$LAST" ]] && known=1
done
(( known )) || die "unknown stage: ${LAST}"


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
