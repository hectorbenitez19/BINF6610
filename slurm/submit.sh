#!/usr/bin/env bash
#-----------------------------------------------------------------------------
# submit.sh — two submissions, thirteen jobs, one barrier.
#
#   bash submit.sh              all twelve samples
#   bash submit.sh 1-2          the first two, to check it works
#-----------------------------------------------------------------------------
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$HERE"
source conf/slurm.env

RANGE=${1:-}

# Slurm will not create this, and a task that cannot open its output file fails
# with nowhere to write why.
mkdir -p logs

ARRAY_ARGS=()
[[ -z "$RANGE" ]] || ARRAY_ARGS=(--array="$RANGE")

# --parsable prints only the job id, so it goes straight into a variable.
ARRAY_ID=$(sbatch --parsable -p "$PARTITION" -A "$ACCOUNT" "${ARRAY_ARGS[@]}" 01_persample.sbatch)
echo "array   ${ARRAY_ID}"

# afterok, not afterany: with afterany a failed sample gives you a cohort matrix
# that is missing a column and says nothing about it.
COHORT_ID=$(sbatch --parsable -p "$PARTITION" -A "$ACCOUNT" \
                   --dependency=afterok:"${ARRAY_ID}" --kill-on-invalid-dep=yes \
                   02_cohort.sbatch)
echo "cohort  ${COHORT_ID}   (waits for ${ARRAY_ID})"
echo
echo "watch:   squeue -u \$USER"
echo "measure: seff ${ARRAY_ID}_1 ; seff ${COHORT_ID}"
