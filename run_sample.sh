#!/usr/bin/env bash
#=============================================================================
# run_sample.sh — stages 0 to 5, for ONE sample.
#
#   bash run_sample.sh <samplesheet.csv> <outdir> <sample_id> [last-stage]
#
# This is what a Slurm array task runs: the scheduler picks the sample, this
# runs the pipeline on it. It calls exactly the same stage functions as
# run_pipeline.sh — which is why the stages live in stages/ and not inside one
# of the two drivers.
#
# It refuses to go past stage 5 on purpose. Stages 6 to 9 need every sample, and
# a single task cannot know whether the other eleven have finished. Deciding
# that is the scheduler's job, not this script's.
#=============================================================================
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

SHEET=${1:?usage: run_sample.sh <samplesheet.csv> <outdir> <sample_id> [last-stage]}
OUT=${2:?usage: run_sample.sh <samplesheet.csv> <outdir> <sample_id> [last-stage]}
SAMPLE=${3:?usage: run_sample.sh <samplesheet.csv> <outdir> <sample_id> [last-stage]}
LAST=${4:-quantify}

source "${HERE}/lib/common.sh"
setup_dirs
for f in "${HERE}"/stages/*.sh; do source "$f"; done

# The sample has to be in the sheet. Without this, a typo runs zero samples
# through every stage and exits 0 -- a task that did nothing and reported success.
[[ -n "$(rows "$SHEET" "$SAMPLE")" ]] || die "no sample '${SAMPLE}' in ${SHEET}"

PER_SAMPLE=(validate qc_raw trim align postprocess quantify)
known=0
for stage in "${PER_SAMPLE[@]}"; do [[ "$stage" == "$LAST" ]] && known=1; done
(( known )) || die "run_sample.sh stops at quantify; '${LAST}' needs the whole cohort"

n=0
for stage in "${PER_SAMPLE[@]}"; do
    log "===== ${SAMPLE} · stage ${n} : ${stage} ====="
    "stage_${stage}"
    [[ "$stage" == "$LAST" ]] && break
    n=$(( n + 1 ))
done
log "${SAMPLE}: done"
