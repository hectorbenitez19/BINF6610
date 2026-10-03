#!/bin/bash
#SBATCH --job-name=w2demo-cohort
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/cohort_%j.out
#SBATCH --export=NONE
#-----------------------------------------------------------------------------
# Stages 6-9, once, after every array task has succeeded.
#
# It runs the WHOLE pipeline with the FULL samplesheet, starting at stage 0.
# Stages 1-5 find their outputs already on disk and skip, which costs about a
# second for all twelve samples -- so there is no "start from stage 6" argument
# to get wrong, and a rerun after a partial failure does the right thing by
# itself.
#
# Measured: 38 seconds, 0.95 cores busy, cgroup high-water 0.77 GB. So -c 2 and
# 4G -- merge, DESeq2 and MultiQC are all single-threaded and none of them is
# large. The first version of this file asked for 4 cores and 8G because nobody
# had measured; it was using a quarter of the cores and reporting 79 % of the
# memory, and the 79 % turned out to be the samplesheet's fastq passing through
# the page cache in stage 0, not anything DESeq2 needed.
#
# --time is 20 minutes rather than the 1.5x of 38 seconds that the measurement
# suggests, on purpose: this job's cost grows with the number of samples and
# with how much MultiQC has to read, and a TIMEOUT here throws away every array
# task that succeeded before it.
#-----------------------------------------------------------------------------
set -euo pipefail
# WHERE THIS MUST BE SUBMITTED FROM, and what happens if you get it wrong.
# $SLURM_SUBMIT_DIR is the directory you were standing in when you typed sbatch,
# so conf/ has to be there. `bash submit.sh` guarantees it by cd-ing here first;
# submitting this file by hand from somewhere else does not. Measured: from /tmp
# the job comes back FAILED 1:0 and there is NO log to read, because
# --output=logs/... is also relative to the submit directory and /tmp/logs does
# not exist. So test for the file and say so, rather than letting `source` fail.
if [[ ! -s "${SLURM_SUBMIT_DIR}/conf/slurm.env" ]]; then
    printf 'no conf/slurm.env in %s\n' "${SLURM_SUBMIT_DIR}" >&2
    printf 'submit from the slurm/ directory, or use: bash submit.sh\n' >&2
    exit 78            # EX_CONFIG
fi
source "${SLURM_SUBMIT_DIR}/conf/slurm.env"

module load miniconda3/25.9.1
source activate "${CONDA_ENV}"

export THREADS="${SLURM_CPUS_PER_TASK}"
export REF_DIR FASTQ_ROOT
export TMPDIR="/tmp/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"
trap 'rm -rf "${TMPDIR}"' EXIT

echo "cohort on $(hostname), ${SLURM_CPUS_PER_TASK} cores"
bash "${PIPELINE_DIR}/run_pipeline.sh" "${SHEET}" "${RUN_ROOT}/run" publish

# The results are small and must outlive the /scratch purge. The BAMs are not
# and stay where they are.
mkdir -p "${KEEP_DIR}"
cp -r "${RUN_ROOT}/run/results" "${KEEP_DIR}/results-${SLURM_JOB_ID}"
echo "kept: ${KEEP_DIR}/results-${SLURM_JOB_ID}"
