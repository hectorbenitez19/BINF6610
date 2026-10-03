#!/bin/bash
#SBATCH --job-name=w2demo-persample
#SBATCH --array=1-12
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=logs/persample_%A_%a.out
#SBATCH --export=NONE
#-----------------------------------------------------------------------------
# Stages 0-5 for ONE sample. Twelve of these run at the same time.
#
# There is no bioinformatics in this file. It decides WHICH sample runs and with
# what resources; the pipeline decides what to do to it.
#
# The numbers above were measured, not guessed. Briefly, and in full in
# ../RESOURCES.md:
#   -c 8     the knee. 16 cores saves 39 s per sample and costs 20 extra
#            core-minutes; 28 cores is SLOWER than 16.
#   --mem=16G  cgroup high-water mark measured twice at 6.76 GB, of which the
#            4.2 GB HISAT2 index is most. 16G is 2.4x, which is the margin for
#            the first sample that is deeper than this cohort.
#   --time     longest of twelve was 3:54 with all twelve running at once.
#-----------------------------------------------------------------------------
set -euo pipefail
# NOT $(dirname "$0"): sbatch copies your script to /var/spool/slurmd and runs
# the COPY, so $0 is there and conf/ is not beside it. SLURM_SUBMIT_DIR is the
# directory you ran sbatch from, and it survives --export=NONE because Slurm
# sets it itself.
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

# The threads Slurm actually gave us, handed to the tools. The pipeline already
# reads THREADS from the environment, so this is the whole change.
export THREADS="${SLURM_CPUS_PER_TASK}"
export REF_DIR FASTQ_ROOT

# Scratch on the node's own disk, and a trap that removes it however this job
# ends. /tmp here is shared with every other job on the node, and a `rm` at the
# bottom of the script does not run when the job is cancelled or times out.
export TMPDIR="/tmp/${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"
trap 'rm -rf "${TMPDIR}"' EXIT

# One task is one sample. The scheduler picks WHICH; run_sample.sh does the work.
SAMPLE=$(awk -F, -v n="${SLURM_ARRAY_TASK_ID}" 'NR == n + 1 { print $1 }' "${SHEET}")

# An array wider than the samplesheet gets an empty name here, and a pipeline
# handed nothing succeeds at every stage. Refuse instead.
[[ -n "${SAMPLE}" ]] || { echo "task ${SLURM_ARRAY_TASK_ID}: no row ${SLURM_ARRAY_TASK_ID} in ${SHEET}" >&2; exit 64; }

echo "task ${SLURM_ARRAY_TASK_ID} on $(hostname): ${SAMPLE}, ${SLURM_CPUS_PER_TASK} cores"
bash "${PIPELINE_DIR}/run_sample.sh" "${SHEET}" "${RUN_ROOT}/run" "${SAMPLE}" quantify
