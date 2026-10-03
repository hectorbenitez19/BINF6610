# shellcheck shell=bash
set -euo pipefail   # already set by the driver that sources this; stated here so the file says so itself
# 9 · publish — say what produced this, in a file anyone can read.
#
# The course's write_manifest.sh writes results/manifest.json. On the cluster its
# platform block records the machine and the Slurm job. `sacct -j <id>` still
# answers months after the files are gone, so the job id is a permanent handle on
# "which run was this". A result you cannot trace to a job id is one you cannot
# defend. RNA-seq is not called in a region, so the fourth argument is left out.
stage_publish() {

    PIPELINE_NAME=rnaseq-de ANNOTATION="${GTF}" \
        bash "${HERE}/lib/write_manifest.sh" "${RES}" "${SHEET}" "${INDEX}"



    log "results in ${RES}:"
    ls -1 "$RES" >&2
}
