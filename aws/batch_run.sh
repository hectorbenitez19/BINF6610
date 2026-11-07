#!/usr/bin/env bash
set -euo pipefail

# Required AWS Batch environment variables
: "${DATA:?set DATA to the S3 input directory}"
: "${RESULTS:?set RESULTS to the S3 results directory}"

# Default reference genome
: "${REF_NAME:=smoke.fa}"

# AWS Batch provides the job ID automatically
job=${AWS_BATCH_JOB_ID:-local}

# Create a temporary working directory
work=$(mktemp -d /tmp/run.XXXXXX)

# Step 1: Download the smoke dataset from S3
echo "[batch_run] 1/3 copying ${DATA}/ to ${work}/data"

aws s3 cp \
    "${DATA%/}/" \
    "${work}/data/" \
    --recursive \
    --no-progress

# Step 2: Run the bioinformatics pipeline
echo "[batch_run] 2/3 running the pipeline"
echo "REF=${REF_NAME} REGION=${REGION:-?} THREADS=${THREADS:-?}"

# Run from the data directory so relative FASTQ paths resolve
cd "${work}/data"

REF="${work}/data/${REF_NAME}" \
bash /pipeline/run_pipeline.sh \
    samplesheet.csv \
    "${work}/out"

# Step 3: Upload the results to S3
echo "[batch_run] 3/3 copying results to ${RESULTS%/}/${job}/"

aws s3 cp \
    "${work}/out/results/" \
    "${RESULTS%/}/${job}/" \
    --recursive \
    --no-progress

echo "[batch_run] done"
