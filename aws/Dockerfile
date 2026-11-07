
FROM mambaorg/micromamba:2.0.5-ubuntu24.04

# Install the pinned bioinformatics tools from Assignment 3
RUN micromamba install -y -n base \
    -c conda-forge -c bioconda \
    bwa=0.7.19 \
    samtools=1.24 \
    bcftools=1.24 \
    gatk4=4.6.2.0 \
    fastqc=0.12.1 \
    fastp=1.3.7 \
    multiqc=1.35 \
    git=2.47.1 \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH

# Add AWS CLI from the official AWS CLI container image
COPY --from=public.ecr.aws/aws-cli/aws-cli:2.37.12 \
    /usr/local/aws-cli /usr/local/aws-cli

# Create a command for the AWS CLI
USER root

RUN ln -s /usr/local/aws-cli/v2/current/bin/aws /usr/local/bin/aws

# Return to the micromamba non-root user
USER $MAMBA_USER

# Copy the entire pipeline into the container
COPY --chown=$MAMBA_USER:$MAMBA_USER . /pipeline/

# Copy the AWS Batch entry script
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    aws/batch_run.sh /usr/local/bin/batch_run.sh
