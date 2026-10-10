process FASTQC {

    tag "${meta.id}"

    container params.containers.fastqc

    cpus 2

    input:
    tuple val(meta), path(reads)

    output:
    path "*_fastqc.zip", emit: zip

    script:
    """
    fastqc \
        --threads ${task.cpus} \
        --outdir . \
        ${reads.join(' ')}
    """
}
