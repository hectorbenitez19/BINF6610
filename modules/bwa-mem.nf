process BWA_MEM {

    tag "${meta.id}"

    container params.containers.bwa

    cpus 2

    input:
    tuple val(meta), path(reads)
    path ref
    path ref_index

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam

    script:
    def readArgs = reads.join(' ')
    def readGroup = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:${meta.id}"

    """
    bwa mem \
        -t ${task.cpus} \
        -R '${readGroup}' \\
        ${ref} \
        ${readArgs} \
        | samtools sort \
            -@ ${task.cpus} \
            -o ${meta.id}.sorted.bam \
            -
    """
}