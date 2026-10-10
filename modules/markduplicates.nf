process MARKDUPLICATES {

    tag "${meta.id}"

    container params.containers.gatk

    cpus 2

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta),
          path("${meta.id}.dedup.bam"),
          path("${meta.id}.dedup.bam.bai"),
          emit: bam

    script:
    """
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${meta.id}.dedup.bam \\
        -M ${meta.id}.duplicate_metrics.txt \\
        --CREATE_INDEX false

    samtools index \\
        -@ ${task.cpus} \\
        ${meta.id}.dedup.bam
    """
}