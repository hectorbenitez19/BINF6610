process FASTP {

    tag "${meta.id}"

    container params.containers.fastp

    cpus 2

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    path "${meta.id}.fastp.json", emit: json

    script:
    def r1 = reads[0]
    def r2 = meta.single_end ? null : reads[1]

    def trimming = meta.single_end
        ? "--in1 ${r1} --out1 ${meta.id}_R1.trim.fastq.gz"
        : "--in1 ${r1} --in2 ${r2} --out1 ${meta.id}_R1.trim.fastq.gz --out2 ${meta.id}_R2.trim.fastq.gz"

    """
    fastp \
        ${trimming} \
        --thread ${task.cpus} \
        --html ${meta.id}.fastp.html \
        --json ${meta.id}.fastp.json
    """
}