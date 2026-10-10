process JOINT_GENOTYPE {

    tag "cohort"

    container params.containers.gatk

    cpus 2

    input:
    path gvcfs
    path indexes
    path ref
    path ref_index
    path ref_dict
    val region

    output:
    path "cohort.vcf.gz", emit: vcf
    path "cohort.vcf.gz.tbi", emit: index

    script:
    def variantArgs = gvcfs.collect { "-V ${it}" }.join(' ')

    """
    gatk CombineGVCFs \\
        -R ${ref} \\
        ${variantArgs} \\
        -O cohort.g.vcf.gz \\
        -L '${region}'

    gatk GenotypeGVCFs \\
        -R ${ref} \\
        -V cohort.g.vcf.gz \\
        -O cohort.vcf.gz \\
        -L '${region}'
    """
}