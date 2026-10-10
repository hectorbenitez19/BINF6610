
process VALIDATE {

    tag "samplesheet"

    container params.containers.tools

    input:
    path samplesheet
    path fastqs
    path ref
    path ref_index
    path ref_dict
    val region

    output:
    path "samplesheet.csv", emit: sheet

    script:
    """

    validate_samplesheet.sh \
        samplesheet.csv \
        ${ref} \
        '${region}'
    """
}
