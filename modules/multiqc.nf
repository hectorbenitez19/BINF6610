process MULTIQC {

    tag "cohort"

    container params.containers.multiqc

    cpus 1

    input:
    path fastqc_reports
    path fastp_reports

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --force
    """
}
