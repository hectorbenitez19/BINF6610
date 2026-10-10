process FILTER {

    tag "cohort"

    container params.containers.gatk

    cpus 2

    input:
    path vcf
    path vcf_index

    output:
    path "cohort.filtered.vcf.gz", emit: vcf
    path "variants.tsv", emit: table

    script:
    """
    gatk VariantFiltration \\
        -V ${vcf} \\
        -O cohort.filtered.vcf.gz \\
        --filter-name "LowQD" \\
        --filter-expression "QD < 2.0" \\
        --filter-name "HighFS" \\
        --filter-expression "FS > 60.0" \\
        --filter-name "LowMQ" \\
        --filter-expression "MQ < 40.0"

    bcftools query \\
        -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\n' \\
        cohort.filtered.vcf.gz \\
        >> variants.tsv
    """
}
