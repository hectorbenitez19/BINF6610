nextflow.enable.dsl = 2

include { VALIDATE } from './modules/validate'
include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { BWA_MEM } from './modules/bwa-mem'
include { MARKDUPLICATES } from './modules/markduplicates'
include { HAPLOTYPECALLER } from './modules/haplotypecaller'
include { JOINT_GENOTYPE } from './modules/joint_genotype'
include { FILTER } from './modules/filter'
include { MULTIQC } from './modules/multiqc'
include { PUBLISH } from './modules/publish'

workflow {

    main:

    // Step 1: Load the samplesheet
    sheet = file(params.samplesheet)

    // Step 2: Create one channel item per sample
    ch_samples = channel.fromPath(sheet)
        .splitCsv(header: true)
        .map { row ->

            // Extract sample metadata
            def meta = [
                id: row.sample_id,
                single_end: row.library_type == 'single'
            ]

            // Resolve FASTQ paths relative to the samplesheet
            def r1 = sheet.parent.resolve(row.r1_fastq)

            // Single-end samples have one FASTQ; paired-end have two
            def reads = meta.single_end
                ? [r1]
                : [r1, sheet.parent.resolve(row.r2_fastq)]

            // Emit sample metadata and reads
            [meta, reads]
        }

        ref = file(params.ref)
        ref_index = files("${params.ref}.*")
        ref_dict = file("${ref.parent}/${ref.baseName}.dict")

        VALIDATE(
        sheet,
        ch_samples.map { meta, reads -> reads }.collect(),
        ref,
        ref_index,
        ref_dict,
        params.region
    )

    // Step 3: Display the channel contents
    ch_checked = ch_samples
        .combine(VALIDATE.out.sheet)
        .map { meta, reads, validated -> [meta, reads] }

    FASTQC(ch_checked)
    FASTP(ch_checked)

    BWA_MEM(
        FASTP.out.reads,
        ref,
        ref_index
    )

    MARKDUPLICATES(
        BWA_MEM.out.bam
    )

    HAPLOTYPECALLER(
    MARKDUPLICATES.out.bam,
    ref,
    ref_index,
    ref_dict,
    params.region
)

JOINT_GENOTYPE(
    HAPLOTYPECALLER.out.gvcf.collect(),
    HAPLOTYPECALLER.out.index.collect(),
    ref,
    ref_index,
    ref_dict,
    params.region
)

FILTER(
    JOINT_GENOTYPE.out.vcf,
    JOINT_GENOTYPE.out.index
)

MULTIQC(
    FASTQC.out.zip.flatten().collect(),
    FASTP.out.json.collect()
)

PUBLISH(
    VALIDATE.out.sheet,
    FILTER.out.vcf
        .mix(FILTER.out.table, MULTIQC.out.report)
        .collect()
)

publish:
      final_results = PUBLISH.out.manifest
        .mix(
            PUBLISH.out.samples,
            FILTER.out.vcf,
            FILTER.out.table,
            MULTIQC.out.report
        )

}

output {
     final_results {
        path '.'
    }

}