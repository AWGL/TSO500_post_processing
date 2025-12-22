include { APP_PIPELINE } from '../modules/tso500_app/app.nf'

workflow ANALYSIS {
    take:
    sample_sheet
    resources
    rna_fastq_ch

    main:
    APP_PIPELINE(
        sample_sheet,
        resources,
        rna_fastq_ch,
    )

    emit:
    results = APP_PIPELINE.out.results
    metrics_output = APP_PIPELINE.out.metrics_output
    bams = APP_PIPELINE.out.bams
}
