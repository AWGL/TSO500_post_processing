include { DEMULTIPLEX } from './workflows/demultiplex.nf'
include { ANALYSIS } from './workflows/analysis.nf'
include { QC } from './workflows/qc.nf'
include { DATABASE } from './workflows/database.nf'

workflow {
    sample_sheet = file(params.sample_sheet, checkIfExists: true)
    run_folder = file(params.run_folder, checkIfExists: true)
    resources = file(params.resources, checkIfExists: true)

    DEMULTIPLEX(
        sample_sheet,
        resources,
        run_folder
    )

    rna_fastq_ch = DEMULTIPLEX.out.rna_fastq_ch

    ANALYSIS(
        sample_sheet,
        resources,
        rna_fastq_ch
    )

    QC(
        rna_fastq_ch,
        ANALYSIS.out.metrics_output
    )

    DATABASE(
        ANALYSIS.out.results
    )
}
