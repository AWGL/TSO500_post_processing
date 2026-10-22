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
        run_folder,
    )

    rna_fastq_ch = DEMULTIPLEX.out.rna_fastq_ch

    ANALYSIS(
        sample_sheet,
        resources,
        rna_fastq_ch,
    )

    QC(
        rna_fastq_ch,
        ANALYSIS.out.metrics_output,
    )

    ntc_bam_ch = ANALYSIS.out.bams.filter { sample_id, _worksheet, _bams ->
        {
            sample_id.contains("NTC")
        }
    }

    DATABASE(
        ANALYSIS.out.results,
        ntc_bam_ch,
        QC.out.sample_qc,
        params.run_id,
        DEMULTIPLEX.out.rna_samples_list,
    )
}
