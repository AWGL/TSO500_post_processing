include { APP_DEMULTIPLEX } from '../modules/tso500_app/demultiplex.nf'
include { FILTER_SAMPLE_SHEET } from '../modules/sample_sheet/filtering.nf'

workflow DEMULTIPLEX {
    take:
    sample_sheet
    resources
    run_folder

    main:                                                                                   
    APP_DEMULTIPLEX(
        sample_sheet,
        resources,
        run_folder,
    )

    FILTER_SAMPLE_SHEET(
        sample_sheet
    )

    // Create a channel of fastqs grouped to the sample ID
    fastq_ch = APP_DEMULTIPLEX.out.fastq
        .flatten()
        .map { filename ->
            {
                def sample_id = filename.name.split('_')[0]
                return tuple(sample_id, filename)
            }
        }
        .groupTuple()

    // Get only RNA fastqs
    rna_fastq_ch = FILTER_SAMPLE_SHEET.out.rna_sample_list
        .splitCsv()
        .combine(fastq_ch, by: 0)

    emit:
    rna_fastq_ch = rna_fastq_ch
}
