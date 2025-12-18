include { FASTQC } from '../modules/qc/fastqc.nf'
include { CREATE_SAMPLE_QC_FILE } from '../modules/qc/create_sample_qc_file.nf'
include { MERGE_QC_FILES } from '../modules/qc/merge_qc_files.nf'

workflow QC {
    take:
    rna_fastq_ch
    metrics_output

    main:
    // Turn a channel like [sample1, [sample1_L001_R2.fastq.gz, sample1_L001_R1.fastq.gz...]]
    // in to [[sample1, sample1_L001_R1.fastq.gz], [sample1, sample1_L001_R2.fastq.gz]] then
    // group by sample_id, lane and read. The same lane and read can have multiple fastq files
    rna_fastqc_ch = rna_fastq_ch
        .transpose()
        .map { sample_id, _worksheet, _type, _referral, fastq ->
            def (_sample_id, lane, read) = (fastq.name =~ /^(.+?)_S\d+_(L\d+)_(R[12])_\d+\.fastq\.gz$/)[0][1..3]
            return tuple(sample_id, lane, read, fastq)
        }
        .groupTuple(by: [0, 1, 2])

    FASTQC(rna_fastqc_ch)

    // Join fastqc results to app analysis metrics on first element of tuple (sample id)
    joined_qc_ch = FASTQC.out.fastqc_summary.join(metrics_output)

    CREATE_SAMPLE_QC_FILE(joined_qc_ch)

    // Group sample by worklist 
    worksheet_qc_channel = CREATE_SAMPLE_QC_FILE.out.sample_qc_file.groupTuple(by: 1)

    MERGE_QC_FILES(worksheet_qc_channel)
}
