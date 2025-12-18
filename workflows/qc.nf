include { FASTQC } from '../modules/qc/fastqc.nf'

workflow QC {
    take:
    rna_fastq_ch
    
    main:
    rna_fastqc_ch = rna_fastq_ch
        .transpose()
        .map { sample_id, _worksheet, _type, _referral, fastq ->
            def (_sample_id, lane, read) = (fastq.name =~ /^(.+?)_S\d+_(L\d+)_(R[12])_\d+\.fastq\.gz$/)[0][1..3]
            return tuple(sample_id, lane, read, fastq)
        }
        .groupTuple(by: [0, 1, 2])

    FASTQC(rna_fastqc_ch)
}

