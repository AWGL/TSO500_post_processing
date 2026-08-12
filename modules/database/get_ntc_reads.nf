process GET_NTC_READS {
    tag "${sample_id}"

    label 'samtools_container'

    input:
    tuple val(sample_id), val(worksheet), path(ntc_bam)

    output:
    tuple val(sample_id), val(worksheet), stdout, emit: ntc_reads

    script:
    """
    samtools view -F4 -c ${ntc_bam} 
    """
}
