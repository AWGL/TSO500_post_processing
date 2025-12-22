process GET_NTC_READS {
    tag "${sample_id}"

    container "community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c"

    input:
    tuple val(sample_id), val(worksheet), path(ntc_bam)

    output:
    tuple val(sample_id), val(worksheet), stdout, emit: ntc_reads

    script:
    """
    samtools view -F4 -c ${ntc_bam} 
    """
}
