process MERGE_SAMPLES_DATABASE {
    cpus 1
    memory 512.MB

    publishDir "${params.output_dir}/Gathered_Results/Database/"

    container "132205776083.dkr.ecr.eu-west-2.amazonaws.com/ghcr/awgl/tso500_post_processing:dev"

    input:
    tuple val(worksheet), stdin

    output:
    path "samples_database_${worksheet}_RNA.csv"

    script:
    """
    cat > samples_database_${worksheet}_RNA.csv 
    """
}

process WRITE_SAMPLE_DB_LINE {
    cpus 1
    memory 512.MB

    container "community.wave.seqera.io/library/samtools:1.23--12d9384dd0649f36"

    input:
    val run_id
    tuple val(worksheet), val(sample_id), val(referral), path(sample_qc_file), val(ntc_sample_id), val(ntc_reads)

    output:
    tuple val(worksheet), stdout, emit: sample_db_entry

    script:
    """
    sample_reads=\$(tail -n1 ${sample_qc_file} | cut -f5)
    echo ${sample_id},${worksheet},TSO500_RNA,${referral},${run_id},GRCh37,\${sample_reads},${ntc_reads}
    """
}
