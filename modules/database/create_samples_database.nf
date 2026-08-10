process MERGE_SAMPLES_DATABASE {
    cpus 1
    memory 512.MB

    publishDir "${params.output_dir}/Gathered_Results/Database/"

    input:
    tuple val(worksheet), stdin

    output:
    path "samples_database_${worksheet}_RNA.csv"

    script:
    """
    echo "sample,worksheet,assay,referral,run,genome_build,sample_reads,ntc_reads" > samples_database_${worksheet}_RNA.csv
    cat >> samples_database_${worksheet}_RNA.csv
    """
}

process WRITE_SAMPLE_DB_LINE {
    cpus 1
    memory 512.MB

    label 'samtools_container'

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
