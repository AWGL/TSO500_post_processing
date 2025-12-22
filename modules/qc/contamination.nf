process CONTAMINATION {
    cpus 1
    memory 512.MB

    container "132205776083.dkr.ecr.eu-west-2.amazonaws.com/ghcr/awgl/tso500_post_processing:dev"
    publishDir "${params.output_dir}"

    input:
    path samples_correct_order
    tuple val(sample_id), val(worksheet), path(fusion_checks, name: "Gathered_Results/Database/")
    path referrals

    output:
    path "contamination-${worksheet}.csv"

    script:
    """
    cat ${samples_correct_order}
    # Move referrals to where the script expects them
    mkdir -p /data/diagnostics/pipelines/TSO500/TSO500_post_processing-main/
    mv RNA_referrals /data/diagnostics/pipelines/TSO500/TSO500_post_processing-main/RNA_referrals

    contamination_TSO500.py ${worksheet} main
    """
}
