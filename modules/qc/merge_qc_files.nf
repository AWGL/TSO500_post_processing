process MERGE_QC_FILES {
    tag "${worksheet}"

    container "132205776083.dkr.ecr.eu-west-2.amazonaws.com/ghcr/awgl/tso500_post_processing:dev"

    publishDir "${params.output_dir}/"

    input:
    tuple val(sample_id), val(worksheet), val(referral), path(qc_files)

    output:
    path "RNA_QC_combined.txt"

    script:
    """
    for file in ${qc_files}; do
        if [[ ! -f RNA_QC_combined.txt ]]; then
            cat \$file > RNA_QC_combined.txt
        else
            cat \$file | tail -n1 >> RNA_QC_combined.txt
        fi
    done
    """
}

