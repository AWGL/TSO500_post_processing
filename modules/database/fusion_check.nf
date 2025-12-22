process FUSION_CHECK {
    container "132205776083.dkr.ecr.eu-west-2.amazonaws.com/ghcr/awgl/tso500_post_processing:dev"

    publishDir "${params.output_dir}/Gathered_Results/Database/"

    input:
    tuple val(worksheet), val(sample_id), val(referral), path(results, name: "sample_results/*"), val(ntc_id), val(ntc_referral), path(ntc_results, name: "ntc_results/*")

    output:
    tuple val(sample_id), val(worksheet), path("${sample_id}_fusion_check.csv"), emit: fusion_check

    script:
    """
    ls sample_results/Results/${sample_id}/
    # AllFusions file isn't made if the app doesnt finish, make a blank one instead
    if [[ ! -f sample_results/Results/${sample_id}/${sample_id}_AllFusions.csv ]]; then
        echo "fusion,exons,reference_reads_1,reference_reads_2,fusion_supporting_reads,left_breakpoint,right_breakpoint,type,in_ntc,spanning_reads,spanning_reads_dedup,split_reads,split_reads_dedup,fusion_caller,fusion_score" > ${sample_id}_fusion_check.csv

    # if app completes properly, format fusions for database upload
    else
        echo "Running fusions2db!"
        fusions2db.py \
            --tsvfile sample_results/Results/${sample_id}/${sample_id}_CombinedVariantOutput.tsv \
            --ntcfile ntc_results/Results/${ntc_id}/${ntc_id}_CombinedVariantOutput.tsv \
            --allfusions sample_results/Results/${sample_id}/${sample_id}_AllFusions.csv \
            --outfile .
    fi
    """
}
