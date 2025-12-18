process CREATE_SAMPLE_QC_FILE {
    tag "${sample_id}"

    container "132205776083.dkr.ecr.eu-west-2.amazonaws.com/ghcr/awgl/tso500_post_processing:dev"

    publishDir "${params.output_dir}/analysis/${sample_id}/"

    input:
    tuple val(sample_id), path(fastqc_summary), val(worksheet), path(metrics_file)

    output:
    tuple val(sample_id), val(worksheet), path("${sample_id}_RNA_QC.txt"), emit: sample_qc_file

    script:
    """
    count_qc_fails() {
        #count how many core FASTQC tests failed
        grep -E "Basic Statistics|Per base sequence quality|Per sequence quality scores|Per base N content" "\$1" | \
        grep -v ^PASS | \
        grep -v ^WARN | \
        wc -l | \
        sed 's/^[[:space:]]*//g'
    }

    # check all fastqc files for any fails
    fastqc_status=PASS

    for report in ${fastqc_summary}; do
    if [ \$(count_qc_fails \$report) -gt 0 ]; then
        fastqc_status=FAIL
    fi
    done

    completed_all_steps=\$(grep COMPLETED_ALL_STEPS ${metrics_file} | cut -f2)

    median_cv_gene_500x=\$(grep "MEDIAN_CV_GENE_500X" ${metrics_file} | cut -f4)
    total_on_target_reads=\$(grep "TOTAL_ON_TARGET_READS" ${metrics_file} | tail -n1 | cut -f4)
    median_insert_size=\$(grep "MEDIAN_INSERT_SIZE" ${metrics_file} | tail -n1 | cut -f4)
    total_pf_reads=\$(grep "TOTAL_PF_READS" ${metrics_file}| tail -n1 | cut -f4)

    # add to sample QC file
    echo -e "Sample\tFastQC\tcompleted_all_steps\tmedian_cv_gene_500x\ttotal_on_target_reads\tmedian_insert_size\ttotal_pf_reads" > ${sample_id}_RNA_QC.txt
    echo -e "${sample_id} \t\$fastqc_status\t\$completed_all_steps\t\$median_cv_gene_500x\t\$total_on_target_reads\t\$median_insert_size\t\$total_pf_reads" >> ${sample_id}_RNA_QC.txt
    """
}