process APP_PIPELINE {
    memory 64.GB
    cpus 16

    container "132205776083.dkr.ecr.eu-west-2.amazonaws.com/tso500_local_app_custom_entrypoint:ruo-2.2.0.12"

    publishDir "${params.output_dir}/analysis/${sample_id}"

    input:
    path sample_sheet, name: "samplesheet/SampleSheet.csv"
    path resources, name: "resources"
    tuple val(sample_id), val(worksheet), val(type), val(referral), path(fastqs, name: "fastq-folder/*")

    output:
    path "Logs_Intermediates/"
    tuple val(sample_id), val(worksheet), val(referral), path("Results/"), emit: results
    tuple val(sample_id), val(worksheet), val(referral), path("Results/MetricsOutput.tsv"), emit: metrics_output
    tuple val(sample_id), val(worksheet), path("Logs_Intermediates/RnaMarkDuplicates/${sample_id}/${sample_id}.bam"), emit: bams

    script:
    """
    echo "Linking resources"
    ln -s `realpath resources` /opt/illumina/resources

    echo "Linking run folder"
    ln -s `realpath run-folder` /opt/illumina/run-folder

    echo "Putting fastqs in subdir"
    mkdir fastq-folder/${sample_id}
    mv fastq-folder/*.fastq.gz fastq-folder/${sample_id}

    echo "Linking fastq folder"
    ln -s `realpath fastq-folder` /opt/illumina/fastq-folder

    echo "Linking SampleSheet"
    ln -s `realpath samplesheet/SampleSheet.csv` /opt/illumina/SampleSheet.csv
    
    echo "Creating inputs.json"
    # Generate input.json 
    echo '{' >> inputs.json
    echo '    "TSO500.workflowVersion": "ruo-2.2.0.12",' >> inputs.json
    echo '    "TSO500.isNovaSeq": true,' >> inputs.json
    echo '    "TSO500.runFolder": "/opt/illumina/run-folder",' >> inputs.json
    echo '    "TSO500.startFromFastq": true,' >> inputs.json
    echo '    "TSO500.fastqFolder": "/opt/illumina/fastq-folder",' >> inputs.json
    echo '    "TSO500.lrmMode": false,' >> inputs.json
    echo '    "TSO500.sampleOrPairIDs": "${sample_id}",' >> inputs.json
    echo '    "TSO500.demultiplex": false' >> inputs.json
    echo '}' >> inputs.json

    echo "Running cromwell"
    java -jar /opt/cromwell/cromwell-36.jar run -i inputs.json /opt/illumina/wdl/TSO500Workflow.wdl

    mv /opt/illumina/analysis-folder/* .
    """
}
