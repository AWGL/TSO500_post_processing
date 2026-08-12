process APP_DEMULTIPLEX {
    memory 64.GB
    cpus 16

    label 'app_container'
    label 'demultiplex_queue'
    
    publishDir "${params.output_dir}/Demultiplex_Output"

    input:
    path(sample_sheet, name: "samplesheet/SampleSheet.csv")
    path(resources, name: "resources")
    path(run_folder, name: "run-folder")

    output:
    path('Logs_Intermediates'), emit: demultiplex_output
    path('Logs_Intermediates/**/*.fastq.gz'), emit: fastq

    script:
    """
    echo "Linking resources"
    ln -s `realpath resources` /opt/illumina/resources

    echo "Linking run folder"
    ln -s `realpath run-folder` /opt/illumina/run-folder

    echo "Linking SampleSheet"
    ln -s `realpath samplesheet/SampleSheet.csv` /opt/illumina/SampleSheet.csv
    
    echo "Creating inputs.json"
    # Generate input.json 
    echo '{' >> inputs.json
    echo '    "TSO500.workflowVersion": "ruo-2.2.0.12",' >> inputs.json
    echo '    "TSO500.isNovaSeq": true,' >> inputs.json
    echo '    "TSO500.runFolder": "/opt/illumina/run-folder",' >> inputs.json
    echo '    "TSO500.startFromFastq": false,' >> inputs.json
    echo '    "TSO500.isSingleLaneMode": true,' >> inputs.json
    echo '    "TSO500.lrmMode": false,' >> inputs.json
    echo '    "TSO500.demultiplex": true' >> inputs.json
    echo '}' >> inputs.json

    echo "Running cromwell"
    java -jar /opt/cromwell/cromwell-36.jar run -i inputs.json /opt/illumina/wdl/TSO500Workflow.wdl
    mv /opt/illumina/analysis-folder/Logs_Intermediates .
    """    
}