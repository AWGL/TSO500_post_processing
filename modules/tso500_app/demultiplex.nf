process APP_DEMULTIPLEX {
    memory 16.GB
    cpus 8
    
    container "132205776083.dkr.ecr.eu-west-2.amazonaws.com/tso500_local_app:ruo-2.2.0.12"

    containerOptions "--entrypoint /usr/bin/env"

    input:
    path(sample_sheet, name: "samplesheet/SampleSheet.csv")
    path(resources, name: "resources")
    path(run_folder, name: "run-folder")

    output:
    path('*.fastq.gz')

    script:
    """
    # Move stuff to /
    ln -s `readlink resources` /opt/illumina/resources
    ln -s `readlink run-folder` /opt/illumina/run-folder
    ln -s `readlink samplesheet/SampleSheet.csv` /opt/illumina/SampleSheet.csv

    mkdir analysis-folder
    rm -rf /opt/illumina/analysis-folder
    ln -sf \$PWD/analysis-folder /opt/illumina/analysis-folder 

    
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

    java -jar /opt/cromwell/cromwell-36.jar run -i inputs.json /opt/illumina/wdl/TSO500Workflow.wdl
    """    
}