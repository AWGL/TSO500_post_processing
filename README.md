# TSO500_post_processing


## Documentation

The Illumina TSO500 manual can be found on Q-Pulse 

`EI-GEN-LocalAppUserGuide_TSO500`

## To run the pipeline

To start the pipeline from demultiplexing, 1_TSO500.sh script should be copied into the run folder (/data/output/results/runid/).
From this folder run the command:

`sbatch --export=raw_data=/data/archive/novaseq/<run_id>/ 1_TSO500.sh` 

The raw data directory must contain the SampleSheet.csv. 

## Duty scientist responsibilites
The duty scientist is responsible for the following tasks:
* Creating the TSO500 samplesheet- see below for requirements
* Signing off the run in autoqc database - for more information refer to the AutoQC sop
* Importing the data into the data into the somatic variant database


## Samplesheet requirements
* The samplesheet must contain the samples in the correct order for the RNA contamination results to be valid
* Every NTC must be named NTC-worksheetid
