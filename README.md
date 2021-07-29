# TSO500_post_processing


##Documentation

The Illumina TSO500 manual can be found on Q-Pulse 

`EI-GEN-LocalAppUserGuide_TSO500`

## To run the pipeline

To run 1_TSO500.sh:

`sbatch --export=raw_data=/data/archive/novaseq/<run_id>/ 1_TSO500.sh` 

##Samplesheet requirements
* The samplesheet must contain the samples in the correct order for the RNA contamination results to be valid
* Every NTC must be named NTC-<worksheetid> 
