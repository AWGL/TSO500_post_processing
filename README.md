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



## Unit tests

Unit tests have been created against the following scripts: `tsv2db.py`, `coverage2json.py`, `fusion_check_with_ntc.py`.

To run all unit tests:
- copy these scripts into the `tests/` folder (relative imports not currently working, will be fixed in future version)
- activate the `TSO500_post_processing` conda environment
- change into tests directory: `cd tests/`
- run `python -m unittest`

To run tests on a specific script, follow the steps above but run `python -m unittest <test_script_name>`
