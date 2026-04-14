# TSO500_post_processing

This pipeline
* Demultiplexes with the Illumina app
* Runs the RNA samples through the Illumina app
* Submits the DNA samples to run through the [somatic_enrichment_nextflow](https://github.com/AWGL/somatic_enrichment_nextflow) pipeline

## Documentation

The Illumina TSO500 manual can be found on Q-Pulse 

`EI-GEN-LocalAppUserGuide_TSO500`

## To run the pipeline

To start the pipeline from demultiplexing, 1_TSO500.sh script should be copied into the run folder (/data/output/results/runid/TSO500/).
From this folder run the command:

`sbatch --export=raw_data=/data/raw/novaseq/<run_id>/ 1_TSO500.sh` 

The raw data directory must contain the SampleSheet.csv. 

## To run a DNA sample

To run a DNA sample through the Illumina app, please run an old release (v1.1.1)

## Duty scientist responsibilites
The duty scientist is responsible for the following tasks:
* Signing off the run in autoqc database and GLIMS (NOTE: the RNA must now be signed off with the DNA)
* Signing off the run in SVD

## Samplesheet requirements
* The samplesheet must contain the samples in the correct order for the RNA contamination results to be valid
* Every NTC must be named NTC-worksheetid


## Unit tests

Unit tests have been created against the following scripts: `fusion_check_with_ntc.py`.

To run all unit tests:
- copy these scripts into the `tests/` folder (relative imports not currently working, will be fixed in future version)
- activate the `TSO500_post_processing` conda environment
- change into tests directory: `cd tests/`
- run `python -m unittest`

To run tests on a specific script, follow the steps above but run `python -m unittest <test_script_name>`


## Adding a new RNA panel

To add a new panel, the following needs to be changed:

**Samplesheet generator:**
- Add the referral to the SampleSheet Generator, including the mapped test directory code

**Pipeline:**
- Add a new file to RNA_referrals named <panel>.txt with the gene names on the panel, one per line

**Somatic variant database:**
- Make a new Panel object in SVD
