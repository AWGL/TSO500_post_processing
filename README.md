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


## Adding a new DNA panel

To add a new panel, the following needs to be changed:

**Samplesheet generator:**
- Add the referral reason to the samplesheet generator (see SOP in Qpulse) and make sure it matches the filename of the new bed files (case sensitive)

**Pipeline:**
- Generate bed files for the new panel
  - The `hotspot_variants/*bed` bed file and `hotspot_coverage/*combined.bed` files are required
  - The `hotspot_coverage/*hotspots.bed` and `hotspot_coverage/*genescreen.bed` files are optional
  - Filenames should be all lowercase
- Make sure that all regions in the new panel are covered in `vendorCaptureBed_100pad_updated.bed` (the bed file file used to generate the depth of coverage file)
- Make sure that any flanking regions are added to the `TSO_extra_padding_chr.interval_list` file - Illumina bed file only goes +/- 2bp so this file contains the extra 3bp to make it +/- 5bp

**Somatic variant database**
- Make a new panel object in the somatic variant database that matches the filename of the new bed files (case sensitive)
- Move the new variants bed file into the `roi/hotspot_variants` folder in the somatic database
