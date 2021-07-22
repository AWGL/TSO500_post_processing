#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --output=%j-%N-1_TSO500.output
#SBATCH --error=%j-%N-1_TSO500.error
#SBATCH --partition=demultiplexing
#SBATCH --cpus-per-task=20

# Description: Demultiplex run using Illumina TSO500 app
# Author:      AWMGS
# Mode:        BY_SAMPLE
# Use:         sbatch within /Output/fastq/run_id directory

version=2.2.0
pipeline_dir=/data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-"$version"


cd $SLURM_SUBMIT_DIR
mkdir Demultiplex_Output

module purge
module load singularity

# catch fails early and terminate
set -euo pipefail

ln -s /data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-2.2.0/trusight-oncology-500-ruo.img .

# make sure to use singularity flag
$pipeline_dir/TruSight_Oncology_500_RUO.sh \
  --resourcesFolder $pipeline_dir/resources \
  --analysisFolder $SLURM_SUBMIT_DIR/Demultiplex_Output \
  --runFolder $raw_data \
  --engine singularity \
  --sampleSheet "$raw_data"/SampleSheet.csv \
  --isNovaSeq \
  --demultiplexOnly

# TODO - make variables files

# TODO - cd into /Output/results/run_id

# mkdir app_output
# cd app_output


# make a list of samples and get correct order of samples for each worksheet

cp "$raw_data"/SampleSheet.csv .

# remove header
sed -n -e '/Sample_ID,Sample_Name/,$p' SampleSheet.csv >> SampleSheet_updated.csv

# 
python "$pipeline_dir"/filter_sample_list.py

# make an empty file for recording completed samples 
> completed_samples.txt

# kick off script 2 for each sample
cat sample_list.txt | while read line; do


    # TODO - make sample folder

    echo kicking off pipeline for $line
    sbatch --export=raw_data="$raw_data",sample_id="$line" 2_TSO500.sh

done

