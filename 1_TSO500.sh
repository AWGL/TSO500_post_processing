#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --output=Demultiplex_Output-%j-%N.out
#SBATCH --error=Demultiplex_Output-%j-%N.err
#SBATCH --partition=demultiplexing
#SBATCH --cpus-per-task=24

# Description: Demultiplex run using Illumina TSO500 app
# Author:      AWMGS
# Mode:        BY_SAMPLE
# Use:         sbatch within /Output/fastq/run_id directory


app_version=2.2.0
app_dir=/data/diagnostics/pipelines/TSO500/illumina_app/TSO500_RUO_LocalApp-"$app_version"

pipeline_version=master
pipeline_dir=/data/diagnostics/pipelines/TSO500/TSO500_post_processing-"$pipeline_version"


cd $SLURM_SUBMIT_DIR
mkdir Demultiplex_Output
mkdir analysis

module purge
module load singularity

# catch fails early and terminate
set -euo pipefail

ln -s "$app_dir"/trusight-oncology-500-ruo.img .

now=$(date +"%T")
echo "Start time: $now" > timings.txt

##############################################################################################
#  Illumina app
##############################################################################################

# make sure to use singularity flag
"$app_dir"/TruSight_Oncology_500_RUO.sh \
  --resourcesFolder "$app_dir"/resources \
  --analysisFolder "$SLURM_SUBMIT_DIR"/Demultiplex_Output \
  --runFolder "$raw_data" \
  --engine singularity \
  --sampleSheet "$raw_data"/SampleSheet.csv \
  --isNovaSeq \
  --demultiplexOnly


##############################################################################################
#  Make variant lists for use downstream
##############################################################################################


# make a list of samples and get correct order of samples for each worksheet

cp "$raw_data"/SampleSheet.csv .

# remove header
sed -n -e '/Sample_ID,Sample_Name/,$p' SampleSheet.csv >> SampleSheet_updated.csv

# 
python "$pipeline_dir"/filter_sample_list.py

# make an empty file for recording completed samples 
> completed_samples.txt


##############################################################################################
#  Kick off script 2
##############################################################################################

# kick off script 2 for each sample
cat sample_list.txt | while read line; do

    echo kicking off pipeline for $line
    sbatch \
      --export=raw_data="$raw_data",sample_id="$line" \
      --output="$line"_2_TSO500-%j-%N.out \
      --error="$line"_2_TSO500-%j-%N.err \
      2_TSO500.sh

done

