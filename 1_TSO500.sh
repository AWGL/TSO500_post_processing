#!/bin/bash

#SBATCH --output=Demultiplex_Output-%j-%N.out
#SBATCH --error=Demultiplex_Output-%j-%N.err
#SBATCH --partition=demultiplexing
#SBATCH --cpus-per-task=24

# Description: Demultiplex run using Illumina TSO500 app and kick off script 2 for each sample
# Use:         from /Output/results/<run_id> directory, run: 
#              sbatch --export=raw_data=/data/archive/novaseq/<run_id> 1_TSO500.sh
# Version:     1.0.9


##############################################################################################
#  Setup
##############################################################################################

# define filepaths for app
app_version=2.2.0
app_dir=/data/diagnostics/pipelines/TSO500/illumina_app/TSO500_RUO_LocalApp-"$app_version"

# define filepaths for post processing
pipeline_version=master
pipeline_dir=/data/diagnostics/pipelines/TSO500/TSO500_post_processing-"$pipeline_version"
pipeline_scripts="$pipeline_dir"/scripts

# setup analysis folders
cd $SLURM_SUBMIT_DIR
mkdir Demultiplex_Output
mkdir analysis

# load singularity and anaconda modules
module purge
module load singularity
. ~/.bashrc
module load anaconda

# catch fails early and terminate
set -euo pipefail

# soft link Illumina app singularity image to current directory
ln -s "$app_dir"/trusight-oncology-500-ruo.img .

# print pipeline start time
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
#  Postprocessing
##############################################################################################

# activate conda env
set +u
conda activate TSO500_post_processing
set -u

# copy samplesheet to current directory
cp "$raw_data"/SampleSheet.csv .

# remove header from samplesheet
sed -n -e '/Sample_ID,Sample_Name/,$p' SampleSheet.csv >> SampleSheet_updated.csv

# make a list of samples and get correct order of samples for each worksheet
python "$pipeline_scripts"/filter_sample_list.py

# create read counts bar chart - runs within it's own conda env
set +u
conda deactivate
conda activate read_count
set -u

python /data/diagnostics/scripts/read_count_visualisation.py

# reactivate main conda env
set +u
conda deactivate
conda activate TSO500_post_processing
set -u

# make an empty file for recording completed samples 
> completed_samples.txt

#run dos2unix on samplesheet ready for cosmic gaps section
dos2unix SampleSheet_updated.csv


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

# deactivate env
set +u
conda deactivate
set -u
