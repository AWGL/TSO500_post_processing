#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --output=%j-%N-2_TSO500.output
#SBATCH --error=%j-%N-2_TSO500.error
#SBATCH --partition=medium
#SBATCH --cpus-per-task=20

# Description: Run Illumina TSO500 app
# Author:      AWMGS
# Mode:        BY_SAMPLE
# Use:         sbatch within run directory, pass in raw_data directory and sample_id

version=2.2.0

SCRATCH_DIR=/localscratch/"$SLURM_JOB_ID"
mkdir -p "$SCRATCH_DIR"
cd "$SCRATCH_DIR"
mkdir $sample_id

module purge
module load singularity

# catch fails early and terminate
set -euo pipefail

pipeline_dir=/data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-"$version"

ln -s /data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-2.2.0/trusight-oncology-500-ruo.img .

# use sampleOrPairIDs flag to run one sample at a time
$pipeline_dir/TruSight_Oncology_500_RUO.sh \
  --analysisFolder $sample_id \
  --resourcesFolder $pipeline_dir/resources \
  --fastqFolder $SLURM_SUBMIT_DIR/Demultiplex_Output \
  --isNovaSeq \
  --sampleSheet "$raw_data"SampleSheet.csv \
  --engine singularity \
  --sampleOrPairIDs $sample_id

cp -r $sample_id $SLURM_SUBMIT_DIR/$sample_id

##############################################################################################
# Run level
##############################################################################################


cd $SLURM_SUBMIT_DIR

# add to completed list once finished
echo $sample_id >> completed_samples.txt

# only run once all samples have finished
expected=$(cat sample_list.txt| wc -l)
complete=$(cat completed_samples.txt | wc -l)


if [ "$complete" -eq "$expected" ]; then

    touch entered_loop.txt
    sbatch --export=raw_data="$raw_data" 3_TSO500.sh

fi

# remove scratch dir
rm -r $SCRATCH_DIR
