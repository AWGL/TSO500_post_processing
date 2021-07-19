#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --output=%j-%N-3_TSO500.output
#SBATCH --error=%j-%N-3_TSO500.error
#SBATCH --partition=medium
#SBATCH --cpus-per-task=20

# Description: ----
# Author:      AWMGS
# Mode:        BY_SAMPLE
# Use:         sbatch within run directory

version=2.2.0

#cd "$SLURM_SUBMIT_DIR"

SCRATCH_DIR=/localscratch/"$SLURM_JOB_ID"
mkdir -p "$SCRATCH_DIR"
cd "$SCRATCH_DIR"
mkdir Gathered_Results

module purge
module load singularity

# catch fails early and terminate
set -euo pipefail

pipeline_dir=/data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-"$version"

# run with gather flag
#mkdir Gathered_Results

samples_to_gather=$(python $pipeline_dir/gather_list.py $SLURM_SUBMIT_DIR/sample_list.txt $SLURM_SUBMIT_DIR)

ln -s /data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-2.2.0/trusight-oncology-500-ruo.img .

# make sure to use singularity flag
$pipeline_dir/TruSight_Oncology_500_RUO.sh \
  --analysisFolder Gathered_Results \
  --resourcesFolder $pipeline_dir/resources \
  --runFolder $raw_data \
  --isNovaSeq \
  --sampleSheet "$raw_data"SampleSheet.csv \
  --engine singularity \
  --gather $samples_to_gather

cp -r Gathered_Results/Results $SLURM_SUBMIR_DIR/Gathered_Results

# run custom analysis
cd $SLURM_SUBMIT_DIR
touch run_complete.txt


rm -r $SCRATCH_DIR
