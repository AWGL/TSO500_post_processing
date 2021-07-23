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

cd "$SLURM_SUBMIT_DIR"
mkdir $sample_id
#cd $sample_id

module purge
module load singularity
. ~/.bashrc
module load anaconda

# catch fails early and terminate
set -euo pipefail

pipeline_dir=/data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-"$version"


# use sampleOrPairIDs flag to run one sample at a time
$pipeline_dir/TruSight_Oncology_500_RUO.sh \
  --analysisFolder "$SLURM_SUBMIT_DIR"/"$sample_id" \
  --resourcesFolder $pipeline_dir/resources \
  --fastqFolder "$SLURM_SUBMIT_DIR"/Demultiplex_Output/Logs_Intermediates/FastqGeneration \
  --isNovaSeq \
  --sampleSheet "$raw_data"/SampleSheet.csv \
  --engine singularity \
  --sampleOrPairIDs $sample_id


#run fastqc


conda activate TSO500_post_processing

fastqc_path=/Output/results/"$run_id"/"$sample_id"/FastQC/

fastq_path=./Demultiplex_Output/Logs_Intermediates/FastqGeneration/"$sample_id"/
cd "$fastq_path"

for fastqPair in $(ls "$sample_id"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do

    #parse fastq filenames
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
    read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)


    fastqc --extract -o "$fastqc_path" "$read1Fastq"
    fastqc --extract -o "$fastqc_path" "$read2Fastq"


    mv "$fastqc_path"/"$sample_id"_S*_"$laneId"_R1_001_fastqc/summary.txt "$fastqc_path"/"$sample_id"_"$laneId"_R1_fastqc.txt
    mv "$fastqc_path"/"$sample_id"_S*_"$laneId"_R2_001_fastqc/summary.txt "$fastqc_path"/"$sample_id"_"$laneId"_R2_fastqc.txt

done

cd "$SLURM_SUBMIT_DIR"

conda deactivate


bash "$pipeline_dir"/depthofcoverage.sh $sample_id $version


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
