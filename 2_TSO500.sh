#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --output=%j-%N-2_TSO500.out
#SBATCH --error=%j-%N-2_TSO500.err
#SBATCH --partition=high
#SBATCH --cpus-per-task=24

# Description: Run Illumina TSO500 app for each sample then run postprocessing steps for RNA only. 
#              Kick off script 3 when all samples completed
# Use:         from /Output/results/<run_id>/TSO500/ directory, for each sample run: 
#              sbatch --export=raw_data=/data/raw/novaseq/<run_id>,sample_id=<sample_id> 2_TSO500.sh
# Version:     1.0.15

##############################################################################################
#  Setup
##############################################################################################

# define filepaths for app
app_version=2.2.0
app_dir=/data/diagnostics/pipelines/TSO500/illumina_app/TSO500_RUO_LocalApp-"$app_version"

# define filepaths for post processing
pipeline_version=main
pipeline_dir=/data/diagnostics/pipelines/TSO500/TSO500_post_processing-"$pipeline_version"
pipeline_scripts="$pipeline_dir"/scripts

# setup analysis folders
cd "$SLURM_SUBMIT_DIR"
output_path="$SLURM_SUBMIT_DIR"/analysis/"$sample_id"
mkdir -p $output_path

# load singularity and anaconda modules
module purge
module load singularity
. ~/.bashrc
module load anaconda

# catch fails early and terminate
set -euo pipefail

##############################################################################################
#  Illumina app
##############################################################################################

# use sampleOrPairIDs flag to run one sample at a time
"$app_dir"/TruSight_Oncology_500_RUO.sh \
  --analysisFolder "$output_path" \
  --resourcesFolder "$app_dir"/resources \
  --fastqFolder "$SLURM_SUBMIT_DIR"/Demultiplex_Output/Logs_Intermediates/FastqGeneration \
  --isNovaSeq \
  --sampleSheet "$raw_data"/SampleSheet.csv \
  --engine singularity \
  --sampleOrPairIDs "$sample_id"


##############################################################################################
#  FastQC
##############################################################################################

# activate conda env
set +u
conda activate TSO500_post_processing
set -u

# make fastqc output folder in the sample folder
fastqc_output="$output_path"/FastQC/
mkdir -p $fastqc_output

# location of sample fastqs
fastq_path="$SLURM_SUBMIT_DIR"/Demultiplex_Output/Logs_Intermediates/FastqGeneration/"$sample_id"/

# run FastQC for each fastq pair
for fastqPair in $(find $fastq_path -name *.fastq.gz -type f -printf "%f\n" | cut -d_ -f1-3 | sort | uniq); do

    #parse fastq filenames
    laneId=$(echo basename "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastq_path""$fastqPair"_R1_*fastq.gz)
    read2Fastq=$(ls "$fastq_path""$fastqPair"_R2_*fastq.gz)

    # run FastQC
    fastqc --extract -o "$fastqc_output" "$read1Fastq"
    fastqc --extract -o "$fastqc_output" "$read2Fastq"

    # rename files
    mv "$fastqc_output"/"$sample_id"_S*_"$laneId"_R1_001_fastqc/summary.txt "$fastqc_output"/"$sample_id"_"$laneId"_R1_fastqc.txt
    mv "$fastqc_output"/"$sample_id"_S*_"$laneId"_R2_001_fastqc/summary.txt "$fastqc_output"/"$sample_id"_"$laneId"_R2_fastqc.txt

done


##############################################################################################
#  Gather QC metrics for sample
##############################################################################################

# function to check FASTQC output
count_qc_fails() {
    #count how many core FASTQC tests failed
    grep -E "Basic Statistics|Per base sequence quality|Per sequence quality scores|Per base N content" "$1" | \
    grep -v ^PASS | \
    grep -v ^WARN | \
    wc -l | \
    sed 's/^[[:space:]]*//g'
}

# check all fastqc files for any fails
fastqc_status=PASS

for report in "$fastqc_output"/"$sample_id"_*_fastqc.txt; do
    if [ $(count_qc_fails $report) -gt 0 ]; then
        fastqc_status=FAIL
    fi
done

# pull out metrics from the Illumina app MetricsOutput.tsv
completed_all_steps=$(grep COMPLETED_ALL_STEPS analysis/"$sample_id"/Results/MetricsOutput.tsv | cut -f2)

# RNA only metrics
median_cv_gene_500x=$(grep "MEDIAN_CV_GENE_500X" analysis/"$sample_id"/Results/MetricsOutput.tsv | cut -f4)
total_on_target_reads=$(grep "TOTAL_ON_TARGET_READS" analysis/"$sample_id"/Results/MetricsOutput.tsv | tail -n1 | cut -f4)
median_insert_size=$(grep "MEDIAN_INSERT_SIZE" analysis/"$sample_id"/Results/MetricsOutput.tsv | tail -n1 | cut -f4)
total_pf_reads=$(grep "TOTAL_PF_READS" analysis/"$sample_id"/Results/MetricsOutput.tsv | tail -n1 | cut -f4)

# add to sample QC file
echo -e "Sample\tFastQC\tcompleted_all_steps\tmedian_cv_gene_500x\ttotal_on_target_reads\tmedian_insert_size\ttotal_pf_reads" > "$output_path"/"$sample_id"_RNA_QC.txt
echo -e "$sample_id\t$fastqc_status\t$completed_all_steps\t$median_cv_gene_500x\t$total_on_target_reads\t$median_insert_size\t$total_pf_reads" >> "$output_path"/"$sample_id"_RNA_QC.txt


##############################################################################################
#  Kick off run level script once all samples have finished
##############################################################################################

# add sample to completed list once finished
echo $sample_id >> completed_samples.txt

# only run once all samples have finished
expected=$(cat samples_correct_order_*_RNA.csv | wc -l)
complete=$(cat completed_samples.txt | wc -l)

# if last sample, kick off script 3
if [ "$complete" -eq "$expected" ]; then
    sbatch --export=raw_data="$raw_data" 3_TSO500.sh
fi

# deactivate env
set +u
conda deactivate
set -u
