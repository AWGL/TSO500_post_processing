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
mkdir Gathered_Results

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


#Create variants and coverage tables in correct format to import to database
for worksheetid in $(cat worksheets_dna.txt); do
	for line in $(cat samples_correct_order_"$worksheetid"_DNA.csv); do
		sample="$(echo "$line" | cut -d, -f1)"
		worksheetid=$(echo "$line" | cut -d, -f2)
		referral=$(echo "$line" | cut -d, -f4)
		python "$pipeline_dir"/tsv2db.py --tsvfile ./Gathered_Results/Results/"$sample"/"$sample"_CombinedVariantOutput.tsv --ntcfile ./Gathered_Results/Results/NTC-"$worksheetid"/NTC-"$worksheetid"_CombinedVariantOutput.tsv --outfile ./Gathered_Results/Results/"$sample"/tsv2db_output.tsv
		python "$pipeline_dir"/coverage2json.py --referral "$referral" --groups_folder "$pipeline_dir"/hotspot_coverage --sample_coverage "$sample"/Coverage_results/ --ntc_coverage NTC-"$worksheetid"/Coverage_results/

	done
done


#create fusions table in correct format to import to database
for worksheetid in $(cat worksheets_rna.txt); do
	for line in $(cat samples_correct_order_"$worksheetid"_RNA.csv); do
		sample="$(echo "$line" | cut -d, -f1)"
		worksheetid=$(echo "$line" | cut -d, -f2)
		python "$pipeline_dir"/fusions_check_with_ntc.py ./Gathered_Results/Results/"$sample"/"$sample"_CombinedVariantOutput.tsv ./Gathered_Results/Results/NTC-"$worksheetid"/NTC-"$worksheetid"_CombinedVariantOutput.tsv ./Gathered_Results/Results/"$sample"/
	done
done

#compileQC report
bash compileQcReport.sh

#Run contamination script
#for worksheetid in $(cat worksheets_rna.txt); do
#	python contamination_TSO500.py "$run_id" "$worksheetid" "$version"
#done


