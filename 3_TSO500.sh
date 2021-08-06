#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --output=%j-%N-3_TSO500.out
#SBATCH --error=%j-%N-3_TSO500.err
#SBATCH --partition=high
#SBATCH --cpus-per-task=20

# Description: ----
# Author:      AWMGS
# Mode:        BY_SAMPLE
# Use:         sbatch within run directory
# Version:     1.0.0


app_version=2.2.0
app_dir=/data/diagnostics/pipelines/TSO500/illumina_app/TSO500_RUO_LocalApp-"$app_version"

pipeline_version=master
pipeline_dir=/data/diagnostics/pipelines/TSO500/TSO500_post_processing-"$pipeline_version"


cd "$SLURM_SUBMIT_DIR"
mkdir Gathered_Results

module purge
module load singularity

# catch fails early and terminate
set -euo pipefail



##############################################################################################
#  Illumina app
##############################################################################################

samples_to_gather=$(python $pipeline_dir/gather_list.py $SLURM_SUBMIT_DIR/sample_list.txt $SLURM_SUBMIT_DIR)

# make sure to use singularity flag
$app_dir/TruSight_Oncology_500_RUO.sh \
  --analysisFolder Gathered_Results \
  --resourcesFolder $app_dir/resources \
  --runFolder $raw_data \
  --isNovaSeq \
  --sampleSheet "$raw_data"/SampleSheet.csv \
  --engine singularity \
  --gather $samples_to_gather



##############################################################################################
#  Generate inputs for variants database - DNA
##############################################################################################

mkdir Gathered_Results/Database

#Create variants and coverage tables in correct format to import to database
for worksheet_id in $(cat worksheets_dna.txt); do

	for line in $(cat samples_correct_order_"$worksheet_id"_DNA.csv); do
		sample="$(echo "$line" | cut -d, -f1)"
		worksheet_id=$(echo "$line" | cut -d, -f2)
		referral=$(echo "$line" | cut -d, -f4)

		python "$pipeline_dir"/tsv2db.py \
                  --tsvfile ./Gathered_Results/Results/"$sample"/"$sample"_CombinedVariantOutput.tsv \
                  --ntcfile ./Gathered_Results/Results/NTC-"$worksheet_id"/NTC-"$worksheet_id"_CombinedVariantOutput.tsv \
                  --outfile ./Gathered_Results/Database/"$sample"_variants.tsv


                echo $line
                if [[ "$referral" != "null" ]]; then
		    python "$pipeline_dir"/coverage2json.py \
                      --referral "$referral" \
                      --groups_folder "$pipeline_dir"/hotspot_coverage/ \
                      --sample_coverage analysis/"$sample"/depth_of_coverage/ \
                      --ntc_coverage analysis/NTC-"$worksheet_id"/depth_of_coverage/ \
                      --outfile ./Gathered_Results/Database/"$sample"_"$referral"_coverage.json
                fi
	done
done


##############################################################################################
#  Generate inputs for variants database - RNA
#  filter fusions by referral type for contamination script
##############################################################################################



#create fusions table in correct format to import to database
for worksheet_id in $(cat worksheets_rna.txt); do

	for line in $(cat samples_correct_order_"$worksheet_id"_RNA.csv); do
		sample="$(echo "$line" | cut -d, -f1)"
		worksheet_id=$(echo "$line" | cut -d, -f2)

                if [[ ! -f ./Gathered_Results/Results/"$sample"/"$sample"_AllFusions.csv ]]; then
                    echo "fusion,exons,reference_reads_1,reference_reads_2,fusion_supporting_reads,left_breakpoint,right_breakpoint,type,in_ntc,spanning_reads,spanning_reads_dedup,split_reads,split_reads_dedup,fusion_caller,fusion_score"> ./Gathered_Results/Database/"$sample"_fusion_check.csv
 
                else

                    python "$pipeline_dir"/fusions_check_with_ntc.py \
                      ./Gathered_Results/Results/"$sample"/"$sample"_CombinedVariantOutput.tsv \
                      ./Gathered_Results/Results/NTC-"$worksheet_id"/NTC-"$worksheet_id"_CombinedVariantOutput.tsv \
                      ./Gathered_Results/Results/"$sample"/"$sample"_AllFusions.csv \
                      ./Gathered_Results/Database/

                fi
                
	done
done



# TODO - move BAMs into gathered results


##############################################################################################
#  QC
##############################################################################################

# compileQC report
# TODO concatenate all *_QC.txt files for DNA and RNA

#Run contamination script
for worksheetid in $(cat worksheets_rna.txt); do
	python "$pipeline_dir"/contamination_TSO500.py "$worksheetid" "$pipeline_version"
done

# move sample log files into their own folders
for sample in $(cat sample_list.txt)
do
    mv "$sample"_2_TSO500*.out analysis/"$sample"
    mv "$sample"_2_TSO500*.err analysis/"$sample"
done



# add timings
now=$(date +"%T")
echo "End time:   $now" >> timings.txt
