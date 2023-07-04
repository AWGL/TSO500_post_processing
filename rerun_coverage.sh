#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --output=%j-%N-rerun_coverage.out
#SBATCH --error=%j-%N-rerun_coverage.err
#SBATCH --partition=high

#Description: Script to re-run coverage calculations for TSO500 App based pipeline when a sample has a new referral. #			Can be just cov2json or from coverage calculator if the analysis folder is deleted. 
#
#Usage: rerun_coverage.sh <RUN_ID> <SAMPLE_ID> <WORKSHEET_ID> <NEW_REFERRAL>

##############################################################
# Setup
##############################################################

#Get variables from command line
run_id=$1
sample=$2
worksheet_id=$3
referral=$4

output_path=/data/output/results/${run_id}

# define filepaths for app
app_version=2.2.0
app_dir=/data/diagnostics/pipelines/TSO500/illumina_app/TSO500_RUO_LocalApp-"$app_version"

# define filepaths for post processing
pipeline_version=master
pipeline_dir=/data/diagnostics/pipelines/TSO500/TSO500_post_processing-"$pipeline_version"
pipeline_scripts="$pipeline_dir"/scripts

# define pipeline variables
minimum_coverage="270 135"
coverage_bed_files_path="$pipeline_dir"/hotspot_coverage
vendor_capture_bed="$pipeline_dir"/vendorCaptureBed_100pad_updated.bed
preferred_transcripts="$pipeline_dir"/preferred_transcripts.txt

#Load anaconda and environment
. ~/.bashrc
module load anaconda

set +u
conda activate TSO500_post_processing
set -u

#Check if analysis folder exists. 
#If it does, check if referral was present at time of run. 
#If not referral at time, or no analysis, start from DOC, otherwise run cov2json.py

analysis_path=${output_path}/analysis

if [ -d "$analysis_path" ];
then
	if [ -f "${analysis_path}/NTC-${worksheet_id}/depth_of_coverage/hotspot_coverage_135x/NTC-${worksheet_id}_${referral}_combined.coverage" ];
	then
		analysis="True"
	else
		analysis="False"
	fi
else
	analysis="False"
fi

#Run cov2json if analysis present, otherwise run from DOC
if [ $analysis = "True" ];
then

	python "$pipeline_scripts"/coverage2json.py \
		--referral "$referral" \
		--groups_folder "$pipeline_dir"/hotspot_coverage/ \
		--sample_coverage "$analysis_path"/"$sample"/depth_of_coverage/ \
		--ntc_coverage analysis/NTC-"$worksheet_id"/depth_of_coverage/ \
		--outfile "$output_path"/Gathered_Results/Database/"$sample"_"$referral"_coverage.json	

else

	bam_path="${output_path}/BAMs/"

	#Do everything for sample and NTC
	for s_id in $sample NTC-${worksheet_id}
	do

		depth_path="${output_path}/reanalysis/${s_id}/depth_of_coverage/"
		mkdir ${output_path}/reanalysis
		mkdir ${output_path}/reanalysis/${s_id}
		mkdir ${output_path}/reanalysis/${s_id}/depth_of_coverage/

		# reheader the bams to local area
		java -jar /Apps/wren/picard/2.21.6/bin/picard.jar AddOrReplaceReadGroups \
			I="$bam_path"/"$s_id".bam \
			O="$bam_path"/"$s_id"_add_rg.bam \
			RGID=4 \
			RGLB=lib1 \
			RGPL=ILLUMINA \
			RGPU=unit1 \
			RGSM=20

		# index new bam
		samtools index "$bam_path"/"$s_id"_add_rg.bam "$bam_path"/"$s_id"_add_rg.bam.bai

		# run depth of coverage
		gatk DepthOfCoverage \
			-I "$bam_path"/"$s_id"_add_rg.bam \
			-L "$vendor_capture_bed" \
			-R "$app_dir"/resources/genomes/hg19_hardPAR/genome.fa \
			-O "$depth_path"/"$s_id"_depth_of_coverage

		# change to tab delimited and remove colon from column 1
		sed 's/:/\t/g' "$depth_path"/"$s_id"_depth_of_coverage \
			| sed 's/,/\t/g' | grep -v 'Locus' \
			| sort -k1,1 -k2,2n | bgzip \
			> "$depth_path"/"$s_id"_depth_of_coverage.gz

		# tabix index depth of coverage file
		tabix \
			-b 2 \
			-e 2 \
			-s 1 \
			"$depth_path"/"$s_id"_depth_of_coverage.gz

		# deactivate env
		set +u
		conda deactivate
		set -u

		# Run coverage calculator for each coverage value
		for min_coverage in $minimum_coverage; do

			# activate coverage calculator conda env
			set +u
			conda activate CoverageCalculatorPy
			set -u

			# set output directory for coverage files
			hscov_outdir=hotspot_coverage_"$min_coverage"x

			# run coverage calculator on each bed file
			for bed_file in "$coverage_bed_files_path"/*.bed; do

				name=$(echo $(basename $bed_file) | cut -d"." -f1)

				python /data/diagnostics/apps/CoverageCalculatorPy/CoverageCalculatorPy-v1.1.0/CoverageCalculatorPy.py \
					-B "$coverage_bed_files_path"/"$name".bed \
					-D "$depth_path"/"$s_id"_depth_of_coverage.gz \
					--depth "$min_coverage" \
					--padding 0 \
					--groupfile "$coverage_bed_files_path"/"$name".groups \
					--outname "$s_id"_"$name" \
					--outdir  "$depth_path"/"$hscov_outdir"/

				# remove header from gaps file
				if [[ $(wc -l < "$depth_path"/"$hscov_outdir"/"$s_id"_"$name".gaps) -eq 1 ]]; then
					# no gaps
					touch "$depth_path"/"$hscov_outdir"/"$s_id"_"$name".nohead.gaps

				else
					# gaps
					grep -v '^#' "$depth_path"/"$hscov_outdir"/"$s_id"_"$name".gaps > "$depth_path"/"$hscov_outdir"/"$s_id"_"$name".nohead.gaps

				fi

				# remove chr from bed file so bed2hgvs works
				cat "$depth_path"/"$hscov_outdir"/"$s_id"_"$name".nohead.gaps | sed 's/^chr//' > "$depth_path"/"$hscov_outdir"/"$s_id"_"$name".nohead_nochr.gaps

				# remove intermediate files
				rm "$depth_path"/"$hscov_outdir"/"$s_id"_"$name".gaps
				rm "$depth_path"/"$hscov_outdir"/"$s_id"_"$name".nohead.gaps

			done


        		# activate bed2hgvs conda env
			set +u
			conda deactivate
			conda activate bed2hgvs
			set -u

			# run on each bed file
			for gaps_file in "$depth_path"/"$hscov_outdir"/*.nohead_nochr.gaps; do

				name=$(echo $(basename $gaps_file) | cut -d"." -f1)
				echo $name

				Rscript /data/diagnostics/apps/bed2hgvs/bed2hgvs-v0.3.0/bed2hgvs.R \
					--bedfile $gaps_file \
					--outname "$name".gaps \
					--outdir "$depth_path"/"$hscov_outdir" \
					--preferred_tx $preferred_transcripts

				# remove intermediate file
				rm "$depth_path"/"$hscov_outdir"/"$name".nohead_nochr.gaps
			done

			# combine all total coverage files
			if [ -f "$depth_path"/"$hscov_outdir"/"$s_id"_coverage.txt ]; then rm "$depth_path"/"$hscov_outdir"/"$s_id"_coverage.txt; fi
			cat "$depth_path"/"$hscov_outdir"/*.totalCoverage | grep "FEATURE" | head -n 1 >> "$depth_path"/"$hscov_outdir"/"$s_id"_coverage.txt
			cat "$depth_path"/"$hscov_outdir"/*.totalCoverage | grep -v "FEATURE" | grep -vP "combined_\\S+_GENE" >> "$depth_path"/"$hscov_outdir"/"$s_id"_coverage.txt

			# deactivate env
			set +u
			conda deactivate
			set -u


			# activate conda env
			set +u
			conda activate TSO500_post_processing
			set -u

			cosmic_tool_path=/data/diagnostics/apps/cosmic_gaps/cosmic_gaps-master

			gaps_file="$depth_path"/"$hscov_outdir"/"$s_id"_"$referral"_hotspots.gaps

			# hotspot gaps file may be missing for some referrals
			if [[ -f $gaps_file ]]
			then

				# only run bedtools intersect for certain referral types
				if [ $referral = "Melanoma" ] ||  [ $referral = "Lung" ] || [ $referral = "Colorectal" ] || [ $referral = "GIST" ] || [ $referral = "breast" ]
				then
					dos2unix $gaps_file

					# find the overlap between the hotspots file and the referral file from cosmic
					bedtools intersect \
						-loj \
						-F 1 \
						-a $gaps_file \
						-b "$cosmic_tool_path"/cosmic_bedfiles/"$referral".bed \
						-wao \
						> "$depth_path"/"$hscov_outdir"/"$s_id"_"$referral"_intersect.txt

				fi

				# filter the output 
				python "$cosmic_tool_path"/filter_table.py \
					--sampleId $s_id \
					--referral $referral \
					--gaps_path "$depth_path"/"$hscov_outdir"/ \
					--bedfile_path "$cosmic_tool_path"/cosmic_bedfiles/
        
			fi
		done
	done

	python "$pipeline_scripts"/coverage2json.py \
		--referral "$referral" \
		--groups_folder "$pipeline_dir"/hotspot_coverage/ \
		--sample_coverage "$output_path"/reanalysis/"$sample"/depth_of_coverage/ \
		--ntc_coverage "$output_path"/reanalysis/NTC-"$worksheet_id"/depth_of_coverage/ \
		--outfile "$output_path"/Gathered_Results/Database/"$sample"_"$referral"_coverage.json
fi

# deactivate env
set +u
conda deactivate
set -u
