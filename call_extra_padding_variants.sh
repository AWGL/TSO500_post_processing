#!/bin/bash

# Illumina app only calls variants +/- 2bp of each exon, we require a larger area (usually +/- 5bp)
# Script re-runs part of illumina app to call variants between +/-2 to +/-5


##############################################################################################
#  Setup
##############################################################################################

# get sample ID from sys args
sample_id=$1

# define filepaths for app
app_version=2.2.0
app_dir=/data/diagnostics/pipelines/TSO500/illumina_app/TSO500_RUO_LocalApp-"$app_version"

# define filepaths for post processing
pipeline_version=master
pipeline_dir=/data/diagnostics/pipelines/TSO500/TSO500_post_processing-"$pipeline_version"

cd "$SLURM_SUBMIT_DIR"

# setup analysis folders and paths for singularity mounts
# folder for output of this script
output_folder="$SLURM_SUBMIT_DIR"/analysis/"$sample_id"/padding
mkdir -p $output_folder

# mount path where data from first pass of Illumina app is stored
app_logs_intermediates="$SLURM_SUBMIT_DIR"/analysis/"$sample_id"/Logs_Intermediates

# mount path for app resources folder
resources="$app_dir"/resources/

# load singularity and anaconda modules
module purge
module load singularity
. ~/.bashrc
module load anaconda

# catch fails early and terminate
set -euo pipefail

# define singularity exec command with folders mapped
SING="singularity exec --bind "$app_logs_intermediates":/mnt/app_result,"$output_folder":/opt/illumina/analysis-folder,"$resources":/opt/illumina/resources "$app_dir"/trusight-oncology-500-ruo.img"

# path to updated ROI - this is only visible inside singularity container
extra_roi_bed_file=/opt/illumina/resources/TSO_extra_padding_chr.interval_list


##############################################################################################
#  Re-run some app steps within singularity, output into output_folder variable above
#  Most setting exactly the same but bed file switched for expanded ROI
##############################################################################################

# re-run variant caller
$SING python /opt/illumina/step_modules/VariantCaller/variant_caller.py \
    --stepName VariantCaller \
    --dotnetPath dotnet \
    --outputFolder /opt/illumina/analysis-folder/Logs_Intermediates/VariantCaller \
    --dsdmFiles /mnt/app_result/Msi/dsdm.json \
    --inputBamFolder /mnt/app_result/StitchedRealigned \
    --inputMetricsFolder /mnt/app_result/CollapsedReads \
    --inputBamFileSuffix ".bam" \
    --outputVcfFileSuffix  ".genome.vcf" \
    --outputFilteredFileSuffix ".filtered.genome.vcf" \
    --piscesPath "/opt/illumina/components/pisces/Pisces.dll" \
    --additionalPiscesArgs "-G /opt/illumina/resources/genomes/hg19_hardPAR -I "$extra_roi_bed_file" -MinVF 0.0001 -SSFilter False -MinMQ 30 -MaxVQ 100 -MinDepthFilter 500 -MinVQ 0 -VQFilter 20 -gVCF True -ReportNoCalls True -CallMNVs True -MaxMNVLength 3 -MaxGapBetweenMNV 1 -Collapse True -ReportRcCounts True -ReportTsCounts True -ThreadByChr True -usestitchedxd true -outputsbfiles false" \
    --psaraPath "/opt/illumina/components/psara/Psara.dll" \
    --additionalPsaraArgs "-roi "$extra_roi_bed_file" -inclusionmodel expand" 

# re-run small variant filter
$SING python "/opt/illumina/step_modules/SmallVariantFilter/small_variant_filter.py" \
    --stepName "SmallVariantFilter" \
    --outputFolder  "/opt/illumina/analysis-folder/Logs_Intermediates/SmallVariantFilter" \
    --dsdmFiles "/opt/illumina/analysis-folder/Logs_Intermediates/VariantCaller/dsdm.json" \
    --inputFolder "/opt/illumina/analysis-folder/Logs_Intermediates/VariantCaller" \
    --inputFileSuffix ".filtered.genome.vcf" \
    --outputFileSuffix  "_SmallVariants.genome.vcf" \
    --dotnetPath "dotnet"\
    --pepePath  "/opt/illumina/components/Pepe/Pepe.dll" \
    --additionalPepeArgs "-genomeDirectory /opt/illumina/resources/genomes/hg19_hardPAR -blacklistBed /opt/illumina/resources/pepe2/pepe_blacklist.bed -baselineFile /opt/illumina/resources/pepe2/pepe_baseline.txt -mode likelihoodratio -depthThreshold 100 -baselineMethod binomial -priorFiles /opt/illumina/resources/pepe2/cosmic_pepe2_prior.vcf -likelihoodRatioNonPriorThreshold 60 -qscoreNonPriorThreshold 60 -qscorePriorThreshold 20 -likelihoodRatioPriorThreshold 20 -backgroundNoiseFrequencyThreshold 0.05 -LoD 0.05  -outputErrorRateTables true"

# re-run phased variants
$SING dotnet /opt/illumina/step_modules/PhasedVariants/PhasedVariantsStepModule.dll \
    --stepName PhasedVariants \
    --componentPath /opt/illumina/components/scylla/Scylla.dll \
    --psaraPath /opt/illumina/components/psara/Psara.dll \
    --runtimePath dotnet \
    --inputBamFolder /mnt/app_result/StitchedRealigned \
    --inputBamFileSuffix .bam \
    --inputVcfFolder /opt/illumina/analysis-folder/Logs_Intermediates/SmallVariantFilter \
    --inputVcfFileSuffix _SmallVariants.genome.vcf \
    --dsdmFiles /opt/illumina/analysis-folder/Logs_Intermediates/SmallVariantFilter/dsdm.json \
    --outputFolder /opt/illumina/analysis-folder/Logs_Intermediates/PhasedVariants \
    --outputFileSuffix .Complex.vcf \
    --additionalComponentArgs '-g /opt/illumina/resources/genomes/hg19_hardPAR -PassingVariantsOnly false --minvq 0 --vqfilter 0 --vffilter 0 --gqfilter 0 --mindpfilter 100' \
    --additionalPsaraArgs '-inclusionmodel=start -roi /opt/illumina/resources/scylla/scylla.interval_list'

# re-run variant matching
$SING dotnet /opt/illumina/step_modules/VariantMatching/VariantMatchingStepModule.dll \
    --stepName VariantMatching \
    --componentPath /opt/illumina/components/yente/Yente.dll \
    --runtimePath dotnet \
    --inputComplexVariantFolder /opt/illumina/analysis-folder/Logs_Intermediates/PhasedVariants \
    --inputComplexVariantFileSuffix .Complex.vcf \
    --inputComplexVariantWhitelistFile /opt/illumina/resources/yente/targeted_complex_variants_list.tsv \
    --inputSmallVariantFolder /opt/illumina/analysis-folder/Logs_Intermediates/SmallVariantFilter \
    --inputSmallVariantFileSuffix _SmallVariants.genome.vcf \
    --dsdmFiles /opt/illumina/analysis-folder/Logs_Intermediates/PhasedVariants/dsdm.json \
    --outputFolder /opt/illumina/analysis-folder/Logs_Intermediates/VariantMatching \
    --outputFileSuffix _MergedSmallVariants.genome.vcf

# re-run contamination
$SING dotnet /opt/illumina/step_modules/Contamination/ContaminationStepModule.dll \
    --stepName "Contamination" \
    --outputFolder "/opt/illumina/analysis-folder/Logs_Intermediates/Contamination" \
    --spoilerPath "/opt/illumina/components/spoiler/Spoiler.dll" \
    --outputFileSuffix ".contamination.json" \
    --inputFileSuffix "_SmallVariants.genome.vcf" \
    --errorFileSuffix "_SmallVariants.genome.vcf.errorRates" \
    --runtimePath "dotnet" \
    --inputFolder "/opt/illumina/analysis-folder/Logs_Intermediates/SmallVariantFilter" \
    --dsdmFiles "/opt/illumina/analysis-folder/Logs_Intermediates/VariantMatching/dsdm.json" \
    --additionalSpoilerParameters "--vaf 0.25 --bedFile /opt/illumina/resources/contamination/Contamination.bed" \

# re-run annotation
$SING python /opt/illumina/step_modules/Annotation/annotation.py \
        --stepName Annotation \
        --outputFolder /opt/illumina/analysis-folder/Logs_Intermediates/Annotation \
        --nirvanaPath /opt/illumina/components/nirvana_3.2.3/Nirvana.dll \
        --additionalNirvanaArgs --disable-recomposition \
        --nirvanaCache /opt/illumina/resources/nirvana_3.2/Cache/26/GRCh37/Both \
        --nirvanaSupplement /opt/illumina/resources/nirvana_3.2/SupplementaryDatabase/3.2.0/GRCh37 \
        --nirvanaRef /opt/illumina/resources/nirvana_3.2/References/6/Homo_sapiens.GRCh37.Nirvana.dat \
        --outputFileSuffix _SmallVariants_Annotated.json.gz \
        --dotNetCorePath dotnet \
        --inputFileSuffix _SmallVariants.genome.vcf \
        --inputFolder /opt/illumina/analysis-folder/Logs_Intermediates/SmallVariantFilter \
        --dsdmFiles /opt/illumina/analysis-folder/Logs_Intermediates/Contamination/dsdm.json

# re-run merged annotation
$SING python /opt/illumina/step_modules/Annotation/annotation.py \
        --stepName MergedAnnotation \
        --outputFolder /opt/illumina/analysis-folder/Logs_Intermediates/MergedAnnotation \
        --nirvanaPath /opt/illumina/components/nirvana_3.2.3/Nirvana.dll \
        --additionalNirvanaArgs --disable-recomposition \
        --nirvanaCache /opt/illumina/resources/nirvana_3.2/Cache/26/GRCh37/Both \
        --nirvanaSupplement /opt/illumina/resources/nirvana_3.2/SupplementaryDatabase/3.2.0/GRCh37 \
        --nirvanaRef /opt/illumina/resources/nirvana_3.2/References/6/Homo_sapiens.GRCh37.Nirvana.dat \
        --outputFileSuffix _MergedVariants_Annotated.json.gz \
        --dotNetCorePath dotnet \
        --inputFileSuffix _MergedSmallVariants.genome.vcf \
        --inputFolder /opt/illumina/analysis-folder/Logs_Intermediates/VariantMatching \
        --dsdmFiles /opt/illumina/analysis-folder/Logs_Intermediates/Annotation/dsdm.json


##############################################################################################
#  Filter JSON output from nirvana into combined variant output file
##############################################################################################

# activate conda env
set +u
conda activate TSO500_post_processing
set -u

# unzip variant JSON
cp "$output_folder"/Logs_Intermediates/MergedAnnotation/"$sample_id"/"$sample_id"_MergedVariants_Annotated.json.gz "$output_folder"
gunzip "$output_folder"/"$sample_id"_MergedVariants_Annotated.json.gz

# format extra variant calls with custom Python script
python "$pipeline_dir"/filter_extra_padding_variants.py "$output_folder"/"$sample_id"_MergedVariants_Annotated.json > "$output_folder"/"$sample_id"_extra_db.tsv

# concatenate extra calls to end of combined output
cat "$SLURM_SUBMIT_DIR"/analysis/"$sample_id"/Results/"$sample_id"/"$sample_id"_CombinedVariantOutput.tsv > "$SLURM_SUBMIT_DIR"/analysis/"$sample_id"/Results/"$sample_id"/"$sample_id"_CombinedVariantOutput_padding.tsv
cat "$output_folder"/"$sample_id"_extra_db.tsv >> "$SLURM_SUBMIT_DIR"/analysis/"$sample_id"/Results/"$sample_id"/"$sample_id"_CombinedVariantOutput_padding.tsv

# remove dsdm files to prevent Gathered results failing downstream - it uses a wildcard and they cause it to match too many folders
rm "$output_folder"/Logs_Intermediates/*/dsdm.json

# deactivate env
set +u
conda deactivate
set -u
