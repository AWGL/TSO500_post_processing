#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --output=%j-%N-TSO500_extra_padding.out
#SBATCH --error=%j-%N-TSO500_extra_padding.err
#SBATCH --partition=high
#SBATCH --cpus-per-task=24

# Description: Run parts of illumina app to call variants between +/-2 to +/-5 (outside of usualy TSO500 ROI)
# Author:      AWMGS
# Mode:        BY_SAMPLE
# Use:         sbatch within run directory, pass in sample_id

version=2.2.0

cd "$SLURM_SUBMIT_DIR"


app_logs_intermediates="$SLURM_SUBMIT_DIR"/analysis/"$sample_id"/Logs_Intermediates
output_folder="$SLURM_SUBMIT_DIR"/analysis/"$sample_id"/padding

mkdir -p $output_folder


module purge
module load singularity
. ~/.bashrc
module load anaconda

# catch fails early and terminate
set -euo pipefail

pipeline_dir=/data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-"$version"
resources="$pipeline_dir"/resources/

# this is only visible inside singularity container
extra_roi_bed_file=/opt/illumina/resources/TSO_extra_padding_chr.interval_list


# define singularity exec command
SING="singularity exec --bind "$app_logs_intermediates":/mnt/app_result,"$output_folder":/opt/illumina/analysis-folder,"$resources":/opt/illumina/resources /data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-2.2.0/trusight-oncology-500-ruo.img"


# run each section in singularity
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



$SING python /opt/illumina/step_modules/Annotation/annotation.py \
        --stepName Annotation \
        --outputFolder /opt/illumina/analysis-folder/Logs_Intermediates/Annotation \
        --nirvanaPath /opt/illumina/components/nirvana_3.2.3/Nirvana.dll \
         \
         \
        --additionalNirvanaArgs --disable-recomposition \
        --nirvanaCache /opt/illumina/resources/nirvana_3.2/Cache/26/GRCh37/Both \
        --nirvanaSupplement /opt/illumina/resources/nirvana_3.2/SupplementaryDatabase/3.2.0/GRCh37 \
        --nirvanaRef /opt/illumina/resources/nirvana_3.2/References/6/Homo_sapiens.GRCh37.Nirvana.dat \
        --outputFileSuffix _SmallVariants_Annotated.json.gz \
        --dotNetCorePath dotnet \
        --inputFileSuffix _SmallVariants.genome.vcf \
        --inputFolder /opt/illumina/analysis-folder/Logs_Intermediates/SmallVariantFilter \
        --dsdmFiles /opt/illumina/analysis-folder/Logs_Intermediates/Contamination/dsdm.json


$SING python /opt/illumina/step_modules/Annotation/annotation.py \
        --stepName MergedAnnotation \
        --outputFolder /opt/illumina/analysis-folder/Logs_Intermediates/MergedAnnotation \
        --nirvanaPath /opt/illumina/components/nirvana_3.2.3/Nirvana.dll \
         \
         \
        --additionalNirvanaArgs --disable-recomposition \
        --nirvanaCache /opt/illumina/resources/nirvana_3.2/Cache/26/GRCh37/Both \
        --nirvanaSupplement /opt/illumina/resources/nirvana_3.2/SupplementaryDatabase/3.2.0/GRCh37 \
        --nirvanaRef /opt/illumina/resources/nirvana_3.2/References/6/Homo_sapiens.GRCh37.Nirvana.dat \
        --outputFileSuffix _MergedVariants_Annotated.json.gz \
        --dotNetCorePath dotnet \
        --inputFileSuffix _MergedSmallVariants.genome.vcf \
        --inputFolder /opt/illumina/analysis-folder/Logs_Intermediates/VariantMatching \
        --dsdmFiles /opt/illumina/analysis-folder/Logs_Intermediates/Annotation/dsdm.json



# filter JSON output from nirvana

set +u
conda activate TSO500_post_processing
set -u

cp "$output_folder"/Logs_Intermediates/MergedAnnotation/"$sample_id"/"$sample_id"_MergedVariants_Annotated.json.gz "$output_folder"
gunzip "$output_folder"/"$sample_id"_MergedVariants_Annotated.json.gz

python "$pipeline_dir"/filter_extra_padding_variants.py "$output_folder"/"$sample_id"_MergedVariants_Annotated.json > "$output_folder"/"$sample_id"_extra_db.tsv


set +u
conda deactivate
set -u

