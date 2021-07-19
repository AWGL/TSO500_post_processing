#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --output=depthofcoverage-%N-%j.output
#SBATCH --error=depthofcoverage-%N-%j.error
#SBATCH --partition=low
#SBATCH --cpus-per-task=1
cd "$SLURM_SUBMIT_DIR"

runid=$1
sampleid=$2

bam_path=/Output/validations/TSO500/"$runid"/"$sampleid"/Logs_Intermediates/StitchedRealigned/"$sampleid"/
vendor_capture_bed=/data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-2.2.0/vendorCaptureBed_100pad_updated.bed
depth_path=


#### reheader the bams to local area


java -jar /Apps/wren/picard/2.21.6/bin/picard.jar AddOrReplaceReadGroups \
I= "$bam_path"/"$sampleid".bam \
O="$bam_path"/"$sampleid"_add_rg.bam \
RGID=4 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=20



#### index new bam
source /home/v.la094519/miniconda3/bin/activate TSO500_post_processing

samtools index "$bam_path"/"$sampleid"_add_rg.bam "$bam_path"/"$sampleid"_add_rg.bam.bai



#### run depth of coverage

gatk DepthOfCoverage \
-I "$bam_path""$sampleid"_add_rg.bam \
-L "$vendor_capture_bed" \
-R /data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-2.2.0/resources/genomes/hg19_hardPAR/genome.fa \
-O "$depth_path"/"$sampleid"_depth_of_coverage


# change to tab delimited and remove colon from column 1
sed 's/:/\t/g' "$depth_path"/"$sampleid"_depth_of_coverage \
| sed 's/,/\t/g' | grep -v 'Locus' \
| sort -k1,1 -k2,2n | bgzip \
> "$depth_path"/"$sampleid"_depth_of_coverage.gz


# tabix file
tabix -b 2 -e 2 -s 1 "$depth_path"/"$sampleid"_depth_of_coverage.gz

conda deactivate
