#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --output=depthofcoverage-%N-%j.output
#SBATCH --error=depthofcoverage-%N-%j.error
#SBATCH --partition=low
#SBATCH --cpus-per-task=1
cd "$SLURM_SUBMIT_DIR"

. ~/.bashrc
module load anaconda


sampleid=$1
version=$2


minimum_coverage=270,135

bam_path=./"$sampleid"/Logs_Intermediates/StitchedRealigned/"$sampleid"/
vendor_capture_bed=/data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-"$version"/vendorCaptureBed_100pad_updated.bed

mkdir ./"$sampleid"/Results/"$sampleid"/depth_results
depth_path=./"$sampleid"/Results/"$sampleid"/depth_results
coverage_bed_files_path=/data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-"$version"/hotspot_coverage



conda activate TSO500_post_processing


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

samtools index "$bam_path"/"$sampleid"_add_rg.bam "$bam_path"/"$sampleid"_add_rg.bam.bai


#### run depth of coverage

gatk DepthOfCoverage \
-I "$bam_path""$sampleid"_add_rg.bam \
-L "$vendor_capture_bed" \
-R /data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-"$version"/resources/genomes/hg19_hardPAR/genome.fa \
-O "$depth_path"/"$sampleid"_depth_of_coverage


# change to tab delimited and remove colon from column 1
sed 's/:/\t/g' "$depth_path"/"$sampleid"_depth_of_coverage \
| sed 's/,/\t/g' | grep -v 'Locus' \
| sort -k1,1 -k2,2n | bgzip \
> "$depth_path"/"$sampleid"_depth_of_coverage.gz


# tabix file
tabix -b 2 -e 2 -s 1 "$depth_path"/"$sampleid"_depth_of_coverage.gz

conda deactivate


#### run coverageCalculator


conda deactivate
conda activate CoverageCalculatorPy



IFS=',' read -r -a COV <<< "${minimum_coverage}"


for min_coverage in "${COV[@]}"; do
    hscov_outdir=hotspot_coverage_"$min_coverage"x

    for bedFile in "$coverage_bed_files_path"/*.bed; do

        name=$(echo $(basename $bedFile) | cut -d"." -f1)
        echo $name


        python /data/diagnostics/apps/CoverageCalculatorPy/CoverageCalculatorPy-v1.1.0/CoverageCalculatorPy.py \
            -B "$coverage_bed_files_path"/"$name".bed \
            -D "$depth_path"/"$sampleid"_depth_of_coverage.gz \
            --depth "$min_coverage" \
            --padding 0 \
            --groupfile "$coverage_bed_files_path"/"$name".groups \
            --outname "$sampleid"_"$name"_output \
            --outdir  "$depth_path"/Coverage_results/"$hscov_outdir"/

    done
done

