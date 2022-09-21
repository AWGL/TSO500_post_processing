#!/bin/bash

#SBATCH --output=TSO500_DNA_nextflow-%j-%N.out
#SBATCH --error=TSO500_DNA_nextflow-%j-%N.err
#SBATCH --partition=high

# Description: #Kick off nextflow script to run TSO500 DNA analysis
# Use: TSO500_DNA_nextflow.sh <PATH TO RAW READS> <PATH AND NAME OF SAMPLES ORDER DNA FILE> <SEQUENCING RUN ID>

#####################################################################
# Set up
#####################################################################

FASTQ_PATH=$1
SAMPLES_ORDER=$2
SEQID=$3

#####################################################################
# Run Command
#####################################################################

#set +u
#conda activate ctDNA
#set -u

nextflow -C /data/diagnostics/pipelines/ctDNA/config/ctDNA/ctDNA-development/ctDNA.config run /data/diagnostics/pipelines/ctDNA/ctDNA-development/ctDNA.nf \
    --fastqs ${FASTQ_PATH}/\*/\*\{R1.fastq.gz,R2.fastq.gz\} \
    --dna_list ${SAMPLES_ORDER} \
    --publish_dir results \
    --sequencing_run ${SEQID} \
    -with-dag ${SEQID}.png \
    -with-report ${SEQID}.html \
    -work-dir work &> pipeline.log

rm -r work
