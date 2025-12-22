process FASTQC {
    tag "${sample_id}_${lane_id}_${read}"

    container "community.wave.seqera.io/library/fastqc:0.11.9--fd0125189547cb9c"

    publishDir "${params.output_dir}/analysis/${sample_id}/FastQC/"

    input:
    tuple val(sample_id), val(lane_id), val(read), path(fastq)

    output:
    path("**")
    tuple val(sample_id), path("${sample_id}_${lane_id}_${read}_fastqc.txt"), emit: fastqc_summary

    script:
    """
    fastqc --extract ${fastq}
    mv */summary.txt ${sample_id}_${lane_id}_${read}_fastqc.txt
    """
}