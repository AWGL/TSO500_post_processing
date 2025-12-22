process FILTER_SAMPLE_SHEET {
    cpus 1
    memory 512.MB

    container "132205776083.dkr.ecr.eu-west-2.amazonaws.com/ghcr/awgl/tso500_post_processing:dev"

    input:
    path sample_sheet

    output:
    path "SampleSheet_updated.csv"
    path "samples_correct_order_*_RNA.csv", emit: rna_sample_list
    path "worksheets_rna.txt", emit: rna_worksheets

    script:
    """
    # remove header from samplesheet
    sed -n -e '/Sample_ID,Sample_Name/,\$p' SampleSheet.csv >> SampleSheet_updated.csv

    # make a list of samples and get correct order of samples for each worksheet
    filter_sample_list.py

    # Run dos2unix here
    dos2unix SampleSheet_updated.csv
    """
}