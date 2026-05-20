process FILTER_SAMPLE_SHEET {
    cpus 1
    memory 512.MB

    input:
    path sample_sheet

    output:
    path "SampleSheet_updated.csv"
    path "samples_correct_order_*_RNA.csv", emit: rna_sample_list
    path "worksheets_rna.txt", emit: rna_worksheets

    script:
    """
    # make a list of samples and get correct order of samples for each worksheet
    filter_sample_list.py

    """
}