#!/usr/bin/env python

#### README
## To run this script you need 4 inputs: './fusions_check_with_ntc.py $tsv_file $ntc_file $all_fusions $output_folder'

## Add tests to GitHub repo - May need re-tweaking with changes I make below...

# This script loads in the sample combined variant output, the NTC combined variant output, and the sample all fusions file.
# Taking the above files, the script will parse these files in pandas and creates a table that compares the sample to the NTC file \
# and returns True if the fusion is also found in the NTC.

# Importing modules needed for this script
import sys
import csv
import pandas as pd
import yaml
import numpy as np
import argparse


# Loads in the csv or tsv file, and extracts the information needed
def load_report(fusions_file):
    with open(fusions_file) as f:
        reader = csv.reader(f, delimiter='\t')
        fusion_report_list = list(reader)
    
    fusion_results_dict = {}
    
    for line in fusion_report_list:
        if line[0] == '':
            pass
        elif line[0].startswith('['):
            try:
                fusion_results_dict[category] = l
            except NameError:
                pass
            category = line[0].strip('[]')
            l = []
        else:
            try:
                l.append(line)
            except NameError:
                pass
    return fusion_results_dict


# Removing unwanted line that contains header info
def remove_unwanted_line(fusion_dict):
    if fusion_dict['Fusions'][0][0].startswith('Gene'):
        fusion_dict['Fusions'].pop(0)
    else:
        pass
    
    if fusion_dict['Splice Variants'][0][0].startswith('Gene'):
        fusion_dict['Splice Variants'].pop(0)
    else:
        pass


# Creating the output table by extracting information from the loaded in files and comparing the NTC to the sample files
def format_fusions(fusion_results_dict, fusion_results_dict_ntc, in_ntc=False):
    l = []
    headers = ['fusion', 'exons', 'reference_reads_1', 'reference_reads_2', 'fusion_supporting_reads', 'left_breakpoint', 'right_breakpoint', 'type', 'in_ntc']
    
    if fusion_results_dict['Fusions'][0][0] == 'NA':
        df = pd.DataFrame(['No fusions called in sample'])
    
    for n, f in enumerate(fusion_results_dict['Fusions']):
        if fusion_results_dict['Fusions'][0][0] == 'NA':
            break
        for j in fusion_results_dict_ntc['Fusions']:
            if j == f:
                line = [f[0], '', f[4], f[5], f[3], f[1], f[2], 'Fusion', True]
                l.append(line)
            else:
                line = [f[0], '', f[4], f[5], f[3], f[1], f[2], 'Fusion', in_ntc]
                l.append(line)
        
    for n, f in enumerate(fusion_results_dict['Splice Variants']):
        if fusion_results_dict['Splice Variants'][0][0] == 'NA':
            break
        for j in fusion_results_dict_ntc['Splice Variants']:
            if j == f:
                line = [f[0], f[1], f[5], '', f[4], f[2], f[3], 'Splice', True]
                l.append(line)
            else:
                line = [f[0], f[1], f[5], '', f[4], f[2], f[3], 'Splice', in_ntc]
                l.append(line)
    
    df = pd.DataFrame.from_records(l, columns=headers)
    return df


# Adding extra columns to the dataframe by extracting info from the all fusions file and then merging these two tables together 
def add_extra_columns(df, df_extra):
    df_extra_cut = df_extra[['Gene A Breakpoint', 'Gene B Breakpoint', 'Alt Pair','Alt Pair Dedup', 'Alt Split', 'Alt Split Dedup', 'Caller', 'Score']]
    df_merged = pd.merge(left=df, right=df_extra_cut, left_on=['left_breakpoint', 'right_breakpoint'], right_on=['Gene A Breakpoint', 'Gene B Breakpoint'], how='left')
    df_merged = df_merged.drop(['Gene A Breakpoint', 'Gene B Breakpoint'], axis=1)
    df_merged_final = df_merged.rename(columns={'Alt Pair': 'spanning_reads', 'Alt Pair Dedup': 'spanning_reads_dedup', 'Alt Split': 'split_reads', 'Alt Split Dedup': 'split_reads_dedup', 'Caller': 'fusion_caller', 'Score': 'fusion_score'})

    cols = ['spanning_reads', 'spanning_reads_dedup', 'split_reads', 'split_reads_dedup']
    df_merged_final[cols] = df_merged_final[cols].fillna(0.0).astype(int)

    return df_merged_final


if __name__ == '__main__':

    # tsv_file = sys.argv[1]
    # ntc_file = sys.argv[2]
    # all_fusions = sys.argv[3]
    # output_folder = sys.argv[4]

    parser = argparse.ArgumentParser()
    parser.add_argument('--tsvfile', '-t', help = 'path to tsv combined variants file')
    parser.add_argument('--ntcfile', '-n', help = 'path to NTC tsv combined variants file')
    parser.add_argument('--allfusions', '-f', help = 'path to tsv all fusions file')
    parser.add_argument('--outfile', '-o', help = 'filename/path for output')
    args = parser.parse_args()


    fusion_results_dict = load_report(args.tsvfile)
    fusion_results_dict_ntc = load_report(args.ntcfile)
    
    sample_id = fusion_results_dict['Analysis Details'][0][1]

    remove_unwanted_line(fusion_results_dict)
    remove_unwanted_line(fusion_results_dict_ntc)

    df = format_fusions(fusion_results_dict, fusion_results_dict_ntc)

    df_extra = pd.read_csv(args.allfusions, comment = '#')
    df_merged_final = add_extra_columns(df, df_extra)

    with open(f'{args.outfile}/{sample_id}_fusion_check.csv', 'w') as f:
    	df_merged_final.to_csv(f, index=False)


