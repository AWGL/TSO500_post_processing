#!/usr/bin/env python

#### README
## To run this script you need 4 inputs: './fusions_check_with_ntc.py $tsv_file $ntc_file $all_fusions $output_folder'

import sys
import csv
import pandas as pd
import yaml
import numpy as np


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


def remove_unwanted_line(fusion_dict):
    if fusion_dict['Fusions'][0][0].startswith('Gene'):
        fusion_dict['Fusions'].pop(0)
    else:
        pass
    
    if fusion_dict['Splice Variants'][0][0].startswith('Gene'):
        fusion_dict['Splice Variants'].pop(0)
    else:
        pass


def format_fusions(fusion_results_dict, in_ntc=False):
    l = []
    headers = ['fusion', 'exons', 'reference_reads_1', 'reference_reads_2', 'fusion_supporting_reads', 'left_breakpoint', 'right_breakpoint', 'type', 'in_ntc']
    
    if fusion_results_dict['Fusions'][0][0] == 'NA':
        df = pd.DataFrame(['No fusions called in sample'])
    
    for n, f in enumerate(fusion_results_dict['Fusions']):
        if fusion_results_dict['Fusions'][0][0] == 'NA':
            break
        line = [f[0], '', f[4], f[5], f[3], f[1], f[2], 'Fusion', in_ntc]
        l.append(line)
        
    for n, f in enumerate(fusion_results_dict['Splice Variants']):
        if fusion_results_dict['Splice Variants'][0][0] == 'NA':
            break
        line = [f[0], f[1], f[5], '', f[4], f[2], f[3], 'Splice', in_ntc]
        l.append(line)
    
    df = pd.DataFrame.from_records(l, columns=headers)
    return df


def create_final_df(df, df_ntc):
    df_final = df.append(df_ntc, ignore_index=True)
    return df_final


def add_extra_columns(df_final, df_extra):
    df_extra_cut = df_extra[['Gene A Breakpoint', 'Gene B Breakpoint', 'Alt Pair','Alt Pair Dedup', 'Alt Split', 'Alt Split Dedup', 'Caller', 'Score']]
    df_merged = pd.merge(left=df_final, right=df_extra_cut, left_on=['left_breakpoint', 'right_breakpoint'], right_on=['Gene A Breakpoint', 'Gene B Breakpoint'], how='left')
    df_merged = df_merged.drop(['Gene A Breakpoint', 'Gene B Breakpoint'], axis=1)
    df_merged_final = df_merged.rename(columns={'Alt Pair': 'spanning_reads', 'Alt Pair Dedup': 'spanning_reads_dedup', 'Alt Split': 'split_reads', 'Alt Split Dedup': 'split_reads_dedup', 'Caller': 'fusion_caller', 'Score': 'fusion_score'})

    cols = ['spanning_reads', 'spanning_reads_dedup', 'split_reads', 'split_reads_dedup']
    df_merged_final[cols] = df_merged_final[cols].fillna(0.0).astype(int)

    return df_merged_final


if __name__ == '__main__':

    tsv_file = sys.argv[1]
    ntc_file = sys.argv[2]
    all_fusions = sys.argv[3]
    output_folder = sys.argv[4]

    fusion_results_dict = load_report(tsv_file)
    fusion_results_dict_ntc = load_report(ntc_file)
    
    sample_id = fusion_results_dict['Analysis Details'][0][1]

    remove_unwanted_line(fusion_results_dict)
    remove_unwanted_line(fusion_results_dict_ntc)

    df = format_fusions(fusion_results_dict)
    df_ntc = format_fusions(fusion_results_dict_ntc, in_ntc=True)

    df_final = create_final_df(df, df_ntc)

    df_extra = pd.read_csv(all_fusions, comment = '#')
    df_merged_final = add_extra_columns(df_final, df_extra)
    



    # df_final["spanning_reads"] = 0
    # df_final["spanning_reads_dedup"] = 0
    # df_final["split_reads"] = 0
    # df_final["split_reads_dedup"] = 0
    # df_final["fusion_caller"] = 0
    # df_final["fusion_score"] = 0

    with open(f'{output_folder}/{sample_id}_fusion_check.csv', 'w') as f:
    	df_merged_final.to_csv(f, index=False)


