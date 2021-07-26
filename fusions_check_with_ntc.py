#!/usr/bin/env python

import sys
import csv
import pandas as pd
import yaml


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
    headers = ['fusion', 'reference_reads', 'fusion_supporting_reads', 'left_breakpoint', 'right_breakpoint', 'type', 'in_ntc']
    
    if fusion_results_dict['Fusions'][0][0] == 'NA':
        df = pd.DataFrame(['No fusions called in sample'])
    
    for n, f in enumerate(fusion_results_dict['Fusions']):
        if fusion_results_dict['Fusions'][0][0] == 'NA':
            break
        line = [f[0], f[5], f[3], f[1], f[2], 'Fusion', in_ntc]
        l.append(line)
        
    for n, f in enumerate(fusion_results_dict['Splice Variants']):
        if fusion_results_dict['Splice Variants'][0][0] == 'NA':
            break
        line = [f[0], f[5], f[3], f[1], f[2], 'Splice', in_ntc]
        l.append(line)
    
    df = pd.DataFrame.from_records(l, columns=headers)
    return df


def create_final_df(df, df_ntc):
    df_final = df.append(df_ntc, ignore_index=True)
    return df_final


if __name__ == '__main__':

    tsv_file = sys.argv[1]
    ntc_file = sys.argv[2]
    output_folder = sys.argv[3]

    fusion_results_dict = load_report(tsv_file)
    fusion_results_dict_ntc = load_report(ntc_file)
    
    sample_id = fusion_results_dict['Analysis Details'][0][1]

    remove_unwanted_line(fusion_results_dict)
    remove_unwanted_line(fusion_results_dict_ntc)

    df = format_fusions(fusion_results_dict)
    df_ntc = format_fusions(fusion_results_dict_ntc, in_ntc=True)

    df_final = create_final_df(df, df_ntc)

    df_final["spanning_reads"] = 0
    df_final["split_reads"] = 0
    df_final["fusion_caller"] = 0
    df_final["fusion_score"] = 0

    with open(f'{output_folder}/{sample_id}_fusion_check.csv', 'w') as f:
    	df_final.to_csv(f, index=False)


