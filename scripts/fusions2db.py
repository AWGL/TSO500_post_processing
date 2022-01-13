"""
This script loads in the sample combined variant output, the NTC combined variant output, and the sample all fusions file.
Taking the above files, the script will parse these files in pandas and creates a table that compares the sample to the NTC file
and returns True if the fusion is also found in the NTC.

To run this script you need 4 inputs: './fusions2db.py $tsv_file $ntc_file $all_fusions $output_folder'

"""

import sys
import csv
import pandas as pd
import numpy as np
import argparse


##############################################################################################
#  Define functions
##############################################################################################

def load_report(fusions_file):
    """
    Loads in the csv or tsv file, and extracts the information needed into a dictionary
    Each section starting with [ is a new level in the dictionary

    """
    # open combined results file and save as list
    with open(fusions_file) as f:
        reader = csv.reader(f, delimiter='\t')
        fusion_report_list = list(reader)

    # define empty dict to store output
    fusion_results_dict = {}

    # split combined file into dict 
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

    # return results dict
    return fusion_results_dict


def remove_unwanted_line(fusion_dict):
    """
    Removing unwanted line that contains header info

    """
    # remove header from fusions section
    if fusion_dict['Fusions'][0][0].startswith('Gene'):
        fusion_dict['Fusions'].pop(0)

    # remove header from splice variants
    if fusion_dict['Splice Variants'][0][0].startswith('Gene'):
        fusion_dict['Splice Variants'].pop(0)


def format_fusions(fusion_results_dict, fusion_results_dict_ntc):
    """
    Creating the output table by extracting information from the loaded in files and comparing the NTC to the sample files

    """
    # create empty list to store output
    l = []

    ## format fusions

    # if no fusions, break from loop
    if fusion_results_dict['Fusions'][0][0] == 'NA':
        pass

    # otherwise loop through fusions
    else: 
        for f in fusion_results_dict['Fusions']:

            # check if fusion is in NTC
            in_ntc = False
            for j in fusion_results_dict_ntc['Fusions']:
                if j == f:
                    in_ntc = True

            # add to list
            line = [f[0], '', f[4], f[5], f[3], f[1], f[2], 'Fusion', in_ntc]
            l.append(line)

    ## format splice variants

    # if no splice variants, break from loop
    if fusion_results_dict['Splice Variants'][0][0] == 'NA':
        pass

    # otherwise loop through splice variants
    else:
        for f in fusion_results_dict['Splice Variants']:

            # check if splice variant is in NTC
            in_ntc = False
            for j in fusion_results_dict_ntc['Splice Variants']:
                if j == f:
                    in_ntc = True

            # add to list
            line = [f[0], f[1], f[5], '', f[4], f[2], f[3], 'Splice', in_ntc]
            l.append(line)

    # define headers and save as dataframe
    headers = ['fusion', 'exons', 'reference_reads_1', 'reference_reads_2', 'fusion_supporting_reads', 'left_breakpoint', 'right_breakpoint', 'type', 'in_ntc']
    df = pd.DataFrame.from_records(l, columns=headers)

    return df


def add_extra_columns(df, df_extra):
    """
    Adding extra columns to the dataframe by extracting info from the all fusions file and then merging these two tables together 

    """
    # cut down columns in extra metrics dataframe
    df_extra_cut = df_extra[['Gene A Breakpoint', 'Gene B Breakpoint', 'Alt Pair','Alt Pair Dedup', 'Alt Split', 'Alt Split Dedup', 'Caller', 'Score']]

    # merge columns on breakpoints
    df_merged = pd.merge(left=df, right=df_extra_cut, left_on=['left_breakpoint', 'right_breakpoint'], right_on=['Gene A Breakpoint', 'Gene B Breakpoint'], how='left')

    # drop duplicated breakpoint column and re-header dataframe
    df_merged = df_merged.drop(['Gene A Breakpoint', 'Gene B Breakpoint'], axis=1)
    df_merged_final = df_merged.rename(columns={'Alt Pair': 'spanning_reads', 'Alt Pair Dedup': 'spanning_reads_dedup', 'Alt Split': 'split_reads', 'Alt Split Dedup': 'split_reads_dedup', 'Caller': 'fusion_caller', 'Score': 'fusion_score'})

    # change float columns to integers
    cols = ['spanning_reads', 'spanning_reads_dedup', 'split_reads', 'split_reads_dedup']
    df_merged_final[cols] = df_merged_final[cols].fillna(0.0).astype(int)

    return df_merged_final


##############################################################################################
#  Main script
##############################################################################################

if __name__ == '__main__':

    # load in arguments from argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--tsvfile', '-t', help = 'path to tsv combined variants file')
    parser.add_argument('--ntcfile', '-n', help = 'path to NTC tsv combined variants file')
    parser.add_argument('--allfusions', '-f', help = 'path to tsv all fusions file')
    parser.add_argument('--outfile', '-o', help = 'filename/path for output')
    args = parser.parse_args()

    # load in both sample fusions dict and NTC fusions dict
    fusion_results_dict = load_report(args.tsvfile)
    fusion_results_dict_ntc = load_report(args.ntcfile)

    # remove unwanted lines from both sample and NTC dicts
    remove_unwanted_line(fusion_results_dict)
    remove_unwanted_line(fusion_results_dict_ntc)

    # create a dataframe of all sample fusions and whether or not they are also called in the NTC
    df = format_fusions(fusion_results_dict, fusion_results_dict_ntc)

    # load in the sample *AllFusions.csv file, which contains extra metrics, removes header which all starts with #
    df_extra = pd.read_csv(args.allfusions, comment = '#')

    # merge the original dataframe with the extra metrics
    df_merged_final = add_extra_columns(df, df_extra)

    # save the outout
    sample_id = fusion_results_dict['Analysis Details'][0][1]
    with open(f'{args.outfile}/{sample_id}_fusion_check.csv', 'w') as f:
        df_merged_final.to_csv(f, index=False)
