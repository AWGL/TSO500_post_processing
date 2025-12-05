#!/usr/bin/env python

### A script to visualise the read counts for runs ###

import pandas as pd
import seaborn as sns
import glob
import sys

def read_count_visualisation(sample_sheet, list_of_demultiplex_stats):
	dfs = {}
	for demultiplex_stat in list_of_demultiplex_stats:
		dfs[demultiplex_stat] = pd.read_csv(demultiplex_stat)
	main_df = dfs[list(dfs.keys())[0]][['SampleID', 'Index', '# Reads']]
	main_df = main_df.rename(columns={'SampleID': 'sample_id', 'Index': 'index', '# Reads': 'reads_lane_{lane}'.format(lane=dfs[list(dfs.keys())[0]]['Lane'][0])})
	main_df['total'] = main_df['reads_lane_{lane}'.format(lane=dfs[list(dfs.keys())[0]]['Lane'][0])]
	for i in list(dfs.keys())[1:]:
		main_df['reads_lane_{lane}'.format(lane=dfs[i]['Lane'][0])] = pd.Series(dfs[i]['# Reads'])
		main_df['total'] += dfs[i]['# Reads']
	 
	main_df = main_df[main_df.sample_id != 'Undetermined']
		
	with open(sample_sheet,'r') as f:
		for num, line in enumerate(f):
		 # check if the current line
		 # starts with "[Data]"
			if line.startswith("[Data]"):
				df_without_header = pd.read_csv(sample_sheet, skiprows=num+1)
	df_without_header = df_without_header[['Sample_ID', 'Sample_Plate']]
	df_without_header = df_without_header.rename(columns={'Sample_ID': 'sample_id', 'Sample_Plate': 'worksheet_id'})
	merged_df = pd.merge(left=main_df, right=df_without_header, left_on=['sample_id'], right_on=['sample_id'], how='left')
	worksheet_ids = merged_df.worksheet_id.unique()
	group_by_worksheet_df = merged_df.groupby(merged_df.worksheet_id)
	for worksheet in worksheet_ids:
		sub_df = group_by_worksheet_df.get_group(worksheet).sort_values(by=['sample_id'])
		g = sns.FacetGrid(sub_df, col="worksheet_id", col_wrap=1, height=10, aspect=1.5)
		g.map_dataframe(sns.barplot, x="sample_id", y="total")  
		g.set(ylim=(0, 200000000))
		g.set_axis_labels("Sample ID", "Total Reads")
		g.set_xticklabels(rotation=45) 
		g.savefig('read_count_{worksheet}.png'.format(worksheet=worksheet))

if __name__ == '__main__':

	demultiplex_stats = glob.glob('Demultiplex_Output/Logs_Intermediates/FastqGeneration/Reports/Lane*/Demultiplex_Stats.csv')
	sample_sheet = 'Demultiplex_Output/Logs_Intermediates/FastqGeneration/SampleSheet_combined.csv'
	try:
		read_count_visualisation(sample_sheet, demultiplex_stats)
	except IndexError:
		print("ERROR! The Illumina app has not demultiplexed the data. You will need to demultiplex with alternate software")
		sys.exit(1)