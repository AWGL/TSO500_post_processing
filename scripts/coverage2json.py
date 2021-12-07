import json
import os
import argparse
import decimal
from decimal import Decimal

import pandas as pd
import numpy as np

"""
- script to format the output of CoverageCalculatorPy (.totalcoverage, .coverage) into a JSON format for importing into the new database
- for use as part of the TSO500 pipeline
- usage: python coverage2json.py -r <referral type> -g <hotspots_coverage folder path> -s <sample Coverage_results folder> -n <ntc Coverage_results folder>
- output: {sampleid}_{referral}_db_coverage.json file
{ 
	"gene1"	:	{
		"average_depth"		:	<integer>,
		"percent_135"		:	<integer>,
		"percent_270"		:	<integer>,
		"average_ntc"		:	<integer>,
		"percent_ntc"		:	<integer>,
		"genescreen_regions"	:	[
			[<chr>, <start>, <end>, "gene1"<gene info>, <avdepth>, <%@270>, <%@135>, <ntc depth>, <%ntc>],
			[<chr>, <start>, <end>, "gene1"<gene info>, <avdepth>, <%@270>, <%@135>, <ntc depth>, <%ntc>]
		],
		"hotspot_regions"	:	[
			[<chr>, <start>, <end>, "gene1"<gene info>, <avdepth>, <%@270>, <%@135>, <ntc depth>, <%ntc>],
			[<chr>, <start>, <end>, "gene1"<gene info>, <avdepth>, <%@270>, <%@135>, <ntc depth>, <%ntc>]
		],
		"gaps_135"	:	[],
		"gaps_270"	:	[
			[<chr>, <start>, <end>, "gene1"<gene info>, <cosmic>],
			[<chr>, <start>, <end>, "gene1"<gene info>, <cosmic>]
		]
	}
	"gene2"	:	{
		"average_depth"		:	<integer>,
		"percent_135"		:	<integer>,
		"percent_270"		:	<integer>,
		"average_ntc"		:	<integer>,
		"percent_ntc"		:	<integer>,
		"genescreen_regions"	:	[
			[<chr>, <start>, <end>, "gene2"<gene info>, <avdepth>, <%@270>, <%@135>, <ntc depth>, <%ntc>],
			[<chr>, <start>, <end>, "gene2"<gene info>, <avdepth>, <%@270>, <%@135>, <ntc depth>, <%ntc>]
		],
		"hotspot_regions"	:	[
			[<chr>, <start>, <end>, "gene2"<gene info>, <avdepth>, <%@270>, <%@135>, <ntc depth>, <%ntc>],
			[<chr>, <start>, <end>, "gene2"<gene info>, <avdepth>, <%@270>, <%@135>, <ntc depth>, <%ntc>]
		],
		"gaps_135"	:	[],
		"gaps_270"	:	[
			[<chr>, <start>, <end>, "gene2"<gene info>, <cosmic>],
			[<chr>, <start>, <end>, "gene2"<gene info>, <cosmic>]
		]
	}

}
"""

#####################################
#####         Functions         #####
#####################################

def np_encoder(object):
	"""
	- function needed to allow json.dump to parse np values correctly
	"""
	if isinstance(object, np.generic):
		return object.item()


def parse_referral_type_files(referral_type, groups_folder):
	"""
	- Function to parse referral type bed and group files
	- input: referral type (eg: 'Thyroid'), pathway to hotspot_coverage folder containing the groups/bed files
	- output:
	1) list of distinct gene names from <referral>_combined.groups file
	2) df of <referral>_genescreen.bed if present, blank if not present for referral type
	3) df of <referral>_hotspots.bed
	4) True or False value for if <referral>_genescreen.bed is present
	"""

	## create zero genescreen counter
	genescreen_count = 0

	dir_list = os.listdir(groups_folder)

	for file in dir_list:

		filepath = os.path.join(groups_folder,file)

		## parse groups file
		if file == f'{referral_type}_combined.groups':
			groups_df = pd.read_csv(filepath, sep = '\t', index_col = False)
			gene_df = groups_df.drop_duplicates()
			gene_list = gene_df['GENE'].values.tolist()

		## parse genescreen.bed file
		elif file == f'{referral_type}_genescreen.bed':
			genescreen_count += 1
			genescreen_df = pd.read_csv(filepath, sep = '\t', names = ['CHR', 'START', 'END', 'INFO'], index_col = False)
			genescreen_df.drop(columns = ['INFO'], inplace = True)

		## parse hotspots.bed file
		elif file == f'{referral_type}_hotspots.bed':
			hotspots_df = pd.read_csv(filepath, sep = '\t', names = ['CHR', 'START', 'END', 'INFO'], index_col = False)
			hotspots_df.drop(columns = ['INFO'], inplace = True)

	## set genescreen variable as not all referral types have genescreen bed/groups
	if genescreen_count > 0:
		genescreen_present = True

	else:
		genescreen_present = False
		## create blank dataframe to return when no genescreen <referral>.bed file present
		genescreen_df = pd.DataFrame()

	return gene_list, genescreen_df, hotspots_df, genescreen_present


def parse_NTC_data(NTC_coverage_folder, referral_type):
	"""
	- function to parse NTC data from .totalCoverage and .coverage files from CoverageCalculatorPy
	- input: path to NTC Coverage_results folder and referral type
	- output: 
	1) df of .totalcoverage gene level (gene id, ntc avg depth)
	2) df of .coverage region level (chr, start, end, meta, ntc avg depth)
	"""

	## parse .totalCoverage and .coverage for NTC_<referral>_combined
	cov_filepath = os.path.join(NTC_coverage_folder,'hotspot_coverage_270x')
	dir_list = os.listdir(cov_filepath)

	## rip sampleid from filename
	sampleid = dir_list[0].split('_')[0]

	for file in dir_list:
		filepath = os.path.join(cov_filepath, file)

		## parse .totalCoverage
		if file == f'{sampleid}_{referral_type}_combined.totalCoverage':
			ntc_gene_df = pd.read_csv(filepath, sep = '\t', index_col = False)

			## drop final row in df as this is the total average across total panel
			ntc_gene_df.drop(ntc_gene_df.tail(1).index, inplace = True)


			## get gene name positioning in 'feature' column (can change dependent on referral type)
			feature_line = ntc_gene_df['FEATURE'][0].split('_')
			gene_pos = len(feature_line) - 2

			## create 'gene' list and add to DF as a column
			gene_list = []
			for line in ntc_gene_df['FEATURE']:
				splitline = line.split('_')
				## if-loop to get around final row in the DF not containing gene name
				if len(splitline) > gene_pos:
					gene_list.append(splitline[gene_pos])
				else:
					gene_list.append()

			ntc_gene_df.insert(0,'GENE', gene_list)

			## drop irrelevant perc_coverage and feature column
			ntc_gene_df.drop(columns = ['PERC_COVERAGE@270', 'FEATURE'], inplace = True)

			ntc_gene_df.rename(columns={'AVG_DEPTH':'NTC_AVG_DEPTH'}, inplace = True)

		## parse .coverage
		elif file == f'{sampleid}_{referral_type}_combined.coverage':
			ntc_region_df = pd.read_csv(filepath, sep = '\t', index_col = False)
			ntc_region_df.drop(columns = ['PERC_COVERAGE@270'], inplace = True)
			ntc_region_df.rename(columns={'AVG_DEPTH':'NTC_AVG_DEPTH'}, inplace = True)

	ntc_gene_df['NTC_AVG_DEPTH'] = ntc_gene_df['NTC_AVG_DEPTH'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))
	ntc_region_df['NTC_AVG_DEPTH'] = ntc_region_df['NTC_AVG_DEPTH'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))

	return ntc_gene_df, ntc_region_df


def parse_sample_data(sample_coverage_folder, referral_type):
	"""
	- function to parse sample data from .totalCoverage and .coverage files from CoverageCalculatorPy
	- input: path to sample Coverage_results folder and referral type
	- output: 
	1) sampleid
	2) df of .totalcoverage gene level (gene id, avg depth, % cov 270, % cov 135)
	3) df of .coverage region level (chr, start, end, meta, avg depth, % cov 270, % cov 135)
	4) df of 135 .gaps (chr, start, end, meta, cosmic)
	5) df of 270 .gaps (chr, start, end, meta, cosmic)
	"""

	## parse 270x
	cov_270_filepath = os.path.join(sample_coverage_folder,'hotspot_coverage_270x')
	dir_list = os.listdir(cov_270_filepath)

	## rip sampleid from filename
	sampleid = dir_list[0].split('_')[0]

	for file in dir_list:
		filepath = os.path.join(cov_270_filepath, file)

		## parse .totalCoverage
		if file == f'{sampleid}_{referral_type}_combined.totalCoverage':
			sample_270_gene_df = pd.read_csv(filepath, sep = '\t', index_col = False)

			## drop final row in df as this is the total average across total panel
			sample_270_gene_df.drop(sample_270_gene_df.tail(1).index, inplace = True)

			## get gene name positioning in 'feature' column (can change dependent on referral type)
			feature_line = sample_270_gene_df['FEATURE'][0].split('_')
			gene_pos = len(feature_line) - 2

			## create 'gene' list and add to DF as a column
			gene_list = []
			for line in sample_270_gene_df['FEATURE']:
				splitline = line.split('_')
				## if loop to get around final row in the DF not containing gene name
				if len(splitline) > gene_pos:
					gene_list.append(splitline[gene_pos])
				else:
					gene_list.append()

			sample_270_gene_df.insert(0,'GENE', gene_list)

			## drop irrelevant feature column
			sample_270_gene_df.drop(columns = ['FEATURE'], inplace = True)

		## parse .coverage
		elif file == f'{sampleid}_{referral_type}_combined.coverage':
			sample_270_region_df = pd.read_csv(filepath, sep = '\t', index_col = False)

		## parse 270 gaps and add cosmic column
		elif file == f'{sampleid}_{referral_type}_hotspots.gaps':
			sample_270_gaps_df = pd.read_csv(filepath, sep = '\t', names = ['CHR', 'START', 'END', 'META'], header = None, index_col = False)
			sample_270_gaps_df['COSMIC'] = 'N/A'

	## parse 135x
	cov_135_filepath = os.path.join(sample_coverage_folder,'hotspot_coverage_135x')
	dir_list = os.listdir(cov_135_filepath)

	## rip sampleid from filename
	sampleid = dir_list[0].split('_')[0]

	for file in dir_list:
		filepath = os.path.join(cov_135_filepath, file)

		## parse .totalCoverage
		if file == f'{sampleid}_{referral_type}_combined.totalCoverage':
			sample_135_gene_df = pd.read_csv(filepath, sep = '\t', index_col = False)

			## drop final row in df as this is the total average across total panel
			sample_135_gene_df.drop(sample_135_gene_df.tail(1).index, inplace = True)

			## get gene name positioning in 'feature' column (can change dependent on referral type)
			feature_line = sample_135_gene_df['FEATURE'][0].split('_')
			gene_pos = len(feature_line) - 2

			## create 'gene' list and add to DF as a column
			gene_list = []
			for line in sample_135_gene_df['FEATURE']:
				splitline = line.split('_')
				gene_list.append(splitline[gene_pos])


			sample_135_gene_df.insert(0,'GENE', gene_list)

			## drop irrelevant feature column
			sample_135_gene_df.drop(columns = ['FEATURE'], inplace = True)

		## parse .coverage
		elif file == f'{sampleid}_{referral_type}_combined.coverage':
			sample_135_region_df = pd.read_csv(filepath, sep = '\t', index_col = False)

		## parse 135 gaps
		elif file == f'{sampleid}_{referral_type}_hotspots.gaps':
			sample_135_gaps_df = pd.read_csv(filepath, sep = '\t', names = ['CHR', 'START', 'END', 'META'], header = None, index_col = False)
			sample_135_gaps_df['COSMIC'] = 'N/A'

	## join 270 and 135 tables
	sample_region_df = pd.merge(sample_270_region_df, sample_135_region_df, how = 'outer', on = ['CHR', 'START', 'END', 'META', 'AVG_DEPTH'])
	sample_gene_df = pd.merge(sample_270_gene_df, sample_135_gene_df, how = 'outer', on = ['GENE', 'AVG_DEPTH'])

	sample_gene_df['AVG_DEPTH'] = sample_gene_df['AVG_DEPTH'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))
	sample_gene_df['PERC_COVERAGE@270'] = sample_gene_df['PERC_COVERAGE@270'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))
	sample_gene_df['PERC_COVERAGE@135'] = sample_gene_df['PERC_COVERAGE@135'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))
	sample_region_df['AVG_DEPTH'] = sample_region_df['AVG_DEPTH'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))
	sample_region_df['PERC_COVERAGE@270'] = sample_region_df['PERC_COVERAGE@270'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))
	sample_region_df['PERC_COVERAGE@135'] = sample_region_df['PERC_COVERAGE@135'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))

	return sampleid, sample_gene_df, sample_region_df, sample_135_gaps_df, sample_270_gaps_df


def create_output_dict(gene_list, main_gene_df, genescreen_region_df, hotspots_region_df, sample_270_gaps_df, sample_135_gaps_df, genescreen_present):
	'''
	- function to collate all dataframes into one nested dictionary to export to JSON format
	- input: gene list, gene level df, genescreen region df, hotspot region df, 270 gaps df, 135 gaps df
	- output: dictionary formatted as JSON specified at top of script
	'''

	output_dict = {}
	for gene in gene_list:
		output_dict[gene] = {}

	for key in output_dict.keys():
		output_dict[key]['average_depth'] = main_gene_df.at[key, 'AVG_DEPTH']
		output_dict[key]['percent_135'] = main_gene_df.at[key, 'PERC_COVERAGE@135']
		output_dict[key]['percent_270'] = main_gene_df.at[key, 'PERC_COVERAGE@270']
		output_dict[key]['average_ntc'] = main_gene_df.at[key, 'NTC_AVG_DEPTH']
		output_dict[key]['percent_ntc'] = main_gene_df.at[key, 'PERC_NTC_DEPTH']

		## genescreen
		if genescreen_present:
			filtered_genescreen_region_list =[]
			genescreen_region_list = genescreen_region_df.values.tolist()
			for item in genescreen_region_list:
				# check if gene is the gene in the 4th column then add
				if key == item[3].split('(')[0]:
					filtered_genescreen_region_list.append(item)
			output_dict[key]['genescreen_regions'] = filtered_genescreen_region_list
		else:
			output_dict[key]['genescreen_regions'] = []

		## hotspots
		filtered_hotspot_region_list = []
		hotspot_region_list = hotspots_region_df.values.tolist()
		for item in hotspot_region_list:
			# check if gene is the gene in the 4th column then add
			if key == item[3].split('(')[0]:
				filtered_hotspot_region_list.append(item)
		output_dict[key]['hotspot_regions'] = filtered_hotspot_region_list

		## check if any gaps and add if there are
		if len(sample_135_gaps_df) > 0:
			filtered_135_gaps_list = []
			sample_135_gaps_list = sample_135_gaps_df.values.tolist()
			for item in sample_135_gaps_list:
				# check if gene is the gene in the 4th column then add
				if key == item[3].split('(')[0]:
					filtered_135_gaps_list.append(item)
			output_dict[key]['gaps_135'] = filtered_135_gaps_list
		else:
			output_dict[key]['gaps_135'] = []
		if len(sample_270_gaps_df) > 0:
			filtered_270_gaps_list = []
			sample_270_gaps_list = sample_270_gaps_df.values.tolist()
			for item in sample_270_gaps_list:
				# check if gene is the gene in the 4th column then add
				if key == item[3].split('(')[0]:
					filtered_270_gaps_list.append(item)
			output_dict[key]['gaps_270'] = filtered_270_gaps_list
		else:
			output_dict[key]['gaps_270'] = []
	
	return output_dict



#########################################
#####           Programme           #####
#########################################


if __name__ == '__main__':


	## args
	parser = argparse.ArgumentParser()
	parser.add_argument('--referral','-r', help = 'Referral type. eg: Thyroid')
	parser.add_argument('--groups_folder','-g', help = 'pathway to hotspots_coverage folder')
	parser.add_argument('--ntc_coverage','-n', help = 'pathway to NTC Coverage_results folder')
	parser.add_argument('--sample_coverage','-s', help = 'pathway to sample Coverage_results folder')
	parser.add_argument('--outfile', '-o', help = 'pathway/file of output')
	args = parser.parse_args()


	## set decimal settings
	decimal.getcontext().rounding = decimal.ROUND_HALF_UP


	### parse referral_type group/bed files
	gene_list, genescreen_df, hotspots_df, genescreen_present = parse_referral_type_files(args.referral, args.groups_folder)

	### parse NTC sample 270x for average depth per gene (.totalCoverage) and per region (.coverage)
	ntc_gene_df, ntc_region_df = parse_NTC_data(args.ntc_coverage, args.referral)


	### parse sample 135x and 270x files
	sampleid, sample_gene_df, sample_region_df, sample_135_gaps_df, sample_270_gaps_df = parse_sample_data(args.sample_coverage, args.referral)


	### create json pieces
	## join NTC data to main df
	main_region_df = pd.merge(sample_region_df, ntc_region_df, how = 'outer', on = ['CHR', 'START', 'END', 'META'])
	main_gene_df = pd.merge(sample_gene_df, ntc_gene_df, how = 'outer', on = ['GENE'])

	## change index of gene df to be the gene name for iloc later
	main_gene_df.set_index('GENE', inplace = True)

	## create and format percent NTC column
	main_region_df['PERC_NTC_DEPTH'] = (main_region_df['NTC_AVG_DEPTH'] / main_region_df['AVG_DEPTH']) * 100
	main_region_df['PERC_NTC_DEPTH'] = main_region_df['PERC_NTC_DEPTH'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))
	main_gene_df['PERC_NTC_DEPTH'] = (main_gene_df['NTC_AVG_DEPTH'] / main_gene_df['AVG_DEPTH']) * 100
	main_gene_df['PERC_NTC_DEPTH'] = main_gene_df['PERC_NTC_DEPTH'].apply(lambda x: None if np.isnan(x) else int((Decimal(str(x)).quantize(Decimal('1')))))

	## seperate main df into hotspot and genescreen if genescreen is present
	if genescreen_present:
		genescreen_region_df = pd.merge(genescreen_df, main_region_df, how = 'left', on = ['CHR', 'START', 'END'])
		hotspots_region_df = pd.merge(hotspots_df, main_region_df, how = 'left', on = ['CHR', 'START', 'END'])
	else:
		hotspots_region_df = main_region_df
		
		## create genescreen_region_df as nothing to allow for exporting and then importing to create dict
		genescreen_region_df = ''


	### create output dict
	output_dict = create_output_dict(gene_list, main_gene_df, genescreen_region_df, hotspots_region_df, sample_270_gaps_df, sample_135_gaps_df, genescreen_present)


	## export dict to JSON
	with open(args.outfile,'w') as f:
		json.dump(output_dict, f, indent = 4, default = np_encoder)
