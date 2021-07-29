import json
import os
import decimal
from decimal import Decimal
import unittest

import pandas as pd
import numpy as np

from coverage2json import parse_referral_type_files, parse_NTC_data, parse_sample_data, create_output_dict


"""
test usage:

type 'python -m unittest' from within a directory which contains the test_coverage2json.py script and test_data folder.

"""


# test NTC values come through

# test genescreen is blank if Melanoma is used (no genescreen files in test_data)
	# can do this at df level and dict level

class test_coverage2json(unittest.TestCase):

	def test_parse_referral_type_genescreen(self):
		"""
		test read of referral type files with genescreen present returns true
		"""
		gene_list, genescreen_df, hotspots_df, genescreen_present = parse_referral_type_files('Tumour', 'test_data/hotspot_coverage')

		self.assertEqual(gene_list[0], 'AR')
		self.assertTrue(genescreen_present)


	def test_parse_referral_type_no_genescreen(self):
		"""
		test read of referral type files with no genescreen present returns false
		"""
		gene_list, genescreen_df, hotspots_df, genescreen_present = parse_referral_type_files('Melanoma', 'test_data/hotspot_coverage')

		self.assertEqual(gene_list[0], 'NRAS')
		self.assertFalse(genescreen_present)


	def test_parse_NTC_data(self):
		"""
		test read of NTC coverage data
		"""
		ntc_gene_df, ntc_region_df = parse_NTC_data('test_data/NTC/depth_of_coverage', 'Melanoma')
		self.assertTrue(len(ntc_region_df) > 0)
		self.assertTrue(len(ntc_gene_df) > 0)


	def test_parse_sample_data(self):
		"""
		test read of sample data
		"""
		sampleid, sample_gene_df, sample_region_df, sample_135_gaps_df, sample_270_gaps_df = parse_sample_data('test_data/sample/depth_of_coverage', 'Melanoma')
		self.assertEqual('sample', sampleid)

		self.assertEqual(sample_gene_df['GENE'][1], 'PTEN')
		self.assertEqual(sample_gene_df['AVG_DEPTH'][1], 510)

		self.assertEqual(sample_region_df['AVG_DEPTH'][1], 780)

		self.assertEqual(sample_135_gaps_df['META'][0], 'NRAS')
		self.assertEqual(sample_270_gaps_df['META'][0], 'NRAS')


	def test_create_output_dict_no_genescreen(self):
		"""
		test creation of output dictionary without genescreen has no genescreen regions in dict
		"""
		gene_list, genescreen_df, hotspots_df, genescreen_present = parse_referral_type_files('Melanoma', 'test_data/hotspot_coverage')
		ntc_gene_df, ntc_region_df = parse_NTC_data('test_data/NTC/depth_of_coverage', 'Melanoma')
		sampleid, sample_gene_df, sample_region_df, sample_135_gaps_df, sample_270_gaps_df = parse_sample_data('test_data/sample/depth_of_coverage', 'Melanoma')
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
		if genescreen_present == True:
			genescreen_region_df = pd.merge(genescreen_df, main_region_df, how = 'left', on = ['CHR', 'START', 'END'])
			hotspots_region_df = pd.merge(hotspots_df, main_region_df, how = 'left', on = ['CHR', 'START', 'END'])
		else:
			hotspots_region_df = main_region_df
			## create genescreen_region_df as nothing to allow for exporting and then importing to create dict
			genescreen_region_df = ''

		output_dict = create_output_dict(gene_list, main_gene_df, genescreen_region_df, hotspots_region_df, sample_270_gaps_df, sample_135_gaps_df, genescreen_present)

		for item in output_dict.values():
			with self.subTest(item = item):
				self.assertEqual(item['genescreen_regions'], [])


	def test_create_output_dict_with_genescreen(self):
		"""
		test creation of output dictionary with genescreen genescreen regions in dict
		"""
		gene_list, genescreen_df, hotspots_df, genescreen_present = parse_referral_type_files('Tumour', 'test_data/hotspot_coverage')
		ntc_gene_df, ntc_region_df = parse_NTC_data('test_data/NTC/depth_of_coverage', 'Tumour')
		sampleid, sample_gene_df, sample_region_df, sample_135_gaps_df, sample_270_gaps_df = parse_sample_data('test_data/sample/depth_of_coverage', 'Tumour')
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
		if genescreen_present == True:
			genescreen_region_df = pd.merge(genescreen_df, main_region_df, how = 'left', on = ['CHR', 'START', 'END'])
			hotspots_region_df = pd.merge(hotspots_df, main_region_df, how = 'left', on = ['CHR', 'START', 'END'])
		else:
			hotspots_region_df = main_region_df
			## create genescreen_region_df as nothing to allow for exporting and then importing to create dict
			genescreen_region_df = ''

		output_dict = create_output_dict(gene_list, main_gene_df, genescreen_region_df, hotspots_region_df, sample_270_gaps_df, sample_135_gaps_df, genescreen_present)
		# print(output_dict)
		self.assertEqual(output_dict['AR']['genescreen_regions'][0][0], 'chrX')

if __name__ == '__main__':
	unittest.main()	