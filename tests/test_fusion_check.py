#!/usr/bin/env python

# from fusions_check_with_ntc_extra_columns import *
from fusions_check_with_ntc import *
#load_report, remove_unwanted_line, format_fusions, create_final_df, add_extra_columns

import unittest
import pathlib as pl

class TestFusions(unittest.TestCase):

	def test_load_report_with_sample(self):
		fusion_results_dict = load_report('testData/21M81042-1_CombinedVariantOutput.tsv')
		self.assertEqual(len(fusion_results_dict.keys()), 7)
		self.assertEqual(len(fusion_results_dict['Fusions'][0]), 6)
		self.assertEqual(len(fusion_results_dict['Splice Variants'][0]), 6)

	# def assertIsFile(self, path):
	# 	if not pl.Path(path).resolve().is_file():
	# 		raise AssertionError("File does not exist: %s" % str(path))

	def test_file_does_not_exist(self):
		path = pl.Path("this_file_does_not_exist.txt")
		self.assertEqual((str(path), path.is_file()), (str(path), False))

	def test_file_exists(self):
		path = pl.Path("testData/21M81042-1_CombinedVariantOutput.tsv")
		self.assertEqual((str(path), path.is_file()), (str(path), True))

	def test_remove_unwanted_line(self):
		fusion_results_dict = load_report('testData/21M81042-1_CombinedVariantOutput.tsv')
		remove_unwanted_line(fusion_results_dict)
		self.assertIsNot(fusion_results_dict['Fusions'][0][0].startswith('Gene'), True)
		self.assertIsNot(fusion_results_dict['Splice Variants'][0][0].startswith('Gene'), True)

	def test_remove_unwanted_line_fail(self):
		# Testing to see if unwanted line is still there and if so then it is passed in if/else statement
		fusion_results_dict = load_report('testData/21M81042-1_CombinedVariantOutput.tsv')
		self.assertIs(fusion_results_dict['Fusions'][0][0].startswith('Gene'), True)
		self.assertIs(fusion_results_dict['Splice Variants'][0][0].startswith('Gene'), True)

	def test_check_correct_df_is_made_from_format_fusions(self):
		fusion_results_dict = load_report('testData/21M81042-1_CombinedVariantOutput.tsv')
		fusion_results_dict_ntc = load_report('testData/NTC-RNA-1_CombinedVariantOutput.tsv')
		remove_unwanted_line(fusion_results_dict)
		df = format_fusions(fusion_results_dict, fusion_results_dict_ntc)
		self.assertEqual(len(df.keys()), 9)
		
	def test_correct_number_of_columns_for_final_merged_dataframe(self):
		all_fusions = 'testData/21M81042-1_AllFusions.csv'
		fusion_results_dict = load_report('testData/21M81042-1_CombinedVariantOutput.tsv')
		fusion_results_dict_ntc = load_report('testData/NTC-RNA-1_CombinedVariantOutput.tsv')
		sample_id = fusion_results_dict['Analysis Details'][0][1]
		remove_unwanted_line(fusion_results_dict)
		remove_unwanted_line(fusion_results_dict_ntc)
		df = format_fusions(fusion_results_dict, fusion_results_dict_ntc)
		df_extra = pd.read_csv(all_fusions, comment = '#')
		df_merged_final = add_extra_columns(df, df_extra)
		self.assertEqual(len(df_merged_final.keys()), 15)
		self.assertEqual(df_merged_final.type.nunique(), 2)


# class ActualTest(TestFusions):
	
# 	def test_file_exists(self):
# 		path = pl.Path("testData/21M81042-1_CombinedVariantOutput.tsv")
# 		self.assertIsFile(path)

# 	def test_file_does_not_exist(self):
# 		path = pl.Path("this_file_does_not_exist.txt")
# 		self.assertIsFile(path)



if __name__ == '__main__':
	unittest.main()

