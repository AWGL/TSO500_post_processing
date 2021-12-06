#!/usr/bin/env python

# from fusions_check_with_ntc_extra_columns import *
from fusions_check_with_ntc import *
#load_report, remove_unwanted_line, format_fusions, create_final_df, add_extra_columns

import unittest
import pathlib as pl


"""

Test usage:

type 'python -m unittest tests/test_fusion_check.py' from within the directory that contains the fusions_check_with_ntc.py 

"""



class TestFusions(unittest.TestCase):

	# Checking that the correct number of columns are present within each dataframe
	def test_load_report_with_sample(self):
		fusion_results_dict = load_report('tests/test_data_fusion_check/21M81042-1_CombinedVariantOutput.tsv')
		self.assertEqual(len(fusion_results_dict.keys()), 7)
		self.assertEqual(len(fusion_results_dict['Fusions'][0]), 6)
		self.assertEqual(len(fusion_results_dict['Splice Variants'][0]), 6)

	# Checking that the correct message appears if file does not exist
	def test_file_does_not_exist(self):
		path = pl.Path("this_file_does_not_exist.txt")
		self.assertEqual((str(path), path.is_file()), (str(path), False))

	# Checking that the if the file exists then the script is executed
	def test_file_exists(self):
		path = pl.Path("tests/test_data_fusion_check/21M81042-1_CombinedVariantOutput.tsv")
		self.assertEqual((str(path), path.is_file()), (str(path), True))

	# Checking that the correct line is removed from the Fusions and Splice Variants lists
	def test_remove_unwanted_line(self):
		fusion_results_dict = load_report('tests/test_data_fusion_check/21M81042-1_CombinedVariantOutput.tsv')
		remove_unwanted_line(fusion_results_dict)
		self.assertIsNot(fusion_results_dict['Fusions'][0][0].startswith('Gene'), True)
		self.assertIsNot(fusion_results_dict['Splice Variants'][0][0].startswith('Gene'), True)

	# Checking that the correct message appears if removing the unwanted line fails
	def test_remove_unwanted_line_fail(self):
		# Testing to see if unwanted line is still there and if so then it is passed in if/else statement
		fusion_results_dict = load_report('tests/test_data_fusion_check/21M81042-1_CombinedVariantOutput.tsv')
		self.assertIs(fusion_results_dict['Fusions'][0][0].startswith('Gene'), True)
		self.assertIs(fusion_results_dict['Splice Variants'][0][0].startswith('Gene'), True)

	# Checking to see if the correct dataframe is created by calculating the complete number of columns
	def test_check_correct_df_is_made_from_format_fusions(self):
		fusion_results_dict = load_report('tests/test_data_fusion_check/21M81042-1_CombinedVariantOutput.tsv')
		fusion_results_dict_ntc = load_report('tests/test_data_fusion_check/NTC-RNA-1_CombinedVariantOutput.tsv')
		remove_unwanted_line(fusion_results_dict)
		df = format_fusions(fusion_results_dict, fusion_results_dict_ntc)
		self.assertEqual(len(df.keys()), 9)

	# Checking that the numbers of columns is correct for the final merged dataframe
	def test_correct_number_of_columns_for_final_merged_dataframe(self):
		all_fusions = 'tests/test_data_fusion_check/21M81042-1_AllFusions.csv'
		fusion_results_dict = load_report('tests/test_data_fusion_check/21M81042-1_CombinedVariantOutput.tsv')
		fusion_results_dict_ntc = load_report('tests/test_data_fusion_check/NTC-RNA-1_CombinedVariantOutput.tsv')

		sample_id = fusion_results_dict['Analysis Details'][0][1]
		remove_unwanted_line(fusion_results_dict)
		remove_unwanted_line(fusion_results_dict_ntc)
		df = format_fusions(fusion_results_dict, fusion_results_dict_ntc)
		df_extra = pd.read_csv(all_fusions, comment = '#')
		df_merged_final = add_extra_columns(df, df_extra)
		self.assertEqual(len(df_merged_final.keys()), 15)
		self.assertEqual(df_merged_final.type.nunique(), 2)

	# Checking to see if in_ntc column is providing correct result
	def test_correct_ntc_result_is_being_displayed(self):
		all_fusions = 'tests/test_data_fusion_check/21M81042-1_AllFusions.csv'
		fusion_results_dict = load_report('tests/test_data_fusion_check/21M81042-1_CombinedVariantOutput.tsv')
		fusion_results_dict_ntc = load_report('tests/test_data_fusion_check/NTC-RNA-1_CombinedVariantOutput.tsv')

		sample_id = fusion_results_dict['Analysis Details'][0][1]
		remove_unwanted_line(fusion_results_dict)
		remove_unwanted_line(fusion_results_dict_ntc)
		df = format_fusions(fusion_results_dict, fusion_results_dict_ntc)
		df_extra = pd.read_csv(all_fusions, comment = '#')
		df_merged_final = add_extra_columns(df, df_extra)
		self.assertEqual(df_merged_final['in_ntc'][0] == True, True)

		# Checking to see if the correct data is being pulled through
	def test_correct_data_is_being_displayed(self):
		all_fusions = 'tests/test_data_fusion_check/21M81042-1_AllFusions.csv'
		fusion_results_dict = load_report('tests/test_data_fusion_check/21M81042-1_CombinedVariantOutput.tsv')
		fusion_results_dict_ntc = load_report('tests/test_data_fusion_check/NTC-RNA-1_CombinedVariantOutput.tsv')

		sample_id = fusion_results_dict['Analysis Details'][0][1]
		remove_unwanted_line(fusion_results_dict)
		remove_unwanted_line(fusion_results_dict_ntc)

		df = format_fusions(fusion_results_dict, fusion_results_dict_ntc)
		df_extra = pd.read_csv(all_fusions, comment = '#')
		df_merged_final = add_extra_columns(df, df_extra)
		# Checking the final result against randomly picked line in original tsv file to check the correct data is being pulled through for each column
		self.assertEqual(df_merged_final['fusion'][0] == fusion_results_dict['Fusions'][0][0], True)
		self.assertEqual(df_merged_final['fusion'][9] == fusion_results_dict['Fusions'][9][0], True)
		self.assertEqual(df_merged_final['exons'][16] == fusion_results_dict['Splice Variants'][0][1], True)
		self.assertEqual(df_merged_final['reference_reads_1'][3] == fusion_results_dict['Fusions'][3][4], True)
		self.assertEqual(df_merged_final['left_breakpoint'][10] == fusion_results_dict['Fusions'][10][1], True)
		self.assertEqual(df_merged_final['right_breakpoint'][15] == fusion_results_dict['Fusions'][15][2], True)


if __name__ == '__main__':
	unittest.main()

