import unittest
import pathlib as pl
import pandas as pd

from fusions2db import load_report, remove_unwanted_line, format_fusions, add_extra_columns


"""

Test usage:

type 'python -m unittest tests/test_fusion_check.py' from within the directory that contains the fusions2db.py 

"""

class TestFusions(unittest.TestCase):

	# define filepaths
	test_data_location = 'test_data_fusion_check'
	test_sample_combined_file = f'{test_data_location}/sample_CombinedVariantOutput.tsv'
	test_sample_all_fusions_file = f'{test_data_location}/sample_AllFusions.csv'
	test_ntc_combined_file = f'{test_data_location}/NTC-RNA_CombinedVariantOutput.tsv'


	def test_load_report_with_sample(self):
		"""
		Checking that the correct number of columns are present within each dataframe

		"""
		fusion_results_dict = load_report(self.test_sample_combined_file)
		self.assertEqual(len(fusion_results_dict.keys()), 7)
		self.assertEqual(len(fusion_results_dict['Fusions'][0]), 6)
		self.assertEqual(len(fusion_results_dict['Splice Variants'][0]), 6)


	def test_file_does_not_exist(self):
		"""
		Checking that the correct message appears if file does not exist

		"""
		path = pl.Path("this_file_does_not_exist.txt")
		self.assertEqual((str(path), path.is_file()), (str(path), False))


	def test_file_exists(self):
		"""
		Checking that the if the file exists then the script is executed

		"""
		path = pl.Path(self.test_sample_combined_file)
		self.assertEqual((str(path), path.is_file()), (str(path), True))


	def test_remove_unwanted_line(self):
		"""
		Checking that the correct line is removed from the Fusions and Splice Variants lists

		"""
		fusion_results_dict = load_report(self.test_sample_combined_file)
		remove_unwanted_line(fusion_results_dict)
		self.assertIsNot(fusion_results_dict['Fusions'][0][0].startswith('Gene'), True)
		self.assertIsNot(fusion_results_dict['Splice Variants'][0][0].startswith('Gene'), True)


	def test_remove_unwanted_line_fail(self):
		"""
		Checking that the correct message appears if removing the unwanted line fails

		"""
		# Testing to see if unwanted line is still there and if so then it is passed in if/else statement
		fusion_results_dict = load_report(self.test_sample_combined_file)
		self.assertIs(fusion_results_dict['Fusions'][0][0].startswith('Gene'), True)
		self.assertIs(fusion_results_dict['Splice Variants'][0][0].startswith('Gene'), True)


	def test_check_correct_df_is_made_from_format_fusions(self):
		"""
		Checking to see if the correct dataframe is created by calculating the complete number of columns

		"""
		fusion_results_dict = load_report(self.test_sample_combined_file)
		fusion_results_dict_ntc = load_report(self.test_ntc_combined_file)
		remove_unwanted_line(fusion_results_dict)
		df = format_fusions(fusion_results_dict, fusion_results_dict_ntc)
		self.assertEqual(len(df.keys()), 9)


	def test_correct_number_of_columns_for_final_merged_dataframe(self):
		"""
		Checking that the numbers of columns is correct for the final merged dataframe

		"""
		all_fusions = self.test_sample_all_fusions_file
		fusion_results_dict = load_report(self.test_sample_combined_file)
		fusion_results_dict_ntc = load_report(self.test_ntc_combined_file)

		sample_id = fusion_results_dict['Analysis Details'][0][1]
		remove_unwanted_line(fusion_results_dict)
		remove_unwanted_line(fusion_results_dict_ntc)
		df = format_fusions(fusion_results_dict, fusion_results_dict_ntc)
		df_extra = pd.read_csv(all_fusions, comment = '#')
		df_merged_final = add_extra_columns(df, df_extra)
		self.assertEqual(len(df_merged_final.keys()), 15)
		self.assertEqual(df_merged_final.type.nunique(), 2)


	def test_correct_ntc_result_is_being_displayed(self):
		"""
		Checking to see if in_ntc column is providing correct result

		"""
		all_fusions = self.test_sample_all_fusions_file
		fusion_results_dict = load_report(self.test_sample_combined_file)
		fusion_results_dict_ntc = load_report(self.test_ntc_combined_file)

		sample_id = fusion_results_dict['Analysis Details'][0][1]
		remove_unwanted_line(fusion_results_dict)
		remove_unwanted_line(fusion_results_dict_ntc)
		df = format_fusions(fusion_results_dict, fusion_results_dict_ntc)
		df_extra = pd.read_csv(all_fusions, comment = '#')
		df_merged_final = add_extra_columns(df, df_extra)
		self.assertEqual(df_merged_final['in_ntc'][0] == True, True)


	def test_correct_data_is_being_displayed(self):
		"""
		Checking to see if the correct data is being pulled through

		"""
		all_fusions = self.test_sample_all_fusions_file
		fusion_results_dict = load_report(self.test_sample_combined_file)
		fusion_results_dict_ntc = load_report(self.test_ntc_combined_file)

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
