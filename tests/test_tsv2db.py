import unittest
import argparse
import decimal
from decimal import Decimal

from tsv2db import parse_ntc_tsv, parse_sample_tsv


"""
test usage:

type 'python -m unittest' from within a directory which contains the test_tsv2db.py script and test_data folder.
"""



class test_tsv2db(unittest.TestCase):

	def test_parse_ntc_tsv_novar(self):
		"""
		test read of ntc tsv with no variants
		"""
		ntc_var_list = parse_ntc_tsv('test_data/NTC_no_var.tsv')
		self.assertEqual(ntc_var_list, ['NA'])


	def test_parse_ntc_tsv_var(self):
		"""
		test read of first 7 columns of variant info correct
		"""
		ntc_var_list = parse_ntc_tsv('test_data/NTC_fake_variants.tsv')
		self.assertEqual(ntc_var_list[0], ['TNFRSF14', 'chr1', '2489805', 'G', 'A', '0.0936', '87'])

		## test read of first 7 columns of cariant info is correct if gene name is missing
		self.assertEqual(ntc_var_list[3][0], '')


	def test_parse_ntc_tsv_novar_blanks(self):
		"""
		test "NA" still output if there are extra blanks left after the first line
		"""
		ntc_var_list = parse_ntc_tsv('test_data/NTC_novar_blanks.tsv')
		self.assertEqual(ntc_var_list, ['NA'])


	def test_parse_sample_tsv_no_var(self):
		"""
		test that when no variants in the sample then output is empty list
		"""
		ntc_var_list = parse_ntc_tsv('test_data/NTC_no_var.tsv')
		sample_var_list = parse_sample_tsv('test_data/12345_CombinedVariantOutput_novar.tsv',ntc_var_list)
		self.assertEqual(sample_var_list, [])

	def test_parse_sample_tsv_no_ntc_var(self):
		"""
		test that when no variants in the NTC then all "in_ntc" column is "False"
		"""
		ntc_var_list = parse_ntc_tsv('test_data/NTC_no_var.tsv')
		sample_var_list = parse_sample_tsv('test_data/12345_CombinedVariantOutput.tsv',ntc_var_list)
		for item in sample_var_list:
			with self.subTest(item = item):
				self.assertEqual(item.split('\t')[12], 'False')


	def test_parse_sample_tsv_var_ntc_fake(self):
		"""
		test that if ntc has a variant then "in_ntc" column is changed to True and 'ntc_vaf', 'ntc_depth', 'ntc_alt_reads' are populated
		"""
		ntc_var_list = parse_ntc_tsv('test_data/NTC_fake_variants.tsv')
		sample_var_list = parse_sample_tsv('test_data/12345_CombinedVariantOutput.tsv',ntc_var_list)
		test_list = sample_var_list[1].split('\t')
		self.assertEqual(test_list[12], 'True')
		self.assertNotEqual(test_list[13], '')
		self.assertNotEqual(test_list[14], '')
		self.assertNotEqual(test_list[15], '')


	def test_parse_sample_tsv_extra_blanks(self):
		"""
		test that extra lines at end of sample tsv file doesnt error script and doesnt put out blanks in the sample var list
		"""
		ntc_var_list = parse_ntc_tsv('test_data/NTC_no_var.tsv')
		sample_var_list = parse_sample_tsv('test_data/12345_CombinedVariantOutput_blanks.tsv',ntc_var_list)
		self.assertNotEqual(sample_var_list[-1], '')





if __name__ == '__main__':
	unittest.main()