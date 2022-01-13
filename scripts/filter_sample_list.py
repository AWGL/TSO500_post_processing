"""
Produce various sample/ worksheet lists for use downstream

"""

import pandas

# load in truncated samplesheet
table = pandas.read_csv('SampleSheet_updated.csv', sep=',')

# save sample list
sample_list = table['Sample_ID']
sample_list.to_csv('sample_list.txt', index=False, header=False)

# filter table to only include columns below
table = table[['Sample_ID', 'Sample_Plate','Sample_Type', 'Description']]

# parse referral from description field
table['Referral_value'] = table['Description'].str.split(';', expand=True)[2]
table['Referral'] = table.Referral_value.str.split('=',expand=True)[1]
table = table.drop(['Referral_value'], axis=1)
table = table.drop(['Description'], axis=1)

# subset all DNA samples
table_dna = table[table['Sample_Type'] == 'DNA']
worksheets_dna = table_dna['Sample_Plate'].unique()

# save a list of DNA worksheets
with open('worksheets_dna.txt','w') as worksheets_file_dna:
	for worksheet in worksheets_dna:
		worksheets_file_dna.write(worksheet)
		worksheets_file_dna.write('\n')

# save one list per worksheet with all samples plus worksheet and referral info
for worksheet in worksheets_dna:
	table_dna_worksheet = table_dna
	table_dna_worksheet = table_dna[table_dna['Sample_Plate'] == worksheet]
	table_dna_worksheet.to_csv(f'samples_correct_order_{worksheet}_DNA.csv', index=False, header=False)

# subset all RNA samples
table_rna = table[table['Sample_Type'] == 'RNA']
worksheets_rna = table_rna['Sample_Plate'].unique()

# save a list of RNA worksheets
with open('worksheets_rna.txt','w') as worksheets_file_rna:
	for worksheet in worksheets_rna:
		worksheets_file_rna.write(worksheet)
		worksheets_file_rna.write('\n')

# save one list per worksheet with all samples plus worksheet and referral info
for worksheet in worksheets_rna:
	table_rna_worksheet = table_rna
	table_rna_worksheet = table_rna[table_rna['Sample_Plate'] == worksheet]
	table_rna_worksheet.to_csv(f'samples_correct_order_{worksheet}_RNA.csv', index=False, header=False)
