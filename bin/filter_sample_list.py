#!/usr/bin/env -S uv run --script

"""
Produce various sample/ worksheet lists for use downstream
include dictionary to translate referral types

"""
import csv
import io

def parse_sample_sheet(input_file):
	# Read the sample sheet file
	with open(input_file) as f:
		lines = f.readlines()

	# Find the index of the line that starts with "[Data]"
	data_index = next(i for i, line in enumerate(lines) if line.startswith("[Data]"))

	# Extract the lines after the "[Data]" line
	data_lines = lines[data_index + 1:]

	# Use csv.DictReader to parse the data lines into a list of dictionaries
	reader = csv.DictReader(io.StringIO("".join(data_lines)))

	return [row for row in reader]


def parse_description(description):
    return dict(pair.split("=") for pair in description.split(";"))


# create referral dictionary - this is only for legacy panels, any new panels should be lower case
referral_dict  = {
    "colorectal": "Colorectal",
    "gist": "GIST",
    "glioma": "Glioma",
    "lung": "Lung",
    "melanoma": "Melanoma",
    "thyroid": "Thyroid",
    "tumour": "Tumour",
    "ntrk": "NTRK",
    "null": "null",
}

# Sets for worksheet IDs and sample list
dna_worksheets = set()
rna_worksheets = set()
samples_list = set()

# Iterate through samples
samples = parse_sample_sheet("SampleSheet.csv")
for line in samples:

	# Get columns we need from sample sheet
	sample_id = line['Sample_ID']
	worksheet = line['Sample_Plate']
	sample_type = line['Sample_Type']
	description = parse_description(line['Description'])

	# Append Sample ID (first element in list) to sample list
	samples_list.add(sample_id)

	# Get referral from Description (tenth element in list), split by ; and get third element
	referral = description['referral']

	# if RNA, update referral based on dictionary
	if sample_type == "RNA" and (referral in referral_dict):
		referral = referral_dict[referral]

	# Add worksheet to set
	if sample_type == "DNA":
		dna_worksheets.add(worksheet)

	elif sample_type == "RNA":
		rna_worksheets.add(worksheet)

	# Write to samples correct order
	with open(f'samples_correct_order_{worksheet}_{sample_type}.csv','a') as samples_correct:
		samples_correct.write(f'{sample_id},{worksheet},{sample_type},{referral}\n')

	# Write any aml referral samples to additional csv
	if referral == "aml":
		with open(f"samples_aml_to_myeloid_{worksheet}_{sample_type}.csv",'a') as samples_aml:
			samples_aml.write(f"{sample_id},myeloid\n")
			

# Write out worksheets to file
with open('worksheets_dna.txt','w') as f:
	for ws in dna_worksheets:
		f.write(ws+"\n")

with open('worksheets_rna.txt','w') as f:
	for ws in rna_worksheets:
		f.write(ws+"\n")

# Write out sample list to file
with open('samples.txt','w') as f:
	for sample in samples_list:
		f.write(sample+"\n")

# Write out a csv of the filtered samplesheet to file
with open('SampleSheet_updated.csv','w') as f:
	writer = csv.writer(f)
	# Write header
	writer.writerow(samples[0].keys())
	# Write rows
	for line in samples:
		writer.writerow(line.values())
