"""
Produce various sample/ worksheet lists for use downstream
include dictionary to translate referral types

"""
import csv

#Open sample sheet
samplesheet = open('SampleSheet_updated.csv','r')

#Open sample_list.txt to write to
samplelist = open('sample_list.txt','w')


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

#Sets for worksheet ids
dna = set()
rna = set()

#Get column by header name instead of header position.
read_samplesheet = csv.DictReader(samplesheet)

#Go through samplesheet until you hit the header lines
for line in read_samplesheet:
	
	#Get columns we need from sample sheet
	sample_id = line["Sample_ID"]
	worksheet = line["Sample_Plate"]
	sample_type = line["Sample_Type"]
	description = line["Description"]

	#Append Sample ID (first element in list) to sample list
	samplelist.write(sample_id+"\n")

	# Split the decription column up, as now additional referrals section
	desc_parts = description.split(";")
	desc_dict = {}
	for part in desc_parts:
		if "=" in part:
			key, value = part.split("=", 1)
			desc_dict[key] = value

	#Get referral from Description (tenth element in list), split by ; and get third element
	referral = desc_dict.get("referral", "null")

	#if RNA, update referral based on dictionary
	if sample_type == "RNA" and (referral in referral_dict):

		referral = referral_dict[referral]
	
	#Add worksheet to set
	if sample_type == "DNA":
		dna.add(worksheet)

	elif sample_type == "RNA":
		rna.add(worksheet)

	#Write to samples correct order
	samplescorrect = open('samples_correct_order_'+worksheet+"_"+sample_type+".csv",'a') 

	samplescorrect.write(sample_id+","+worksheet+","+sample_type+","+referral+"\n")
	
	samplescorrect.close()

	#Write any aml referral samples to additional csv
	if referral == "aml":
		samplesaml = open("samples_aml_to_myeloid_"+worksheet+"_"+sample_type+".csv",'a')
		
		samplesaml.write(sample_id+",myeloid\n")

		samplesaml.close()

	# Get additional referrals to make one csv per referral
	additional_referrals = desc_dict.get("additional_referrals", "")
	if additional_referrals:
		for add_ref in additional_referrals.split(","):
			add_ref = add_ref.strip()
			if add_ref:
				add_ref_file = open(
					"samples_additional_referral_"+add_ref+"_"+worksheet+"_"+sample_type+".csv", 'a'	
				)
				add_ref_file.write(sample_id+","+add_ref+"\n")
				add_ref_file.close()

#Write out worksheets to file
with open('worksheets_dna.txt','w') as f:
	for ws in dna:
		f.write(ws+"\n")

with open('worksheets_rna.txt','w') as f:
	for ws in rna:
		f.write(ws+"\n")

samplesheet.close()
samplelist.close()

