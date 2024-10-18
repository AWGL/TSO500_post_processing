"""
Produce various sample/ worksheet lists for use downstream
include dictionary to translate referral types

"""

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

#Go through samplesheet until you hit the header lines
for line in samplesheet:
	
	#Remove new line character
	line = line.strip()

	#Skip if header line
	if line.startswith('Sample'):
		
		next

	else:
		
		#Split line into list
		line = line.split(",")

		#Get columns we need from sample sheet
		sample_id = line[0]
		worksheet = line[2]
		sample_type = line[7]
		description = line[9]

		#Append Sample ID (first element in list) to sample list
		samplelist.write(sample_id+"\n")

		#Get referral from Description (tenth element in list), split by ; and get third element
		referral = description.split(";")[2]
		referral = referral.split("=")[1]

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
			samplesaml = open('samples_aml_to_myeloid.csv",'a')
			
			samplesaml.write(sample_id+",myeloid\n")

			samplesaml.close()
			

#Write out worksheets to file
with open('worksheets_dna.txt','w') as f:
	for ws in dna:
		f.write(ws+"\n")

with open('worksheets_rna.txt','w') as f:
	for ws in rna:
		f.write(ws+"\n")

samplesheet.close()
samplelist.close()

