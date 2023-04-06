"""
Produce various sample/ worksheet lists for use downstream
include dictionary to translate referral types

"""

#Open sample sheet
SAMPLESHEET = open('SampleSheet_updated.csv','r')

#Open sample_list.txt to write to
SAMPLELIST = open('sample_list.txt','w')


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
for line in SAMPLESHEET:
	
	#Remove new line character
	line = line.strip()

	#Skip if header line
	if line.startswith('Sample'):
		
		next

	else:
		
		#Split line into list
		line = line.split(",")

		#Append Sample ID (first element in list) to sample list
		SAMPLELIST.write(line[0]+"\n")

		#Get referral from Description (tenth element in list), split by ; and get third element
		referral = line[9].split(";")[2]
		referral = referral.split("=")[1]

		#if RNA, update referral based on dictionary
		if line[7]== "RNA" and (referral in referral_dict):
	
			referral = referral_dict[referral]
		
		#Add worksheet to set
		if line[7] == "DNA":
			dna.add(line[2])

		elif line[7] == "RNA":
			rna.add(line[2])

		#Write to samples correct order
		SAMPLESCORRECT = open('samples_correct_order_'+line[2]+"_"+line[7]+".csv",'a')
		
		SAMPLESCORRECT.write(line[0]+","+line[2]+","+line[7]+","+referral+"\n")
		
		SAMPLESCORRECT.close()

#Write out worksheets to file
with open('worksheets_dna.txt','w') as f:
	for ws in dna:
		f.write(ws+"\n")

with open('worksheets_rna.txt','w') as f:
	for ws in rna:
		f.write(ws+"\n")

SAMPLESHEET.close()
SAMPLELIST.close()

