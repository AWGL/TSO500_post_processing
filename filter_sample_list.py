import pandas

table=pandas.read_csv("SampleSheet_updated.csv", sep=",")
sample_list=table['Sample_ID']
sample_list.to_csv("sample_list.txt", index=False, header=False)


print(table)
table=table[["Sample_ID", "Sample_Plate","Sample_Type", "Description"]]

table["Referral_value"]=table["Description"].str.split(';', expand=True)[2]


table['Referral'] = table.Referral_value.str.split("=",expand=True)[1]
print(table)
table=table.drop(["Referral_value"], axis=1)
table=table.drop(["Description"], axis=1)

table_rna=(table[table["Sample_Type"]=="RNA"])
table_dna=(table[table["Sample_Type"]=="DNA"])


worksheets_rna=table_rna["Sample_Plate"].unique()
worksheets_dna=table_dna["Sample_Plate"].unique()


#print a list of worksheets to a file
worksheets_file_dna=open('worksheets_dna.txt','w')

for worksheet in worksheets_dna:
     worksheets_file_dna.write(worksheet)
     worksheets_file_dna.write('\n')
worksheets_file_dna.close()


#print a list of worksheets to a file
worksheets_file_rna=open('worksheets_rna.txt','w')

for worksheet in worksheets_rna:
     worksheets_file_rna.write(worksheet)
     worksheets_file_rna.write('\n')
worksheets_file_rna.close()


for worksheet in worksheets_rna:
	table_rna_worksheet=table_rna
	table_rna_worksheet=(table_rna[table_rna["Sample_Plate"]==worksheet])
	table_rna_worksheet.to_csv("samples_correct_order_"+worksheet+"_RNA.csv", index=False, header=False)

for worksheet in worksheets_dna:
	table_dna_worksheet=table_dna
	table_dna_worksheet=(table_dna[table_dna["Sample_Plate"]==worksheet])
	table_dna_worksheet.to_csv("samples_correct_order_"+worksheet+"_DNA.csv", index=False, header=False)