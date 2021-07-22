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

worksheets=table["Sample_Plate"].unique()

#print a list of worksheets to a file
worksheets_file=open('worksheets.txt','w')

for worksheet in worksheets:
     worksheets_file.write(worksheet)
     worksheets_file.write('\n')
worksheets_file.close()

for worksheet in worksheets:
	table_worksheet=(table[table["Sample_Plate"]==worksheet])

	table_worksheet_RNA=(table_worksheet[table_worksheet["Sample_Type"]=="RNA"])
	table_worksheet_RNA.to_csv("samples_correct_order_"+worksheet+"_RNA.csv", index=False)


	table_worksheet_DNA=(table_worksheet[table_worksheet["Sample_Type"]=="DNA"])
	table_worksheet_DNA.to_csv("samples_correct_order_"+worksheet+"_DNA.csv", index=False)