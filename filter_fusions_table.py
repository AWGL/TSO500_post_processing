import sys
import pandas

sample=sys.argv[1]
version=sys.argv[2]

#get list of genes in tumour referral
tumour_referral_file=pandas.read_csv("/data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-"+version+"/RNA_referrals/Tumour.txt", sep="\t")
gene_list=list(tumour_referral_file['Genes'])

#get fusion table for sample
fusions_table=pandas.read_csv("./Gathered_Results/Database/"+sample+"_fusion_check.csv", sep=",")

#Create a fusions file for every gene on the tumour panel


for gene in gene_list:


        fusions_table_gene=fusions_table[fusions_table["fusion"].str.contains(gene)]

        fusions_table_gene.to_csv("./Gathered_Results/Results/"+sample+"/"+gene+"_fusions.csv", index=False)


