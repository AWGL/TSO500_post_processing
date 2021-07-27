import sys
import pandas
from collections import OrderedDict

worksheetid=sys.argv[1]
version=sys.argv[2]

referrals_path="/data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-"+version+"/RNA_referrals"


#need to determine if the samples are RNA or need samples_correct_order_file- with second column for dna or rna
samples_file = pandas.read_csv("./samples_correct_order_"+worksheetid+"_RNA.csv", sep=",", names=["Sample", "Worksheet", "Type", "Referral"])


sampleList=samples_file["Sample"].tolist()
len_sample_list=len(sampleList)

#contamination dictionary
contamination=["No"]*len_sample_list
contamination_dict=OrderedDict(zip(sampleList, contamination))

#contamination referral dictionary
contamination_referral=["No"]*len_sample_list
contamination_referral_dict=OrderedDict(zip(sampleList, contamination_referral))

#get a list of genes in the tumour referral type to be able to later filter by genes on the panel
tumour_referral_file=pandas.read_csv("/data/diagnostics/pipelines/TSO500_RUO_LocalApp/TSO500_RUO_LocalApp-"+version+"/RNA_referrals/Tumour.txt", sep="\t")
gene_list=list(tumour_referral_file['Genes'])
tumour_referrals="|".join(gene_list)


sample_number=0

for sample in sampleList:
    print(sample)
    NTC_in_sample= ("NTC" in sample.upper())

    if (NTC_in_sample==False):

        sample_table=samples_file[samples_file['Sample']==sample]
        referral_table=sample_table.filter(items=['Referral'])
        referral=referral_table.iloc[0,0]

	# read in fusions file and duplicate the gene for splice variants in the fusions column
        all_fusions_table= pandas.read_csv("./Gathered_Results/Database/"+sample+"_fusion_check.csv", sep=",")
        splice_fusion=all_fusions_table[(all_fusions_table.type=="Splice")]
        if (len(splice_fusion)!=0):
            splice_fusion['fusion'] = splice_fusion['fusion'].apply(lambda x: x+"-"+x)
            fusion=all_fusions_table[(all_fusions_table.type=="Fusion")]
            report=splice_fusion.append(fusion)
        else:
            report=all_fusions_table

        #get list of genes in the referral for the sample to be used to filter file with all fusions by referral later
        referral_file=pandas.read_csv(referrals_path+"/"+referral+".txt", sep="\t")
        gene_list=list(referral_file['Genes'])
        len_gene_list=len(gene_list)

        #set default contamination to "No"
        contamination="No"
        contamination_next="No"
        contamination_previous="No"
        contamination_referral="No"
        contamination_referral_next="No"
        contamination_referral_previous="No"

        #get the sampleId of the sample before and after in the sampleList
        if (sample_number!=0):
            sample_previous_number=sample_number-1
            sample_previous=sampleList[sample_previous_number]

	    # read in fusions file for previous sample and duplicate the gene for splice variants in the fusions column
            all_fusions_table_previous= pandas.read_csv("./Gathered_Results/Database/"+sample_previous+"_fusion_check.csv", sep=",")
            splice_fusion_previous=all_fusions_table_previous[(all_fusions_table_previous.type=="Splice")]
            if ((len(splice_fusion_previous))!=0):
                splice_fusion_previous['fusion'] = splice_fusion_previous['fusion'].apply(lambda x: x+"-"+x)
                fusion_previous=all_fusions_table_previous[(all_fusions_table_previous.type=="Fusion")]
                report_previous=splice_fusion_previous.append(fusion_previous)
            else:
                report_previous=all_fusions_table_previous

 
        if (sample_number!=(len_sample_list-1)):
            sample_next_number=sample_number+1
            sample_next=sampleList[sample_next_number]

	    # read in fusions file for next sample and duplicate the gene for splice variants in the fusions column
            all_fusions_table_next= pandas.read_csv("./Gathered_Results/Database/"+sample_next+"_fusion_check.csv", sep=",")
            splice_fusion_next=all_fusions_table_next[(all_fusions_table_next.type=="Splice")]
            if ((len(splice_fusion_next))!=0):
                splice_fusion_next['fusion'] = splice_fusion_next['fusion'].apply(lambda x: x+"-"+x)
                fusion_next=all_fusions_table_next[(all_fusions_table_next.type=="Fusion")]
                report_next=splice_fusion_next.append(fusion_next)
            else:
                report_next=all_fusions_table_next


        fusion_list=[]
        fusion_list_previous=[]
        fusion_list_next=[]
        report=report[report["fusion"].str.contains(tumour_referrals)]
        if (len(report)>0):
            fusion_list1=report["fusion"].tolist()

            #get a list of fusions with the genes the alternate way round
            report["fusion"]=report["fusion"].astype(str)
            report["gene1"]=report["fusion"].str.split('-|/', expand=True)[0]
            report["gene2"]=report["fusion"].str.split('-|/', expand=True)[1]
            report["fusion2"]=report["gene2"].str.cat(report["gene1"], sep="-")
            fusion_list2=report["fusion2"].tolist()

            #create list containing the fusions with genes in both orders
            fusion_list=fusion_list1+fusion_list2


            #compare the list of fusions in the current sample with those in the sample before
            if (sample_number!=0):
                report_previous=report_previous[report_previous["fusion"].str.contains(tumour_referrals)]
                if (len(report_previous)>0):
                    
                    #TODO- filter by the tumour config file instead of making it hard coded 
                    report_previous=report_previous[report_previous["fusion"].str.contains(tumour_referrals)]
                    if (len(report_previous)>0):
                        fusion_list_previous1=report_previous["fusion"].tolist()

                        #get a list of fusions with the genes the alternate way round for the sample before
                        report_previous["fusion"]=report_previous["fusion"].astype(str)
                        report_previous["gene1"]=report_previous["fusion"].str.split('-|/', expand=True)[0]
                        report_previous["gene2"]=report_previous["fusion"].str.split('-|/', expand=True)[1]
                        report_previous["fusion2"]=report_previous["gene2"].str.cat(report_previous["gene1"], sep="-")
                        fusion_list_previous2=report_previous["fusion2"].tolist()

                        #create list containing the fusions with genes in both orders for sample before
                        fusion_list_previous=fusion_list_previous1+fusion_list_previous2

                        #compare the fusion list with the fusion list for the sample before. Label contamination as "yes" if any fusions match.                            
                        for fusion_current in fusion_list:
                            for fusion_previous in fusion_list_previous:
                                if (fusion_current==fusion_previous):
                                    contamination_previous="Yes"
                                    contamination="Yes"

            #compare the list of fusions in the current sample with those in the sample after
            if (sample_number<(len_sample_list-1)):
                report_next=report_next[report_next["fusion"].str.contains(tumour_referrals)]
                if (len(report_next)>0):
                    report_next=report_next[report_next["fusion"].str.contains(tumour_referrals)]
                    if (len(report_next)>0):
                        fusion_list_next1=report_next["fusion"].tolist()

                        #get a list of fusions with the genes the alternate way round for the sample after 
                        report_next["fusion"]=report_next["fusion"].astype(str)
                        report_next["gene1"]=report_next["fusion"].str.split('-|/', expand=True)[0]
                        report_next["gene2"]=report_next["fusion"].str.split('-|/', expand=True)[1]
                        report_next["fusion2"]=report_next["gene2"].str.cat(report_next["gene1"], sep="-")
                        fusion_list_next2=report_next["fusion2"].tolist()

                        #create list containing the fusions with genes in both orders for the sample after
                        fusion_list_next=fusion_list_next1+fusion_list_next2

                        #compare the fusion list for the sample with the fusion list for the sample after. Label contamination as "yes" if any fusions match.
                        for fusion_current in fusion_list:
                            for fusion_next in fusion_list_next:
                                if (fusion_current==fusion_next):
                                    contamination_next="Yes"
                                    contamination="Yes"


        #compare the fusions within the sample to the fusions in sample before and after for the REFERRAL

        #loop through all the genes in the referral
        if (len_gene_list!=0):

            for gene in gene_list:


	        # read in fusions file for sample referral and duplicate the gene for splice variants in the fusions column

                all_fusions_table_referral= pandas.read_csv("./Gathered_Results/Database/"+sample+"_fusion_check.csv", sep=",")
                all_fusions_table_referral=all_fusions_table_referral[all_fusions_table_referral["fusion"].str.contains(gene)]

                splice_fusion_referral=all_fusions_table_referral[(all_fusions_table_referral.type=="Splice")]
                if ((len(splice_fusion_referral))!=0):
                    splice_fusion_referral['fusion'] = splice_fusion_referral['fusion'].apply(lambda x: x+"-"+x)
                    fusion_referral=all_fusions_table_referral[(all_fusions_table_referral.type=="Fusion")]
                    report_referral=splice_fusion.append(fusion_referral)
                else:
                    report_referral=all_fusions_table_referral

                if (len(report_referral)>0):
                    fusion_list_referral1=report_referral["fusion"].tolist()

                    #get a list of fusions with the genes the alternate way round
                    report_referral["fusion"]=report_referral["fusion"].astype(str)
                    report_referral["gene1"]=report_referral["fusion"].str.split('-|/', expand=True)[0]
                    report_referral["gene2"]=report_referral["fusion"].str.split('-|/', expand=True)[1]
                    report_referral["fusion2"]=report_referral["gene2"].str.cat(report_referral["gene1"], sep="-")
                    fusion_list_referral2=report_referral["fusion2"].tolist()

                    #create list containing the fusions with genes in both orders for sample
                    fusion_list_referral=fusion_list_referral1+fusion_list_referral2

        
                    if (sample_number!=0): 
                        print(sample_previous)

	                # read in fusions file for sample referral and duplicate the gene for splice variants in the fusions column
                        all_fusions_table_referral_previous= pandas.read_csv("./Gathered_Results/Database/"+sample_previous+"_fusion_check.csv", sep=",")
                        all_fusions_table_referral_previous=all_fusions_table_referral_previous[all_fusions_table_referral_previous["fusion"].str.contains(gene)]
                        splice_fusion_referral_previous=all_fusions_table_referral_previous[(all_fusions_table_referral_previous.type=="Splice")]
                        if ((len(splice_fusion_referral_previous))!=0):
                            splice_fusion_referral_previous['fusion'] = splice_fusion_referral_previous['fusion'].apply(lambda x: x+"-"+x)
                            fusion_referral_previous=all_fusions_table_referral_previous[(all_fusions_table_referral_previous.type=="Fusion")]
                            report_referral_previous=splice_fusion.append(fusion_referral_previous)
                        else:
                            report_referral_previous= all_fusions_table_referral_previous

                        if (len(report_referral_previous)>0):
                                fusion_list_referral_previous1=report_referral_previous["fusion"].tolist()

                                #get a list of fusions with the genes the alternate way round for the sample before
                                report_referral_previous["fusion"]=report_referral_previous["fusion"].astype(str)
                                report_referral_previous["gene1"]=report_referral_previous["fusion"].str.split('-|/', expand=True)[0]
                                report_referral_previous["gene2"]=report_referral_previous["fusion"].str.split('-|/', expand=True)[1]
                                report_referral_previous["fusion2"]=report_referral_previous["gene2"].str.cat(report_referral_previous["gene1"], sep="-")
                                fusion_list_referral_previous2=report_referral_previous["fusion2"].tolist()

 
                                #create list containing the fusions with genes in both orders for sample before
                                fusion_list_referral_previous=fusion_list_referral_previous1+fusion_list_referral_previous2

                                #compare the fusion list for the sample with fusion list for the sample before. Label contamination as "yes" if any fusions match.
                                for fusion_current_referral in fusion_list_referral:
                                    for fusion_previous_referral in fusion_list_referral_previous:
                                        if (fusion_current_referral==fusion_previous_referral):
                                            contamination_referral="Yes"

                    if (sample_number<(len_sample_list-1)):              

	                # read in fusions file for sample referral and duplicate the gene for splice variants in the fusions column
                        all_fusions_table_referral_next= pandas.read_csv("./Gathered_Results/Database/"+sample_next+"_fusion_check.csv", sep=",")
                        all_fusions_table_referral_next=all_fusions_table_referral_next[all_fusions_table_referral_next["fusion"].str.contains(gene)]

                        splice_fusion_referral_next=all_fusions_table_referral_next[(all_fusions_table_referral_next.type=="Splice")]
                        if ((len(splice_fusion_referral_next))!=0):
                            splice_fusion_referral_next['fusion'] = splice_fusion_referral_next['fusion'].apply(lambda x: x+"-"+x)
                            fusion_referral_next=all_fusions_table_referral_next[(all_fusions_table_referral_next.type=="Fusion")]
                            report_referral_next=splice_fusion.append(fusion_referral_next)                
                        else:
                            report_referral_next=all_fusions_table_referral_next


                        if (len(report_referral_next)>0):
                                fusion_list_referral_next1=report_referral_next["fusion"].tolist()

                                #get a list of fusions with the genes the alternate way round for the sample after
                                report_referral_next["fusion"]=report_referral_next["fusion"].astype(str)
                                report_referral_next["gene1"]=report_referral_next["fusion"].str.split('-|/', expand=True)[0]
                                report_referral_next["gene2"]=report_referral_next["fusion"].str.split('-|/', expand=True)[1]
                                report_referral_next["fusion2"]=report_referral_next["gene2"].str.cat(report_referral_next["gene1"], sep="--")
                                fusion_list_referral_next2=report_referral_next["fusion2"].tolist()

                                #create list containing the fusions with genes in both orders for sample after
                                fusion_list_referral_next=fusion_list_referral_next1+fusion_list_referral_next2

                                #compare the fusion list for the sample with the fusion list for the sample after. Label contamination as "yes" if any fusions match.
                                for fusion_current_referral in fusion_list_referral_next:
                                    for fusion_next_referral in fusion_list_referral_next:
                                        if (fusion_next_referral==fusion_next_referral):
                                            contamination_referral="Yes"


        #Replace the contamination values in the contamination dictionaries

        if (contamination_dict[sample]=="No"):
            contamination_dict[sample]=contamination
        if (sample_number<(len_sample_list-1)):
            if (contamination_dict[sample_next]=="No"):
                contamination_dict[sample_next]=contamination_next
        if (sample_number!=0):
            if (contamination_dict[sample_previous]=="No"):
                contamination_dict[sample_previous]=contamination_previous

        if (contamination_referral_dict[sample]=="No"):
            contamination_referral_dict[sample]=contamination_referral




        #calculate the level of contamination in the NTC. If the number of fusions for the panel genes is less 0, contamination is "No". Otherwise contamination is "Yes".
        elif (NTC_in_sample==True):
            ntc_fusions=pandas.read_csv("./Gathered_Results/Database/"+sample+"_fusion_check.csv")

            ntc_splice_fusions=ntc_fusions[(ntc_fusions.type=="Splice")]
            if ((len(ntc_splice_fusions))!=0):
                ntc_splice_fusions['fusion'] = ntc_splice_fusions['fusion'].apply(lambda x: x+"-"+x)
                fusions_ntc=ntc_fusions[(ntc_fusions.type=="Fusion")]
                ntc_report=ntc_splice_fusions.append(fusions_ntc)    



            if (len(ntc_report)>0):
                ntc_report=ntc_report[ntc_report["fusion"].str.contains(tumour_referrals)]
                len_ntc_report=len(ntc_report)
                if (len_ntc_report>0):
                    contamination_dict[sample]="Yes"
                    contamination_referral_dict[sample]="Yes"
 		
        #Increase the sample_number by 1   
        sample_number=sample_number+1


#Convert the contamination dictionaries to a dataframe and output to a csv file
contamination_panel_dataframe=pandas.DataFrame(list(contamination_dict.items()), columns=["Sample", "Contamination"])
contamination_referral_dataframe=pandas.DataFrame(list(contamination_referral_dict.items()), columns=["Sample", "Contamination_referral"])
contamination_dataframe=pandas.merge(contamination_panel_dataframe, contamination_referral_dataframe, on="Sample")
contamination_dataframe.to_csv("contamination-"+worksheetid+".csv", index=False)












