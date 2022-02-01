import pandas
import argparse
import decimal
from decimal import Decimal
    
parser=argparse.ArgumentParser()
parser.add_argument('--sampleId', required=True)
parser.add_argument('--referral', required=True)
parser.add_argument('--gaps_path', required=True)
parser.add_argument('--bedfile_path', required=True)
args=parser.parse_args()

sampleId=args.sampleId
referral=args.referral
gaps_path=args.gaps_path
bedfile_path=args.bedfile_path


#load file
file=pandas.read_csv(gaps_path+sampleId+"_"+referral+"_intersect.txt", sep='\t', names=['Chr','Start', 'End', 'Info', 'Chr_cosmic','Start_cosmic', 'End_cosmic', '7', '8','9', '10', 'Counts'])

if (referral=="Melanoma" or referral=="Lung" or referral=="Glioma" or referral=="Colorectal"):

    if (file.shape[0]!=0):

        #remove regions where there is no overlap
        file=file[file['Start_cosmic'] != -1]

    if (file.shape[0]!=0):

        #only keep certain columns
        file[['Gene','Ignore']] = file.Info.str.split("(",expand=True,)
        file=file.filter(items=['Chr', 'Start', 'End','Info', 'Gene', 'Counts'])

        #combine the rows for the same region and add the counts column for these rows
        file['Counts']=file['Counts'].apply(lambda x: int(x))
        grouped_file= file.groupby(['Chr', 'Start', 'End', 'Info', 'Gene'], as_index=False).aggregate({'Counts': 'sum'})

        #Create a dictionary of the number of variants in the cosmic file for each of the genes 
        file=pandas.read_csv(bedfile_path+referral+".bed", sep='\t', names=['Chr','Start', 'End', 'Gene', 'ENST','ENSP','Info', 'Counts'])
        genes=file['Gene']
        genes=list(set(genes))

        gene_counts_dict = {}
        for gene in genes:
            gene_counts=file[file["Gene"]==gene]
            gene_counts['Counts']=gene_counts['Counts'].apply(lambda x: int(x))
            Total = gene_counts['Counts'].sum()
            gene_counts_dict[gene] = Total

        #Add the total counts and percentage column for each of the gaps
        grouped_file['Total_counts']=grouped_file['Gene'].apply(lambda x: gene_counts_dict.get(x))
        grouped_file['Percentage']=(grouped_file['Counts']/grouped_file['Total_counts'])*100

        #Round the Percentage column to 2dp
        decimal.getcontext().rounding=decimal.ROUND_HALF_UP
        grouped_file['Percentage']=grouped_file['Percentage'].apply(lambda x: float((Decimal(str(x)).quantize(Decimal('0.01')))))

        #output table to csv file
        grouped_file=grouped_file.filter(items=['Chr', 'Start', 'End','Info', 'Gene', 'Counts', 'Percentage'])
        grouped_file.to_csv(gaps_path+sampleId+"_"+referral+"_output.csv", sep=',', index=False)

    else:
        file=file.filter(items=['Chr', 'Start', 'End','Info', 'Gene', 'Counts', 'Percentage' ])
        file.to_csv(gaps_path+sampleId+"_"+referral+"_output.csv", sep=',', index=False)

else:

    file=pandas.read_csv(gaps_path+sampleId+"_"+referral+"_hotspots.gaps", sep='\t', names=['Chr','Start', 'End', 'Info'])
    file['Counts']='N/A'
    file['Percentage']='N/A'
    file.to_csv(gaps_path+sampleId+"_"+referral+"_output.csv", sep=',', index=False)


