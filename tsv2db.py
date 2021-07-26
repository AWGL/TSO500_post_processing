import argparse
import decimal
from decimal import Decimal

"""
- script to format small variants part of variant.tsv files from the illumina localapp output to new database friendly format
- For use as part of TSO500 pipeline
- usage: python tsv2db.py -t <tsvfile path> -n <NTC tsv file> -o <output path>
- output: tab delimited table with columns: gene, chr, pos, ref, alt, vaf, depth, hgvs_p, hgvs_c, consequence, exon, alt_reads, in_ntc, ntc_vaf, ntc_depth, ntc_alt_reads
"""


###############################
###        Functions        ###
###############################

def parse_ntc_tsv(ntcfile):
	"""
	- returns a list of list of variant information (chr, start, end, meta, vaf, depth) for ntc or ['NA'] if no variants found
	"""

	## create base variables
	start_parsing = False
	rowcounter = 0
	ntc_var_list = []

	tsv = open(ntcfile, 'r')
	for tsv_line in tsv:

		## look for dataset_title before starting to parse the subsequent lines
		if start_parsing == False:
			if tsv_line.startswith('[Small Variants]'):
				start_parsing = True
				pass

			else:
				pass

		else:

			## treat first line differently as it is a header line
			if rowcounter == 0:
				rowcounter += 1

			else:
				if tsv_line.strip() == 'NA':
					ntc_var_list.append('NA')
					break

				else:
					var_id = tsv_line.split('\t')[0:7]
					ntc_var_list.append(var_id)

	tsv.close()

	return ntc_var_list


def parse_sample_tsv(tsvfile, ntcvarlist):
	"""
	- returns a list of variant informaiton for sample in format:
	'gene', 'chr', 'pos', 'ref', 'alt', 'vaf', 'depth', 'hgvs_p', 'hgvs_c', 'consequence', 'exon', 'alt_reads', 'in_ntc', 'ntc_vaf', 'ntc_depth'. 'ntc_alt_reads'
	- ntc_vaf and ntc_depth are only populated if in_ntc = TRUE
	"""

	## create base variables
	start_parsing = False
	rowcounter = 0
	var_list = []

	sampletsv = open(tsvfile, 'r')
	for tsv_line in sampletsv:

		## look for dataset_title before starting to parse the subsequent lines
		if start_parsing == False:
			if tsv_line.startswith('[Small Variants]'):
				start_parsing = True
				pass

			else:
				pass
		else:

			## treat first line differently as it is a header line
			if rowcounter == 0:
				rowcounter += 1

			else:

				## parse small variants data
				tsv_line1 = tsv_line.replace('\n','')
				tsv_line_list = tsv_line1.split('\t')

				## check if variant is in NTC file
				var_id_list = tsv_line_list[0:5]
				in_ntc_count = 0

				for item in ntcvarlist:
					
					## if var is in ntc var, create counter and then 
					if item[0:5] == var_id_list:
						in_ntc_count += 1
						ntc_vaf = item[5]
						ntc_depth = item[6]
						ntc_allele_percent = float(ntc_vaf) * 100
						ntc_alt_reads = (int(ntc_depth) / 100) * ntc_allele_percent
						ntc_alt_reads_rounded = int((Decimal(str(ntc_alt_reads)).quantize(Decimal('1'))))

				## if in ntc, then give out ntc vaf, depth, and alt read count
				if in_ntc_count > 0:
					in_ntc = f'True\t{ntc_vaf}\t{ntc_depth}\t{ntc_alt_reads_rounded}'

				else:
					in_ntc = 'False'

				## get alt_reads
				## special case to exclude blank lines at end of tsv
				if len(tsv_line_list) > 5:
					allele_percent = float(tsv_line_list[5]) * 100
					alt_reads = (int(tsv_line_list[6]) / 100) * allele_percent
					alt_reads_rounded = int((Decimal(str(alt_reads)).quantize(Decimal('1'))))

					## append info in correct order
					var_info = f'{tsv_line1}\t{alt_reads_rounded}\t{in_ntc}'

					## add line list into main list for exporting later
					var_list.append(var_info)
	sampletsv.close()

	return var_list


###############################
###        Programme        ###
###############################


if __name__ == '__main__':

	## args
	parser = argparse.ArgumentParser()
	parser.add_argument('--tsvfile','-t', help = 'pathway to tsv combined variants file')
	parser.add_argument('--ntcfile','-n', help = 'pathway to NTC tsv combined variants file')
	parser.add_argument('--outfile','-o', help = 'filename/path for output')
	args = parser.parse_args()

	## set Decimal settings
	decimal.getcontext().rounding = decimal.ROUND_HALF_UP

	### parse NTC file
	ntc_var_list = parse_ntc_tsv(args.ntcfile)
	print(ntc_var_list)

	### parse sample file
	var_list = parse_sample_tsv(args.tsvfile, ntc_var_list)


	### export data
	## create header line
	outfile_headers = ['gene', 'chr', 'pos', 'ref', 'alt', 'vaf', 'depth', 'hgvs_p', 'hgvs_c', 'consequence', 'exon', 'alt_reads', 'in_ntc', 'ntc_vaf', 'ntc_depth', 'ntc_alt_reads']
	outfile_headers_string = '\t'.join(outfile_headers)

	## export information to file
	with open(args.outfile,'w',newline='') as file:
		file.write(outfile_headers_string + '\n')

		for line in var_list:
			file.write(line + '\n')