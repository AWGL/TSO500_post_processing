import sys
from pysam import VariantFile
import decimal

# set decimal to always round down
decimal.getcontext().rounding = decimal.ROUND_DOWN

# load in VCF and make pysam object. Must be gzipped and tabix indexed
input_vcf = sys.argv[1]
vcf = VariantFile(input_vcf)

# list to store output
output_list = []

# list of the variants that we want to ignore the blacklist on, plus associated annotations
# Add to this list if any more variants need removing from the blacklist
variants_to_keep = [
    {
        'chr': 'chr20',
        'pos': 31022441,
        'gene': 'ASXL1',
        'ref': 'A',
        'alt': 'AG',
        'hgvs_c': 'NM_015338.5:c.1934dup',
        'hgvs_p': 'NP_056153.2:p.(Gly646TrpfsTer12)',
        'csq': 'frameshift_variant',
        'exon': '13/17',
        'cutoff': 0.1,
    },
    {
        'chr': 'chr10',
        'pos': 89720653,
        'gene': 'PTEN',
        'ref': 'C',
        'alt': 'CA',
        'hgvs_c': 'NM_000314.6:c.808dup',
        'hgvs_p': 'NP_000305.3:p.(Met270AsnfsTer28)',
        'csq': 'frameshift_variant',
        'exon': '8/9',
        'cutoff': 0,
    },
    {
        'chr': 'chr10',
        'pos': 89720679,
        'gene': 'PTEN',
        'ref': 'C',
        'alt': 'T',
        'hgvs_c': 'NM_000314.6:c.830C>T',
        'hgvs_p': 'NP_000305.3:p.(Thr277Ile)',
        'csq': 'missense_variant',
        'exon': '8/9',
        'cutoff': 0,
    },
]

# search VCF for each of the variants we'd like to keep
for var in variants_to_keep:

    # loop through all variants at the position (might be multiple as each alt has a new line)
    for rec in vcf.fetch(var['chr'], var['pos'] -1, var['pos']):

        # get all filters for variant
        filters = rec.filter.keys()

        # remove blacklist from filters list
        if 'Blacklist' in filters:
            filters.remove('Blacklist')

        # if no more filters (i.e. blacklist was the only filter), then add variant to report
        if len(filters) == 0:

            # only pull through specific base changes
            if var['ref'] == rec.ref and var['alt'] == rec.alts[0]:
                base_changes_pass = True
            else:
                base_changes_pass = False

            # pull out VAF, make decimal and round to 4dp
            vaf = decimal.Decimal(rec.samples[0]['VF'][0]).quantize(decimal.Decimal('0.0001'))

            # some should only be pulled through above a certain VAF, defined in config above. If no cutoff then config will be 0
            if vaf > var['cutoff']:
                vaf_pass = True
            else:
                vaf_pass = False

            # if both filters pass then add to output
            if vaf_pass and base_changes_pass:

                output_list.append([
                    var['gene'],
                    rec.chrom, 
                    str(rec.pos), 
                    rec.ref, 
                    rec.alts[0],
                    str(vaf),
                    str(rec.samples[0]['DP']),
                    var['hgvs_p'],
                    var['hgvs_c'],
                    var['csq'],
                    var['exon'],
                ])

# print to screen - redirect to output within pipeline        
for var in output_list:
    print('\t'.join(var))
