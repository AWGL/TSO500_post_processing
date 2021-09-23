import sys
import json

# filter JSON output from nirvana to give output ready for variant database upload (minus the NTC checks)

# open JSON
with open(sys.argv[1], 'r') as f:
    y = json.load(f)

# make empty list for results output
out_list = []

# loop through each variant and pull out info if the filter is pass (based on structure of JSON)
for var in y['positions']:
    if 'PASS' in var['filters']:

        # info common to variant (regarless of transcript)
        chr = var['chromosome']
        pos = str(var['position'])
        ref = var['refAllele']

        vaf = var['samples'][0]['variantFrequencies'][0]
        depth = var['samples'][0]['totalDepth']
        alt_reads = int(vaf * depth)

        # could be multiallelic so loop through each option
        for variant in var['variants']:

            alt = variant['altAllele']

            # transcript specific info - use canonical refseq
            for transcript in variant['transcripts']:
                if transcript['source'] == 'RefSeq' and 'isCanonical' in transcript.keys():
                    if transcript['isCanonical']:

                        # keys are only present in JSON if they're not empty
                        if 'hgnc' in transcript.keys():
                            gene = transcript['hgnc']
                        else:
                            gene = ''

                        if 'hgvsp' in transcript.keys():
                            hgvs_p = transcript['hgvsp']
                        else:
                            hgvs_p = ''

                        if 'hgvsc' in transcript.keys():
                            hgvs_c = transcript['hgvsc']
                        else:
                            hgvs_c = ''

                        if 'consequence' in transcript.keys():
                            csq = ':'.join(transcript['consequence'])
                        else:
                            csq = ''

                        if 'exon' in transcript.keys():
                            exon = transcript['exon']
                        else:
                            exon = ''

        # add to variant list
        out_list.append(
            [gene, chr, pos, ref, alt, str(vaf), str(depth), hgvs_p, hgvs_c, csq, exon]
        )


# print to screen - redirect to output within pipeline
for var in out_list:
    print('\t'.join(var))
