#!/bin/python3

# Script to change the gene lengths column in a featureCounts counts file
# Made by Alejandro Baars as part of his General Research Profile
# Utrecht University 2023

# IMPORT PACKAGES
import re, argparse as ap
from itertools import islice


# SET ARGUMENTS
cgl = ap.ArgumentParser(
	prog = 'Change gene length',
	usage = 'change_gene_length.py --input [path/to/input.counts] --output [path/to/output.counts] --annotation [path/to/annotation.gff3]',
	description = 'Change the values in the "length" column in a counts file from the gene lengths to the mRNA transcript lengths.')

cgl.add_argument('--version', '-v',
	action = 'version',
	version = '0.2.0',
	help = 'Prints the script version and then exits.')
cgl.add_argument('--input', '-i',
	action = 'store',
	required = True,
	help = '<Required> Provide the path and filename of your counts file.')
cgl.add_argument('--output', '-o',
	action = 'store',
	required = True,
	help = '<Required> Provide the path and filename for your output file.')
cgl.add_argument('--annotation', '-a',
	action = 'store',
	required = True,
	help = '<Required> Provide the path and filename of your annotation file (GFF v3 format only)')
cgl.add_argument('--header',
	action = 'store',
	default = 2,
	help = 'Specify how many lines are header lines (default = 2). They are parsed straight into the output.')



# PROCESS ARGUMENTS
args = cgl.parse_args()

input = str(args.input)
output = str(args.output)
annot = str(args.annotation)
header = int(args.header)



# READ ANNOTATION
## read file as list of lines & keep lines of interest
annotation = {}									# dictionary with gene/transcript names and lengths
print('Reading annotation file ' + annot)

with open(annot, 'r') as afile:
	for line in afile:
		if not line.startswith('#'):			# filter out header/comment lines
			line_split = line.split('\t')
			if line_split[2] == 'exon':			# keep only exon lines

				# calculate exon length
				start = int(line_split[3])
				stop = int(line_split[4])
				length = stop - start + 1		# both positions are 1-based

				# get exon ID
				attributes = line_split[8].split(';')
				for att in attributes:
					if att.startswith('Parent='):
						att_name = att.split('=')
						exonID = att_name[1]

				# add to annotation dictionary
				if exonID in annotation:				# if the entry exists already (from previous exons) we
					tmp_length = annotation[exonID]		# sum the values and add the sum to the dictionary
					new_length = tmp_length + length
					annotation[exonID] = new_length
				else:									# if the entry does not exist, simply create it
					annotation[exonID] = length



# READ AND PROCESS COUNTS
print('Reading counts file ' + input)
outlist = []							# append output data here
with open(input, 'r') as ifile:
	data = ifile.readlines()

file_length = len(data) - header
cnt = 1

print('No. counts found: ' + str(file_length))

for line in data:

	if cnt <= header:				# appends the first 2 lines (as they are header lines)
		outlist.append(line)
	else:
		line_split = line.split('\t')
		transcript_id = line_split[0]

		# filter the right lines from the annotation
		if transcript_id in annotation:
			length = annotation[transcript_id]
		else:
			length = line_split[5]				# keep the old value if a new one cannot be found

		# exchange the old length for the new one
		line_split[5] = str(length)
		newline = '\t'.join(line_split)
		outlist.append(newline)

	cnt += 1



# WRITE DATA TO OUTPUT
print('Writing output to ' + output)
with open(output, 'w') as ofile:
	for line in outlist:
		ofile.write(line)

print('Gene lengths successfully changed')
