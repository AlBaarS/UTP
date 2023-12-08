#!/bin/python3

# Script to rename geneIDs in difExp.R output.csv files
# Made by Alejandro Baars as part of his General Research Profile
# Utrecht University 2023

# IMPORT PACKAGES
import os, sys, argparse as ap



# SET ARGUMENTS
rgi = ap.ArgumentParser(
	prog = 'Rename geneIDs',
	usage = 'rename_geneIDs.py --input [path/to/table.csv] --annotation [path/to/annotation.gff3] --output [path/to/output.csv]'
		'--incol [integer] --feature [string] --attribute [string] --separator [separator] --header [integer] --merge',
	description = 'Rename geneIDs in a table (e.g. output from difExp.R) to an attribute in an annotation file (gff3 format)')

rgi.add_argument('--input', '-i',
	action = 'store',
	required = True,
	help = '<Required> Specify the path and name of your input table.')
rgi.add_argument('--annotation', '-a',
	action = 'store',
	required = True,
	help = '<Required> Specify the path and name of your annotation file.')
rgi.add_argument('--output', '-o',
	action = 'store',
	required = True,
	help = '<Required> Specify the path and name of your output file.')
rgi.add_argument('--incol',
	action = 'store',
	default = '1',
	help = 'Specify the column in the input file containing the geneIDs (default = 1).')
rgi.add_argument('--feature',
	action = 'store',
	default = 'mRNA',
	help = 'Specify the feature type in which the current IDs can be found (default = "mRNA").')
rgi.add_argument('--attribute',
	action = 'store',
	default = 'gene',
	help = 'Specify the attribute to be used for name changing (default = "gene"). Needs to match <attribute>=<name> in the gff file.')
rgi.add_argument('--separator',
	action = 'store',
	default = ',',
	help = 'Specify the separator to be used in the input and output file (default = ","). Specify "none" if the file contains a single column.')
rgi.add_argument('--header',
	action = 'store',
	default = 0,
	help = 'Specify how many lines should be skipped (default = 0).')
rgi.add_argument('--merge',
	action = 'store_true',
	default = False,
	help = 'Specify if the changed name should be replacing the old one (e.g. G0001.rna-XM001 to Gene0001) or be merged (to Gene0001|G0001.rna-XM001) (default = False).')



# PROCESS ARGUMENTS
args = rgi.parse_args()

input = str(args.input)
annot = str(args.annotation)
output = str(args.output)
incol = int(args.incol)
feature = str(args.feature)
attribute = str(args.attribute)
separator = str(args.separator)
header = bool(args.header)
merge = bool(args.merge)

if separator == 'tab':
  separator = '\t'



# FUNCTIONS
## functions to remove prefixes and suffixes from strings
def remove_prefix(text, prefix):
	if text.startswith(prefix):
		return text[len(prefix):]
	else:
		return text

def remove_suffix(text, suffix):
	if text.endswith(suffix):
		return text[:-len(suffix)]
	else:
		return text



# READ ANNOTATION
## read file as list of lines & keep lines of interest
annotation = []
with open(annot, 'r') as afile:
	for line in afile:
		if not line.startswith('#'):			# filter out header/comment lines
			line_split = line.split('\t')
			if line_split[2] == feature:		# keep lines of interest
				annotation.append(line_split)	# add line to list



# READ & PROCESS INPUT
outlist = []
with open(input, 'r') as ifile:
	t = 0
	for line in ifile:
		if t >= header:
			# define the element name that we wish to change
			if separator != 'none':
				line_split = line.split(separator)
				newline_split = line_split
				name = line_split[(incol - 1)]
			else:
				name = remove_suffix(line, '\n')
			# search for match in annotation
			for a in annotation:
				if name in a[-1]:
					match = a[-1]
					attributes = match.split(';')
					# change the name
					for atr in attributes:
						if str(attribute + '=') in atr:
							geneID = remove_prefix(atr, str(attribute + '='))
							if merge:
								newname = str(geneID + '|' + name)
								if separator != 'none':
									newline_split[(incol - 1)] = newname
								else:
									newline = newname
							else:
								if separator != 'none':
									newline_split[(incol - 1)] = geneID
								else:
									newline = geneID

			# sometimes a gene name was already included due to it being a pseudogene, or a different RNA name.
			# these cases are fixed/filtered respectively here
			if separator != 'none':
				if newline_split[(incol - 1)].startswith('gene-'):
					name = remove_prefix(newline_split[(incol - 1)], 'gene-')
					if merge:
						newname = str(name + '|' + name)
						newline_split[(incol - 1)] = newname
					else:
						newline_split[(incol - 1)] = name
				elif newline_split[(incol - 1)].startswith('rna-'):
					continue
				newline = separator.join(newline_split)
				outlist.append(newline)
			else:
				if newline.startswith('gene-'):
					name = remove_prefix(newline, 'gene-')
					newline = name
				elif newline.startswith('rna-'):
					continue
				outlist.append(str(newline + '\n'))
		else:
			outlist.append(line)
		t += 1


# WRITE TO OUTPUT FILE
with open(output, 'w') as ofile:
	for line in outlist:
		ofile.write(line)
