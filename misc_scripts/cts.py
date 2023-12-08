#!/bin/python3

# Script to calculate transcriptome size based on an annotation file
# Made by Alejandro Baars as part of his General Research Profile
# Utrecht University 2023

# IMPORT PACKAGES
import os, sys, argparse as ap

# SET ARGUMENTS
cts = ap.ArgumentParser(
	prog = 'Calculate transcriptome size',
	usage = 'cts.py --input [BED/GFF] --element [gene/transcript/exon]',
	description = 'Calculate the size of a transcriptome given an annotation file. Returns the number '
		'of bases covered by the specified element')

cts.add_argument('--version', '-v',
	action = 'version',
	version = '0.1.0')
cts.add_argument('--input', '-i',
	action = 'store',
	required = True,
	help = '<Required> Specify an annotation file to perform the calculations with.')
cts.add_argument('--element', '-e',
	action = 'store',
	default = 'gene',
	help = 'Specify the element in the annotation file to be used for transcriptome calculation (default: \'gene\'). Works only with GFF/GTF/GFF3 files.')



# PROCESS ARGUMENTS

args = cts.parse_args()

input = str(args.input)
element = str(args.element)



# READ FILE

with open(input, 'r') as file:
	data = file.readlines()



# CALCULATE TRANSCRIPTOME

t_list = []

if input.endswith(('gff','gtf','gff3')):
	for line in data:
		if not line.startswith('#'):
			line_split = line.split('\t')
			if line_split[2] == element:
				start = int(line_split[3])
				end = int(line_split[4])
				length = end - start + 1 #both start and end positions are 1-based
				t_list.append(length)
elif input.endswith('bed'):
	for line in data:
		if not line.startswith('#'):
			line_split = line.split('\t')
			start = int(line_split[1])
			end = int(line_split[2])
			length = end - start #bed-format specifies that start is 0-based and the end is 1-based
			t_list.append(length)
else:
	sys.exit('ERROR: File suffix not recognized. This script supports BED, GFF, GFF3 and GTF files only.')

t_size = sum(t_list)



# RETURN OUTPUT

print('Total transcriptome size:\t' + str(t_size) + 'bp')
print('Number of ' + element + 's present:\t' + str(len(t_list)))
print('Largest ' + element + ':\t' + str(max(t_list)) + 'bp')
print('Smallest ' + element + ':\t' + str(min(t_list)) + 'bp')
