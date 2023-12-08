#!/bin/python3

# Script to remove specific columns from a featureCounts counts file
# Made by Alejandro Baars as part of his General Research Profile
# Utrecht University 2023

# IMPORT PACKAGES
import os, sys, argparse as ap, pandas as pd



# SET ARGUMENTS
rmc = ap.ArgumentParser(
	prog = 'Remove counts',
	usage = 'remove_counts.py --input [/path/and/countsfile] --output [/path/and/filename] --samples [samples (numerical)]',
	description = 'Remove columns from a file with featureCounts .counts output. Keeps the original '
		'intact and creates a copy with the desired columns removed at the specified location.')

rmc.add_argument('--version', '-v',
	action = 'version',
	version = '0.1.0')
rmc.add_argument('--input', '-i',
	action = 'store',
	required = True,
	help = '<Required> Specify the path and filename of your input file.')
rmc.add_argument('--output', '-o',
	action = 'store',
	required = True,
	help = '<Required> Specify the path and filename for your output file.')
rmc.add_argument('--samples', '-s',
	nargs = '+',
	action = 'store',
	required = True,
	help = '<Required> Specify the sample(s) which need to be removed. Specify by number, separated by a whitespace (eg: -s 1 2 5)')



# PROCESS ARGUMENTS
args = rmc.parse_args()

input = str(args.input)
output = str(args.output)
samples = list(args.samples)

# convert samples to integers
samples = [eval(i) for i in samples]

# check if input file and output dir exist
if not os.path.isfile(input):
	sys.exit('ERROR: Cannot find file ' + input)

outdir = os.path.dirname(output)
if not os.path.isdir(outdir):
	os.mkdir(outdir)
	print('Could not find output directory ' + outdir + '. Outdir created.')



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



# PROCESSING
# read data
indata = pd.read_table(input, sep = '\t', comment = '#')

# remove specified columns
columns = [x + 5 for x in samples]
outdata = indata.drop(indata.columns[columns], axis = 1)

# OUTPUT
# write to output file
outdata.to_csv(output, sep = '\t', index = False)

# EPILOGUE
print('Counts successfully processed.')
