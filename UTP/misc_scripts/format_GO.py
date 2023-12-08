#!/usr/bin/env python3

# GO-term formatter for GO-term enrichment analysis
# Made by Alejandro Baars as part of his General Research Profile
# Utrecht University 2023

# importing packages
import os, sys, glob, argparse as ap, datetime

# define script
gof = ap.ArgumentParser(
	prog = 'GO-term table formatter',
	usage = 'format_GO.py --input [FILE] --annotation [FILE] --outdir [DIR] --outname [STR] '
		'--GO_ID [INT] --evidence_str [STR] [STR1 STR2 STRN...] --evidence_col [INT] '
		'--gene_ID [INT] --cut_GO_ID [STR] --gene_names [FILE]',
	description = 'Script to properly format GO-terms and annotation files to install '
		'as an annotation package. For columns, the first column is considered 1.')

# add arguments to the script
gof.add_argument('--version', '-v',
	action = 'version',
	version = 'a0.1.1')
gof.add_argument('--input','-i',
	action = 'store',
#	required = True,
	help = '<Required> Input for your table. Accepts a single file.')
gof.add_argument('--annotation','-a',
	action = 'store',
#	required = True,
	help = '<Required> Input for your annotation file (GFF/GFF3/GTF/BED).')
gof.add_argument('--outdir', '-o',
	action = 'store',
#	required = True,
	help = '<Required> Specify the output directory.')
gof.add_argument('--outname', '-n',
	action = 'store',
	default = 'GO_table.csv',
	help = 'Specify the name of your output table. Default = GO_table.csv.')
gof.add_argument('--in_delim',
	action = 'store',
	default = ',',
	help = 'Specify the delimiter of your input table. Default = , (assumes csv files).')
gof.add_argument('--GO_ID',
	action = 'store',
#	required = True,
	help = '<Required> Specify the column in which the GO-term ID\'s can be found.')
gof.add_argument('--evidence_col',
	action = 'store',
	default = 0,
	help = 'Specify the column in which the evidence codes can be found.')
gof.add_argument('--evidence_str',
	nargs = '+',
	action = 'store',
	default = ' ',
	help = 'If there is no column with evidence codes, submit one or more codes manually. They will be put on all entries.')
gof.add_argument('--gene_ID',
	action = 'store',
#	required = True,
	help = '<Required> Specify the column in which the gene ID\'s can be found.')
gof.add_argument('--cut_GO_ID',
	action = 'store',
	default = ' ',
	help = 'If the GO-term is accompanied by a description, or if there are multiple GO-terms in the same column, specify the character '
		'with which they can be separated.')
gof.add_argument('--skip_lines',
	action = 'store',
	default = 0,
	help = 'Specify whether lines (e.g. header lines) should be skipped. Default = 0.')
gof.add_argument('--gene_names',
	action = 'store',
	default = ' ',
	help = 'Specify a csv table with gene IDs and corresponding gene names (optional)')

# PROCESS ARGUMENTS
namespace = gof.parse_args()

## IO
input = str(namespace.input)
annot = str(namespace.annotation)
outdir = str(namespace.outdir)
outname = str(namespace.outname)
gene_names = str(namespace.gene_names)

## columns
col_go_id = int(namespace.GO_ID)
col_gene_id = int(namespace.gene_ID)
col_evidence = int(namespace.evidence_col)
if len(namespace.evidence_str) > 1:
	lst_evidence = list(namespace.evidence_str)
	str_evidence = ';'.join(lst_evidence)
else:
	lst_evidence = list(namespace.evidence_str)	# the input is a string by default
	str_evidence = lst_evidence[0]

## misc
cut_go_id = str(namespace.cut_GO_ID)
in_delim = str(namespace.in_delim)
skip_lines = int(namespace.skip_lines)

## check if submitted files and directory exist
if not os.path.isfile(input):
	sys.exit('ERROR: Cannot find submitted file ' + input + '. Exiting.')
if not os.path.isfile(annot):
	sys.exit('ERROR: Cannot find submitted file ' + annot + '. Exiting.')
if gene_names != ' ':
	if not os.path.isfile(gene_names):
		sys.exit('ERROR: Cannot find submitted file ' + gene_names + '. Exiting.')
if not os.path.isdir(outdir):
	print('Output directory ' + outdir + ' not found, creating it for you')
	os.mkdir(outdir)

# check if at least 1 thing is submitted for evidence
if col_evidence == 0:
	if str_evidence == ' ':
		sys.exit('ERROR: The scripts needs either evidence_col or evidence_str')



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

## function to read and re-organize the annotation data
def get_annotation(input_annotation, gene_names):
	## read and check input data
	with open(input_annotation, 'r') as file:
		input_lines = file.readlines()
	annot_extension = input_annotation.split('.')[-1]
	if annot_extension not in ['gff', 'gtf', 'gff3', 'bed']:
		sys.exit('ERROR: Input annotation file ' + input_annotation + ' has an unknown extension. Exiting.')
	if gene_names != ' ':
		with open(gene_names, 'r') as gfile:
			gene_names_list = gfile.readlines()
		test_line = gene_names_list[5]
		if len(test_line.split(',')) != 2:
			sys.exit('ERROR: Incorrect number of columns in ' + gene_names + '. Two columns, separated by a comma, '
				' are required. Exiting.')

	## create empty output lists
	gene_symbol_list = []
	gene_chroms_list = []

	## extract necessary data per annotation format
	counter = 0
	if annot_extension in ['gff', 'gtf', 'gff3']:
		header = False
		for line in input_lines:
			if line.startswith('#'):			# check headers for info if a header is present
				header = True
				if annot_extension == 'gff':
					if 'gff-version 3' in line:		# sometimes gff3 uses the extension gff, so this is to check
						annot_extension = 'gff3'	# and make sure the right tools are used
					elif 'gff-version 2' in line:	# gff v2 = gtf
						annot_extension = 'gtf'
					else:
						print('Warning: GFF version cannot be inferred from header. Assuming GFF version 3.')
						annot_extension = 'gff3'
			else:
				if header == False:
					print('Warning: GFF version cannot be inferred from header. Assuming GFF version 3.')
					annot_extension = 'gff3'
				line_split = line.split('\t')
				prefix_suffix = ['\n','ID=','Id=','id=','\"',
					'gene-','-gene','gene|','|gene','gene_id',
					'transcript-','-transcript','transcript|','|transcript','transcript_id']
#				if counter <= 12:
#					print(line_split) #tmpprint
				if line_split[2] in ['gene', 'transcript']:
#					if counter <= 12:
#						print(line_split[2]) #tmpprint
					chrom = line_split[0]
					attributes = line_split[8]
					attributes_split = attributes.split(';')
					# take the correct data based on the gff format
					if annot_extension == 'gff3':
						gene_ID = attributes_split[0]
						if not gene_ID.startswith('ID='):
							sys.exit('ERROR: Attribute found ID ' + gene_ID + ' does not begin with "ID=".\n'
								'Check your GFF file\'s format!')
						## start formatting the string
						for ps in prefix_suffix:
							gene_ID = remove_prefix(gene_ID, ps)
							gene_ID = remove_suffix(gene_ID, ps)
					elif annot_extension == 'gtf':
						gene_ID = attributes_split[0]
						for ps in prefix_suffix:
							gene_ID = remove_prefix(gene_ID, ps)
							gene_ID = remove_suffix(gene_ID, ps)
					# pair the gene name to the geneID, if names were submitted
					if gene_names != ' ':
						for gn in gene_names_list:
							gns = gn.split(',')
							if gene_ID in gns:
								gene_name = gns
					else:
						gene_name = gene_ID

					# store info in output lists
					gene_symbol_line = '\t'.join([gene_ID, gene_ID, gene_name])
					if gene_symbol_line.endswith('\n'):
						remove_suffix(gene_symbol_line, '\n')
					gene_symbol_line = gene_symbol_line + '\n'
					gene_symbol_list.append(gene_symbol_line)
					gene_chroms_line = '\t'.join([gene_ID, chrom])
					if gene_chroms_line.endswith('\n'):
						remove_suffix(gene_chroms_line, '\n')
					gene_chroms_line = gene_chroms_line + '\n'
					gene_chroms_list.append(gene_chroms_line)
				counter += 1 #tmpprint

	else:
		for line in input_lines:
			if not line.startswith('#'):
				line_split = line.split('\t')
				if len(line_split) < 4:
					sys.exit('ERROR: BED files need the "name" field in the 4th column. Exiting.')
				chrom = line_split[0]
				gene_ID = line_split[3]
				# pair the gene name to the geneID, if names were submitted
				if gene_names != ' ':
					for gn in gene_names_list:
						gns = gn.split(',')
						if gene_ID in gns:
							gene_name = gns
				else:
					gene_name = gene_ID

				# store info in output lists
				gene_symbol_line = '\t'.join([gene_ID, gene_ID, gene_name])
				gene_symbol_line = gene_symbol_line + '\n'
				gene_symbol_list.append(gene_symbol_line)
				gene_chroms_line = '\t'.join([gene_ID, chrom])
				gene_chroms_line = gene_chroms_line + '\n'
				gene_chroms_list.append(gene_chroms_line)
	return gene_symbol_list, gene_chroms_list

## function to get the GO-term table
def get_go(input_file, col_evidence, str_evidence, col_go_id, col_gene_id):
	## read input data
	with open(input_file, 'r') as file:
		input_lines = file.readlines()

	## set up empty list and counter
	output_lines = []
	cnt = 0	#counter

	## go through the input line by line
	for line in input_lines:
		if cnt >= skip_lines:	# skips lines before a specified number
			line_cut = line.split(in_delim)
			## assign the right data to the right variable
			if cut_go_id != ' ':	# in this case, the GO ID's contain multiple sub-columns which must be separated
				go_id_tmp = line_cut[(col_go_id - 1)]
				go_id_cut = go_id_tmp.split(cut_go_id)  # cutting with the submitted delimiter
				for elem in go_id_cut:	# looping through the elements in the GO ID field
					if elem.startswith('GO:'):	# and selecting only those that start with 'GO:' to remove any comments/summaries
						go_id = elem	# assign the first value
						if col_evidence != 0:	# take value from column if int
							evidence = line_cut[(col_evidence - 1)]
						else:					# take value from input as str
							evidence = str_evidence
						gene_id = line_cut[(col_gene_id - 1)]
					newline = ','.join([gene_id, go_id, evidence])
			else:	# in this case, the GO ID is a single value
				go_id = line_cut[(col_go_id - 1)]
				if col_evidence != 0:	# take value from column if int
					evidence = line_cut[(col_evidence - 1)]
				else:					# take value from input as str
					evidence = str_evidence
				gene_id = line_cut[(col_gene_id - 1)]
				newline = ','.join([gene_id, go_id, evidence])
			if not newline.endswith('\n'):
				newline = newline + '\n'
			output_lines.append(newline)
		cnt += 1
	return output_lines



# PROCESSING DATA
## get name and chromosome tables
tables_tuple = get_annotation(annot, gene_names)
table_sym = tables_tuple[0]
table_chrom = tables_tuple[1]

## get GO table
table_GO = get_go(input, col_evidence, str_evidence, col_go_id, col_gene_id)

## tmpprint
#print(table_GO[0:5])
#print(table_sym[0:5])
#print(table_chrom[0:5])



# WRITE OUTPUT
## general check
if not outdir.endswith('/'):
	outdir = outdir + '/'

## GO-terms
outfile_GO = outdir + outname + '.GO.reorganized.csv'

with open(outfile_GO, 'w') as out_GO:
	for outline in table_GO:
		out_GO.write(outline)

## Symbols
outfile_sym = outdir + outname + '.symbols.csv'

with open(outfile_sym, 'w') as out_sym:
	for outline in table_sym:
		out_sym.write(outline)

## Chromosomes
outfile_chrom = outdir + outname + '.chroms.csv'

with open(outfile_chrom, 'w') as out_chrom:
	for outline in table_chrom:
		out_chrom.write(outline)
