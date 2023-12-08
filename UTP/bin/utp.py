#!/usr/bin/env python3

# The uniform transciptomics pipeline
# Made by Alejandro Baars as part of his General Research Profile
# Utrecht University 2023

# importing packages
import os, sys, glob, argparse as ap, datetime

#print(sys.version) #tmpprint

utp = ap.ArgumentParser(
	prog = "Uniform transcriptomics pipeline",
	usage = "utp.py --reads [READS] --genome [GENOME(S)] --annotation [ANNOTATION] "
		"--outdir_bam [OUTDIR_BAM] --outdir_counts [OUTDIR_COUNTS] --tmpdir [TMPDIR] [arguments]",
	description = "Pipeline to run gene expression analysis")

# arguments
utp.add_argument('--version','-v',
	action = 'version',
	version = 'UTP v0.1.5')
utp.add_argument('--verbose','-V',
	action = 'store_true',
	default = False,
	help = 'If specified, prints messages as to what it is currently doing. Off by default.')
utp.add_argument('--reads','-r',
	nargs = '+',
	action = 'store',
	required = True,
	help = '<Required> Input for your reads. Accepts a file, list of files, directory, or list of directories.')
utp.add_argument('--genome','-g',
	nargs = '+',
	action = 'store',
	required = True,
	help = '<Required> Input for your genome(s). Accepts a file, list of files, directory, or list of directories. '
		'Cannot use compressed genome files.')
utp.add_argument('--annotation','-a',
	nargs = '+',
	action = 'store',
	required = True,
	help = '<Required> Input for your annotation files. Accepts a file, list of files, directory, or list of directories. '
		'Cannot use compressed annotation files.')
utp.add_argument('--threads','-t',
	action = 'store',
	default = 2,
	help = 'Specify the number of threads used by programmes that can use multithreading. Default = 2.')
utp.add_argument('--single_end',
	action = 'store_true',
	default = False,
	help = 'The pipeline assumes the reads are paired-end. If you use single-end reads, include this tag.')
utp.add_argument('--outdir_bam','-ob',
	action = 'store',
	required = True,
	help = '<Required> Specifies your output directory for the alignment (BAM) files.')
utp.add_argument('--outdir_counts', '-oc',
	action = 'store',
	required = True,
	help = '<Required> Specifies your output directory for the gene counts.')
utp.add_argument('--tmpdir',
	action = 'store',
	required = True,
	help = '<Required> Specify a temporary directory to put intermediate files. Files will be cleared when they are no longer needed.')
utp.add_argument('--keep_tmp_data',
	action = 'store_true',
	default = False,
	help = 'Keep the data written to the specified temporary directory. By default, this directory is cleared.')
utp.add_argument('--copy_tmp_data',
	action = 'store',
	default = 'no',
	help = 'Copy the data written in the specified temporary directory to a directory specified here instead of deleting it.')
utp.add_argument('--old_picard',
	action = 'store_true',
	default = False,
	help = 'Use the old Picard syntax, in case you use an old installation.')
utp.add_argument('--config','-c',
	action = 'store',
	default = '~/UTP/config/utp.config',
	help = 'Specify a (customized) config file to use rather than the default one.')

# process arguments
namespace = utp.parse_args()

## input data
input_reads = list(namespace.reads)
input_genome = list(namespace.genome)
input_annotation = list(namespace.annotation)

## directories & configuration
outdir = str(namespace.outdir_bam)
outcnt = str(namespace.outdir_counts)
tmpdir = str(namespace.tmpdir)
config = str(namespace.config)
threads = int(namespace.threads)
tmpcopy = str(namespace.copy_tmp_data)

## internal optional arguments
verbose = bool(namespace.verbose)
keep_tmp_data = bool(namespace.keep_tmp_data)
single_end = bool(namespace.single_end)
old_picard = bool(namespace.old_picard)



# set functions

## this function disentangles the various methods to provide input for the reads and genomes
def input_untangler(input, mode):
	output = []
	# determine which files to filter:
	if mode == 'reads':
		extensions = ['.fastq','.fastq.gz','.fastq.bz2','.fastq.zip',
			'.fq','.fq.gz','.fq.bz2','.fq.zip',
			'.fasta','.fasta.gz','.fasta.bz2','.fasta.zip']
	elif mode == 'genome':
		extensions = ['.fasta','.fa','.fna']
	elif mode == 'annotation':
		extensions = ['.gtf','.gff','.gff3','.bed']
	# check if the input consists of 1 entry or a list of entries
	if len(input) == 1:
		# check if the input is a directory or a file
		if os.path.isdir(input[0]):
			dir = input[0]
			for file in os.listdir(dir):
				if file.endswith(tuple(extensions)):
					if 'unpaired' not in file:
						file_out = os.path.join(dir, file)
						output.append(file_out)
					else:
						if verbose:
							print('Readfile omitted: ' + file)
		elif os.path.isfile(input[0]):
			file = input[0]
			if str(file).endswith(tuple(extensions)):
				if 'unpaired' not in file:
					output.append(file)
		else:
			print('UTP Error: Files not found; aborting.')
			sys.exit()
	elif len(input) > 1:
		for i in input:
			if os.path.isdir(i):
				dir = i
				for file in os.listdir(dir):
					if file.endswith(tuple(extensions)):
						if 'unpaired' not in file:
							file_out = os.path.join(dir, file)
							output.append(file_out)
			elif os.path.isfile(i):
				file = i
				if file.endswith(tuple(extensions)):
					if 'unpaired' not in file:
						output.append(file)
	# returning a properly formatted list of all input files with path
	return output

## this function merges the genome files and writes it to a new file
def genome_merger(genomes, annotation, dir):
	# setting empty variables to store the information in
	out_file = []
	out_IDs = []
	# going through the data one genome file at a time
	for file in genomes:
		if verbose:
			print('Reading ' + file)
		fname = os.path.basename(file)
		fp = open(file, 'r')
		sname = os.path.basename(file).split('.')[0]
		out_IDs.append(sname)
		out_genome_name = dir + sname + '.genome.fasta'
		out_genome = open(out_genome_name, 'w')
		chrom = []
		for line in fp:
			if '>' in line:
				nline = line.strip('\n')
				nline = nline.split(' ')[0]
				nline = nline + '@' + fname
				out_file.append(nline)
				out_genome.write(str(nline + '\n'))
				name_line = nline.strip('>')
				chrom.append(name_line)
			else:
				out_file.append(line)
				out_genome.write(str(line))
#		out_IDs[file] = chrom
		fp.close()
		# finding the correct annotation(s) for the genome
		annot_genome = str()
		for a in annotation:
			if sname in os.path.basename(a):
				annot_genome = str(a)
				if verbose:
					print('Matched genome ' + fname + ' to annotation ' + os.path.basename(a) + '.')
			# if no match was found: exit
			if not any(sname in an for an in annotation):
				print('UTP Error: Could not find match for ' + fname + '.')
				sys.exit()
		# go through the matched genome file(s)
		ap = open(annot_genome, 'r')
		extension = annot_genome.split('.')[-1]
		out_annot_name = dir + sname + '.annotation.' + extension
		annot_out = open(out_annot_name, 'w')
		for line in ap:
			if '#' in line:
				annot_out.write(line)
			else:
				for c in chrom:
					annot_query = c.split('@')[0]
					if annot_query in line:
						anline = line.replace(annot_query, c, 1)
						annot_out.write(anline)
		annot_out.close()
	# write merged fasta lines to file
	out_genomes_merged = open(dir + 'genomes_merged.fasta', 'w')
	for line in out_file:
		if '>' in line:
			out_genomes_merged.write(line + '\n')
		else:
			out_genomes_merged.write(line)
	out_genomes_merged.close()
	return out_IDs

## this function reads the config file and stores the options for each programme in a dictionary
def config_reader(file):
	if os.path.isfile(file):
		cfile = {}
		with open(file, 'r') as f:
			for line in f:
				if '#' not in line:
					if line != '\n':
						line = line.strip('\n')
						line_split = line.split('~')
						if len(line_split) == 1:
							line_split.append('')
						cfile[str(line_split[0])] = str(line_split[1])
		return cfile
	else:
		print('FqC Error: Cannot find config file: ' + file)
		sys.exit()

## this function detects the current time and date and outputs it in a string format
def timestamp(mode):
	time_stamp = datetime.datetime.now()
	# taking and re-formatting the necessary data from the timestamp
	t_yrs = str(time_stamp.year)
	t_mon = str(time_stamp.month)
	if len(t_mon) == 1:
		t_mon = '0' + t_mon
	t_day = str(time_stamp.day)
	if len(t_day) == 1:
		t_day = '0' + t_day
	t_hrs = str(time_stamp.hour)
	if len(t_hrs) == 1:
		t_hrs = '0' + t_hrs
	t_min = str(time_stamp.minute)
	if len(t_min) == 1:
		t_min = '0' + t_min
	t_sec = str(time_stamp.second)
	if len(t_sec) == 1:
		t_sec = '0' + t_sec
	# output data based on mode
	if mode == 'full':
		timestamp = str(t_yrs + '-' + t_mon + '-' + t_day + '_' + t_hrs + ':' + t_min + ':' + t_sec)
	elif mode == 'time':
		timestamp = str(t_hrs + ':' + t_min + ':' + t_sec)
	elif mode == 'date':
		timestamp = str(t_yrs + '-' + t_mon + '-' + t_day)
	return timestamp


## this function removes a given suffix, if present, from a string
def remove_suffix(text, suffix):
	if text.endswith(suffix):
		return text[:-len(suffix)]
	else:
		return text
## this script was written for an older version of Python (3.6.3), hence the native .removesuffix() method isn't an option

# input processing

## print start
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' UTP run initiated')

## neatly organize the input reads, genomes, and annotations
config = os.path.realpath(config)
input_reads_files = input_untangler(input_reads, 'reads')
input_genome_files = input_untangler(input_genome, 'genome')
input_annotation_files = input_untangler(input_annotation, 'annotation')

## set merge
if len(input_genome_files) > 1:
	merge = True
else:
	merge = False


## collect read filenames
input_reads_names = []
for file in input_reads_files:
	rfname = os.path.basename(file)
	rfname = rfname.split('.')[0]
	input_reads_names.append(rfname)

## check if the input reads are compressed
if input_reads_files[0].endswith('.gz'):
	reads_compressed = True
	reads_decomp_method = 'gunzip -c'
elif input_reads_files[0].endswith('.bz2'):
	reads_compressed = True
	reads_decomp_method = 'bunzip2 -c'
elif input_reads_files[0].endswith('.zip'):
	reads_compressed = True
	reads_decomp_method = 'unzip -c'
else:
	reads_compressed = False

## check if the submitted outdir and tmpdir are existing directories
## bam outdir
if os.path.isdir(outdir):
	if not outdir.endswith('/'):
		outdir = outdir + '/'
else:
	print('UTP Error: UTP cannot find output directory ' + outdir)
	sys.exit()

if os.listdir(outdir) != []:
	print('UTP Warning: Output directory not empty. Existing files may be overwritten.')

## counts outdir
if os.path.isdir(outcnt):
	if not outcnt.endswith('/'):
		outcnt = outcnt + '/'
else:
	print('UTP Error: UTP cannot find output directory ' + outcnt)
	sys.exit()

if os.listdir(outcnt) != []:
	print('UTP Warning: Output directory not empty. Existing files may be overwritten.')

## tmp dir
if os.path.isdir(tmpdir):
	if not tmpdir.endswith('/'):
		tmpdir = tmpdir + '/'
else:
	print('UTP Error: UTP cannot find tmp directory ' + tmpdir)
	sys.exit()

# creating unique tmp-dir for this run
time_stamp = timestamp('full')

tmpdir = str(tmpdir + 'UTP_tmpdir_' + time_stamp + '_tmp/')
os.mkdir(tmpdir)

## reading config (returns a dictionary with key(type) and value(arguments/path)
config = os.path.realpath(config)
if os.path.isfile(config):
	config_dict = config_reader(config)
elif os.path.isdir(config):
	config_dict = config_reader(str(config + 'utp.config'))
else:
	print('UTP Error: Cannot find config file at ' + config + '. Exiting...')
	sys.exit()

# list input files
if verbose:
	print('Read file(s) found:')
	for r in input_reads_files:
		print(r)
	print('Genome file(s) found:')
	for g in input_genome_files:
		print(g)
	print('Annotation file(s) found:')
	for a in input_annotation_files:
		print(a)
if len(input_genome_files) != len(input_annotation_files):
	print('UTP Error: number of genome files != number of annotation files, exiting.')
	sys.exit()



# creating empty lists and making necessary temporary directories
input_alignment_files = []

tmpdir_bam = tmpdir + 'bamfiles/'
tmpdir_bam_rmdpl = tmpdir + 'bamfiles_duplicates_removed/'
tmpdir_genome = tmpdir + 'genome_index/'
tmpdir_STAR = tmpdir + 'STAR_files/'
tmpdir_split = tmpdir + 'bamfiles_split/'

os.mkdir(tmpdir_bam)
os.mkdir(tmpdir_bam_rmdpl)
os.mkdir(tmpdir_genome)
if merge:
	os.mkdir(tmpdir_split)



# merge genome files (Python)
# merge genomes
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Merging/Converting submitted genomes')
input_genome_names = genome_merger(input_genome_files, input_annotation_files, tmpdir)

# convert annotation if necessary
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Merging/Converting submitted annotation')
input_annotation_files = glob.glob(str(tmpdir + '*.annotation.*'))
for anfile in input_annotation_files:
	if not anfile.endswith('.gtf'):
		if verbose:
			time_stamp = timestamp('time')
			print(time_stamp + ' Converting ' + anfile + ' to GTF format.')
		annot_name = os.path.basename(anfile)
		annot_name = tmpdir_genome + annot_name
		annot_name = remove_suffix(annot_name, '.bed')
		annot_name = remove_suffix(annot_name, '.gff')
		annot_name = remove_suffix(annot_name, '.gff3')
		annot_name = annot_name + '.gtf'
		gffread_cmd = str(config_dict['gffread_executable'] +
			' -o ' + annot_name +
			' --gtf ' +
			' -F ' +
			config_dict['gffread_arguments'] +
			' ' + anfile)
		os.system(gffread_cmd)
	else:
		os.system('cp ' + anfile + ' ' + tmpdir_genome)
# create list of annotation files in gtf format
input_annotation_files_gtf = glob.glob(str(tmpdir_genome + '*.gtf'))

# select correct files for alignment
if merge:
	input_genome = str(tmpdir + 'genomes_merged.fasta')
	input_annotation = str(tmpdir_genome + 'annotation_merged.gtf')
else:
	input_genome = glob.glob(str(tmpdir + '*.genome.fasta'))[0]
	input_annotation = input_annotation_files_gtf[0]

## make file and append annotation
if merge:
	os.system(str('touch ' + input_annotation))
	with open(input_annotation, 'a') as ip:
		for file in input_annotation_files_gtf:
			with open(file, 'r') as f:
				lines = f.readlines()
				for line in lines:
					if '#' not in line:
						ip.write(line)

# completing merging step
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Merging/Converting complete')



# perform alignment (STAR), remove duplicates (Picard), and perform QC (Picard)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Starting genome indexing')

## set STAR index command
STAR_index_cmd = str(config_dict['STAR_executable'] +
	' --runMode genomeGenerate' +
	' --runThreadN ' + str(threads) +
	' --genomeDir ' + tmpdir_genome +
	' --genomeFastaFiles ' + input_genome +
	' --sjdbGTFfile ' + input_annotation +
	' ' + config_dict['STAR_index_arguments'])
os.system(STAR_index_cmd)

if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Finished genome indexing. Starting alignment.')

## set empty lists to store filenames for downstream QC
bam_raw_list = []
bam_rmdpl_list = []

## quality control with Picard CollectAlignmentSummaryMetrics
## bam_raw_list = list of raw bam output from STAR
## bam_rmdpl_list = list of bam files with duplicates removed
## bam_final_list = bamfiles separated by genome
## making output directories for the QC
outdir_bamqc_raw = str(outdir + 'bamqc_raw/')
outdir_bamqc_rmdpl = str(outdir + 'bamqc_dupl-removed/')
outdir_bamqc_final = str(outdir + 'bamqc_final/')

## make directories if they don't exist
if not os.path.exists(outdir_bamqc_raw):
	os.mkdir(outdir_bamqc_raw)
if not os.path.exists(outdir_bamqc_rmdpl):
	os.mkdir(outdir_bamqc_rmdpl)
if not os.path.exists(outdir_bamqc_final):
	os.mkdir(outdir_bamqc_final)

## set Picard syntax arguments
if old_picard:
	picard_I = ' I='
	picard_O = ' O='
	picard_M = ' M='
	picard_R = ' R='
else:
	picard_I = ' -I '
	picard_O = ' -O '
	picard_M = ' -M '
	picard_R = ' -R '



## perform STAR alignment, remove duplicate reads, and perform QC
if single_end:
	for rfile in input_reads_files:
		basename_bam = os.path.basename(rfile)
		basename_bam = basename_bam.split('.')[0]
		tmpname_bam = tmpdir_bam + basename_bam #not + '.bam' as STAR requires an output prefix rather than a name
		strname_bam = tmpname_bam + 'Aligned.sortedByCoord.out.bam' #STAR-assigned name
		rmdname_bam = tmpdir_bam_rmdpl + basename_bam + '.bam'
		# collect bamnames for downstream QC
		bam_raw_list.append(strname_bam)
		bam_rmdpl_list.append(rmdname_bam)
		# add decompression flag if the reads are compressed
		if reads_compressed:
			read_decomp = ' --readFilesCommand ' + reads_decomp_method
		else:
			read_decomp = ' '
		# set STAR alignment command
		STAR_align_cmd = str(config_dict['STAR_executable'] +
			' --runThreadN ' + str(threads) +
			' --genomeDir ' + tmpdir_genome +
			' --readFilesIn ' + rfile +
			' --outTmpDir ' + tmpdir_STAR +
			' --outFileNamePrefix ' + tmpname_bam +
			' --outSAMtype BAM SortedByCoordinate ' +
			config_dict['STAR_align_arguments'] +
			' ' + read_decomp)
		# run STAR
		os.system(STAR_align_cmd)
		os.system('rm -r ' + tmpdir_STAR)
		# collect alignment metrics for the raw bam output from STAR
		bamfile_raw_basename = os.path.basename(strname_bam)
		bamfile_raw_out = outdir_bamqc_raw + bamfile_raw_basename + '.info.txt'
		picard_casm_raw_cmd = str(config_dict['picard_executable'] + ' CollectAlignmentSummaryMetrics' +
			picard_I + strname_bam +
			picard_O + bamfile_raw_out +
			picard_R + tmpdir + 'genomes_merged.fasta' +
			' ' + config_dict['picard_casm_arguments'])
		if verbose:
			print(picard_casm_raw_cmd)
		os.system(picard_casm_raw_cmd)
		# set Picard MarkDuplicates command
		picard_rmdpl_cmd = str(config_dict['picard_executable'] +
			' MarkDuplicates' +
			picard_I + strname_bam +
			picard_O + rmdname_bam + 
			picard_M + rmdname_bam + '.info.txt' +
			' ' + config_dict['picard_mdpl_arguments'])
		# run Picard MarkDuplicates
		os.system(picard_rmdpl_cmd)
		# collect alignment metrics for the bam files with duplicates removed
		bamfile_rmdpl_basename = os.path.basename(rmdname_bam)
		bamfile_rmdpl_out = outdir_bamqc_rmdpl + bamfile_rmdpl_basename + '.info.txt'
		picard_casm_rmdpl_cmd = str(config_dict['picard_executable'] + ' CollectAlignmentSummaryMetrics' +
			picard_I + rmdname_bam +
			picard_O + bamfile_rmdpl_out +
			picard_R + tmpdir + 'genomes_merged.fasta' +
			' ' + config_dict['picard_casm_arguments'])
		if verbose:
			print(picard_casm_rmdpl_cmd)
		os.system(picard_casm_rmdpl_cmd)
		# remove raw STAR output if tmp data not saved
		if not keep_tmp_data:
			if tmpcopy == 'no':
				os.system('rm ' + strname_bam)
		# index bam file
		samtools_index_rmdpl_cmd = str(config_dict['samtools_executable'] +
			' index' +
			' ' + rmdname_bam +
			' ' + config_dict['samtools_index_arguments'])
		os.system(samtools_index_rmdpl_cmd)
		# add filename to list
		input_alignment_files.append(rmdname_bam)
else:
	r1files = []
	r2files = []
	for rfile in input_reads_files:
		if 'R1' in rfile:
			r1files.append(rfile)
		elif 'R2' in rfile:
			r2files.append(rfile)
		else:
			print('UTP Error: Cannot find "R1" or "R2" in file name ' + rfile + ', removing this file from analysis.')
	# sorting the input lists
	r1files.sort()
	r2files.sort()
	if len(r1files) == len(r2files):
		index = range(0, (len(r1files)))
		for i in index:
			r1file = r1files[i]
			r2file = r2files[i]
			basename_bam = os.path.basename(r1file)
			basename_bam = basename_bam.split('.')[0]
			basename_bam = basename_bam.replace('R1','out')
			tmpname_bam = tmpdir_bam + basename_bam + '.' #not + '.bam' because STAR requires a prefix instead of a filename
			strname_bam = tmpname_bam + 'Aligned.sortedByCoord.out.bam' #STAR-assigned name
			rmdname_bam = tmpdir_bam_rmdpl + basename_bam + '.bam'
			# collect bamnames for downstream QC
			bam_raw_list.append(strname_bam)
			bam_rmdpl_list.append(rmdname_bam)
			# add decompression flag if the reads are compressed
			if reads_compressed:
				read_decomp = ' --readFilesCommand ' + reads_decomp_method
			else:
				read_decomp = ' '
			# set STAR alignment command
			STAR_align_cmd = str(config_dict['STAR_executable'] +
				' --runThreadN ' + str(threads) +
				' --genomeDir ' + tmpdir_genome +
				' --readFilesIn ' + r1file + ' ' + r2file +
				' --outTmpDir ' + tmpdir_STAR +
				' --outFileNamePrefix ' + tmpname_bam +
				' --outSAMtype BAM SortedByCoordinate ' +
 				config_dict['STAR_align_arguments'] +
				' ' + read_decomp)
			# run STAR
			os.system(STAR_align_cmd)
			os.system('rm -r ' + tmpdir_STAR)
			# collect alignment metrics for the raw bam output from STAR
			bamfile_raw_basename = os.path.basename(strname_bam)
			bamfile_raw_out = outdir_bamqc_raw + bamfile_raw_basename + '.info.txt'
			picard_casm_raw_cmd = str(config_dict['picard_executable'] + ' CollectAlignmentSummaryMetrics' +
				picard_I + strname_bam +
				picard_O + bamfile_raw_out +
				picard_R + tmpdir + 'genomes_merged.fasta' +
				' ' + config_dict['picard_casm_arguments'])
			if verbose:
				print(picard_casm_raw_cmd)
			os.system(picard_casm_raw_cmd)
			# set Picard MarkDuplicates command
			picard_rmdpl_cmd = str(config_dict['picard_executable'] +
				' MarkDuplicates' +
				picard_I + strname_bam +
				picard_O + rmdname_bam +
				picard_M + rmdname_bam + '.info.txt' +
				' ' + config_dict['picard_mdpl_arguments'])
			# run Picard MarkDuplicates
			os.system(picard_rmdpl_cmd)
			# collect alignment metrics for the bam files with duplicates removed
			bamfile_rmdpl_basename = os.path.basename(rmdname_bam)
			bamfile_rmdpl_out = outdir_bamqc_rmdpl + bamfile_rmdpl_basename + '.info.txt'
			picard_casm_rmdpl_cmd = str(config_dict['picard_executable'] + ' CollectAlignmentSummaryMetrics' +
				picard_I + rmdname_bam +
				picard_O + bamfile_rmdpl_out +
				picard_R + tmpdir + 'genomes_merged.fasta' +
				' ' + config_dict['picard_casm_arguments'])
			if verbose:
				print(picard_casm_rmdpl_cmd)
			os.system(picard_casm_rmdpl_cmd)
			# remove raw STAR output if not keep tmp data
			if not keep_tmp_data:
				if tmpcopy == 'no':
					os.system('rm ' + strname_bam)
			# index bam file
			samtools_index_rmdpl_cmd = str(config_dict['samtools_executable'] +
				' index' +
				' ' + rmdname_bam +
				' ' + config_dict['samtools_index_arguments'])
			os.system(samtools_index_rmdpl_cmd)
			# add filename to list
			input_alignment_files.append(rmdname_bam)

if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Finished read alignment')



# setting empty list to store bam names in
bam_final_list = []

# splitting bamfiles (if necessary)
no_genomes = len(input_genome_names)
number = 1
if merge:
	if verbose:
		time_stamp = timestamp('time')
		print(time_stamp + ' Started separating merged genome alignments')
	# finding the chromosomes that correspond to the right genome
	for genome_name in input_genome_names:
		genome_file_list = []
		genome_chrom_list = []
		with open(str(tmpdir_genome + 'chrNameLength.txt'), 'r') as genome_index_STAR:
			for line in genome_index_STAR:
				if genome_name in line:
					genome_file_list.append(line)
					gcl = line.split('\t')
					genome_chrom_list.append(gcl[0])
		# write genome file for re-heading
		genome_file = str(tmpdir + genome_name + '.genome')
		with open(genome_file, 'w') as gbf:
			for line in genome_file_list:
				gbf.write(line)
		# convert gtf to bam with bedtools bedtobam
		## find correct annotation file
		for an in input_annotation_files_gtf:
			if genome_name in an:
				annotation_file_gtf = an
				annotation_file_base = an.strip('.gtf')
				annotation_file_bam = annotation_file_base + '.bam'
				bedtools_bedtobam_cmd = str(config_dict['bedtools_executable'] + ' bedtobam' +
				' -i ' + annotation_file_gtf +
				' -g ' + genome_file +
				' > ' + annotation_file_bam)
				if verbose:
					print(bedtools_bedtobam_cmd)
				os.system(bedtools_bedtobam_cmd)
		# split bamfiles per scaffold/chromosome
		for bam_rmdpl in bam_rmdpl_list:
			## get the sample name
			bam_rmdpl_basename = os.path.basename(bam_rmdpl)
			bam_rmdpl_basename = bam_rmdpl_basename.split('.')[0]
			sam_split_name = str(tmpdir_split + bam_rmdpl_basename + '.' + genome_name + '.sam')
			bam_final_name = str(outdir + bam_rmdpl_basename + '.' + genome_name + '.bam')
			bam_final_list.append(bam_final_name)
			## print header of bam file made with bedtools bedtobam
			samtools_print_header_cmd = str(config_dict['samtools_executable'] + ' view' +
				' -H' +
				' ' + annotation_file_bam +
				' > ' + sam_split_name)
			## set samtools command to select reads belonging to the respective genome
			samtools_view_genome_select_cmd = str(config_dict['samtools_executable'] + ' view' +
				' ' + config_dict['samtools_view_arguments'] +
				' ' + bam_rmdpl +
				' ' + ' '.join(genome_chrom_list) +
				' >> ' + sam_split_name)
			## convert subsequent sam file to bam (outputs to pipe)
			samtools_samtobam_cmd = str(config_dict['samtools_executable'] + ' view' +
				' -b ' + sam_split_name)
			## set samtools sort command (input from pipe, output to pipe)
			samtools_sort_cmd = str(config_dict['samtools_executable'] + ' sort' +
				' -@ ' + str(int(threads - 1)) + # flag is to specify additional threads, not the total
				' -T ' + tmpdir + 'tmpbam' + # writes temporary files to the used tmpdir
				' ' + str(config_dict['samtools_sort_arguments']))
			## run commands
			if verbose:
				print(samtools_print_header_cmd)
			os.system(samtools_print_header_cmd)
			if verbose:
				print(samtools_view_genome_select_cmd)
			os.system(samtools_view_genome_select_cmd)
			if verbose:
				print(samtools_samtobam_cmd + ' | ' + samtools_sort_cmd + ' > ' + bam_final_name)
			os.system(samtools_samtobam_cmd + ' | ' + samtools_sort_cmd + ' > ' + bam_final_name)
			# delete bam and sam file if tmp data is deleted
			if number == no_genomes: # only deletes the bam file when on the last genome
				if not keep_tmp_data:
					if tmpcopy == 'no':
						os.system('rm ' + bam_rmdpl)
						os.system('rm ' + sam_split_name)
			# collect alignment metrics for the genome separated bam files
			bamfile_final_basename = os.path.basename(bam_final_name)
			bamfile_final_out = outdir_bamqc_final + bamfile_final_basename + '.info.txt'
			refname = bamfile_final_basename.split('.')[1] #second element of the name is the genome name
			refname = tmpdir + refname + '.genome.fasta'
			picard_casm_final_cmd = str(config_dict['picard_executable'] + ' CollectAlignmentSummaryMetrics' +
				picard_I + bam_final_name +
				picard_O + bamfile_final_out +
				picard_R + refname +
				' ' + config_dict['picard_casm_arguments'])
			if verbose:
				print(picard_casm_final_cmd)
			os.system(picard_casm_final_cmd)
		# increase number to track if it is on the last genome
		number += 1
	if verbose:
		time_stamp = timestamp('time')
		print(time_stamp + ' Finished seperating alignment per genome')
else:
	for rmd_bam in input_alignment_files:
		os.system('mv ' + rmd_bam + ' ' + outdir)
		rmd_bam_basename = os.path.basename(rmd_bam)
		bam_final_list.append(str(outdir + rmd_bam_basename))



# index bam files
samtools_index_cmd = str(
	'for file in ' + outdir + '*.bam ; do ' +
	config_dict['samtools_executable'] + ' index ' +
	config_dict['samtools_index_arguments'] +
	' $file ; done')
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Indexing output BAM files')
os.system(samtools_index_cmd)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Finished indexing of output BAM files')



# perform multiQC on Picard CollectAlignmentSummaryMetrics output
if merge:
	if verbose:
		time_stamp = timestamp('time')
		print(time_stamp + ' Running MultiQC to summarize alignment quality')
	## set filenames for reports
	multiqc_raw_name = outdir_bamqc_raw + 'bamqc_raw.summary.txt'
	multiqc_rmdpl_name = outdir_bamqc_rmdpl + 'bamqc_rmdpl.summary.txt'
	multiqc_final_name = outdir_bamqc_final + 'bamqc_final.summary.txt'
	## set commands for multiqc
	multiqc_raw_cmd = str(config_dict['multiqc_executable'] + ' ' + config_dict['multiqc_arguments'] +
		' -n ' +  multiqc_raw_name + ' ' + outdir_bamqc_raw)
	multiqc_rmdpl_cmd = str(config_dict['multiqc_executable'] + ' ' + config_dict['multiqc_arguments'] +
		' -n ' +  multiqc_rmdpl_name + ' ' + outdir_bamqc_rmdpl)
	multiqc_final_cmd = str(config_dict['multiqc_executable'] + ' ' + config_dict['multiqc_arguments'] +
		' -n ' +  multiqc_final_name + ' ' + outdir_bamqc_final)
	## run multiqc
	if verbose:
		print(multiqc_raw_cmd)
	os.system(multiqc_raw_cmd)
	if verbose:
		print(multiqc_rmdpl_cmd)
	os.system(multiqc_rmdpl_cmd)
	if verbose:
		print(multiqc_final_cmd)
	os.system(multiqc_final_cmd)
	if verbose:
		time_stamp = timestamp('time')
		print(time_stamp + ' Finished MultiQC')
else:
	if verbose:
		time_stamp = timestamp('time')
		print(time_stamp + ' Running MultiQC to summarize alignment quality')
	## set filenames for reports
	multiqc_raw_name = outdir_bamqc_raw + 'bamqc_raw.summary.txt'
	multiqc_rmdpl_name = outdir_bamqc_rmdpl + 'bamqc_rmdpl.summary.txt'
	## set commands for multiqc
	multiqc_raw_cmd = str(config_dict['multiqc_executable'] + ' ' + config_dict['multiqc_arguments'] +
		' -n ' +  multiqc_raw_name + ' ' + outdir_bamqc_raw)
	multiqc_rmdpl_cmd = str(config_dict['multiqc_executable'] + ' ' + config_dict['multiqc_arguments'] +
		' -n ' +  multiqc_rmdpl_name + ' ' + outdir_bamqc_rmdpl)
	## run multiqc
	if verbose:
		print(multiqc_raw_cmd)
	os.system(multiqc_raw_cmd)
	if verbose:
		print(multiqc_rmdpl_cmd)
	os.system(multiqc_rmdpl_cmd)
	if verbose:
		time_stamp = timestamp('time')
		print(time_stamp + ' Finished MultiQC')



# running featureCounts
## if merged, use list input_annotation_files_gtf
## if not, use single file input_annotation
## set if the reads are paired-end or single-end
if single_end:
	read_end = ' '
else:
	read_end = ' -p '
## construct featureCounts command
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Running featureCounts to generate gene expression counts')
if merge:
	for annot_file_cnt in input_annotation_files_gtf:
		annot_file_basename = os.path.basename(annot_file_cnt)
		annot_file_basename = annot_file_basename.split('.')[0]
		# select alignment files per genome
		cnt_subset = []
		for bam_file_cnt in bam_final_list:
			if annot_file_basename in bam_file_cnt.split('.')[1]:
				cnt_subset.append(bam_file_cnt)
		# run featurecounts on the subset of bamfiles
		featurecounts_cmd = str(config_dict['featurecounts_executable'] +
			read_end +
			' ' + config_dict['featurecounts_arguments'] +
			' -f' +
			' -T ' + str(threads) +
			' --tmpDir ' + tmpdir +
			' -a ' + annot_file_cnt +
			' -o ' + outcnt + annot_file_basename + '.counts ' +
			' '.join(cnt_subset))
		if verbose:
			print(featurecounts_cmd)
		os.system(featurecounts_cmd)
else:
	annot_file_basename = os.path.basename(input_annotation)
	annot_file_basename = annot_file_basename.split('.')[0]
	# run featurecounts on the bamfiles
	featurecounts_cmd = str(config_dict['featurecounts_executable'] +
		read_end +
		' ' + config_dict['featurecounts_arguments'] +
		' -T ' + str(threads) +
		' --tmpDir ' + tmpdir +
		' -a ' + input_annotation +
		' -o ' + outcnt + annot_file_basename + '.counts ' +
		' '.join(bam_final_list))
	if verbose:
		print(featurecounts_cmd)
	os.system(featurecounts_cmd)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Finished featureCounts')



# running multiQC on the featureCounts summaries
multiqc_featureCounts_cmd = str(config_dict['multiqc_executable'] + ' ' + config_dict['multiqc_arguments'] +
	' -n ' + outcnt + 'featureCounts_summary ' + outcnt)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Running MultiQC on featureCounts summary reports')
	print(multiqc_featureCounts_cmd)
os.system(multiqc_featureCounts_cmd)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Finished MultiQC')



# handle tmp directory data
if not tmpcopy == 'no':
	os.system('mv ' + tmpdir + ' ' + tmpcopy)
elif not keep_tmp_data:
	time_stamp = timestamp('time')
	print(time_stamp + ' Cleaning temporary directory')
	os.system('rm -rf ' + tmpdir)

# exit pipeline
time_stamp = timestamp('time')
print(time_stamp + ' UTP run complete')
