#!/usr/bin/env python3

# Cleaner script, part of the uniform transciptomics pipeline
# Made by Alejandro Baars as part of his General Research Profile
# Utrecht University 2023

# importing packages
import os, sys, glob, argparse as ap, datetime

FqC = ap.ArgumentParser(
	prog="Fastq Cleaner script",
	usage="FqC.py --reads [READS] --outdir [OUTDIR] --tmpdir [TMPDIR] [arguments]",
	description="Pipeline to clean and perform quality control on Illumina sequencing reads.")

# arguments
FqC.add_argument('--version','-v',
	action = 'version',
	version = 'FqC v0.1.1')
FqC.add_argument('--verbose','-V',
	action = 'store_true',
	default = False,
	help = 'If specified, prints messages as to what it is currently doing. Off by default.')
FqC.add_argument('--reads','-r',
	nargs = '+',
	action = 'store',
	required = True,
	help = '<Required> Input for your reads. Accepts a file, list of files, directory, or list of directories.')
FqC.add_argument('--threads','-t',
	action = 'store',
	default = 2,
	help = 'Specify the number of threads used by programmes that can use multithreading. Default = 2.')
FqC.add_argument('--single_end',
	action = 'store_true',
	default = False,
	help = 'The pipeline assumes the reads are paired-end. If you use single-end reads, include this tag.')
FqC.add_argument('--outdir','-o',
	action = 'store',
	required = True,
	help = '<Required> Specify your output directory. All your output will be stored here.')
FqC.add_argument('--config','-c',
	action = 'store',
	default = ' ',
	help = 'Specify a (customized) config file to use rather than the default one.')
FqC.add_argument('--tmpdir',
	action = 'store',
	required = True,
	help = '<Required> Specify a directory for intermediate files.')
FqC.add_argument('--keep_tmp_data',
	action = 'store_true',
	default = False,
	help = 'Keep the data written to the specified temporary directory. By default, this directory is cleared.')
FqC.add_argument('--copy_tmp_data',
	action = 'store',
	default = 'no',
	help = 'Copy the data written in the specified temporary directory to a directory specified here instead of deleting it.')

# process arguments
namespace = FqC.parse_args()

## input data
input_reads = list(namespace.reads)

## directories & configuration
outdir = str(namespace.outdir)
config = str(namespace.config)
tmpdir = str(namespace.tmpdir)
threads = int(namespace.threads)

## internal optional arguments
verbose = bool(namespace.verbose)
single_end = bool(namespace.single_end)
keep_tmp_data = bool(namespace.keep_tmp_data)
tmpcopy = str(namespace.copy_tmp_data)

## config
if config == ' ':
	# tmpprint
	print('config = None')
	utp_path = os.path.dirname(os.path.realpath(__file__))
	utp_path = utp_path.strip('/bin/')
	config = '/' + utp_path + '/config/utp.config'
	# tmpprint
	print(utp_path)
	print(config)



# set functions

## this function disentangles the various methods to provide input for the reads and genomes
## this function is the same as in the UTP script, so it does contain some redundant code
def input_untangler(input):
	output = []
	extensions = ['.fastq','.fastq.gz','.fastq.zip','.fq','.fq.gz','.fq.zip','.fastqsanger','fastqsanger.gz','fastqsanger.zip']
	# check if the input consists of 1 entry or a list of entries
	if len(input) == 1:
		# check if the input is a directory or a file
		if os.path.isdir(input[0]):
			dir = input[0]
			for file in os.listdir(dir):
				if file.endswith(tuple(extensions)):
					file_out = os.path.join(dir, file)
					output.append(file_out)
		elif os.path.isfile(input[0]):
			file = input[0]
			if str(file).endswith(tuple(extensions)):
				output.append(file)
		else:
			print('FqC Error: Files not found; aborting.')
			sys.exit()
	elif len(input) > 1:
		for i in input:
			if os.path.isdir(i):
				dir = i
				for file in os.listdir(dir):
					if file.endswith(tuple(extensions)):
						file_out = os.path.join(dir, file)
						output.append(file_out)
			elif os.path.isfile(i):
				file = i
				if file.endswith(tuple(extensions)):
					output.append(file)
	# returning a properly formatted list of all input files with path
	return output

## this function reads the config file and stores the options for each programme in a dictionary
def config_reader(file):
	if os.path.isfile(file):
		cfile = {}
		file = os.path.realpath(file)
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

# pipeline start
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' FqC initialized')

# input processing

## neatly organize the input reads
input_reads_files = input_untangler(input_reads)

## check if the submitted directories exists
if os.path.isdir(outdir):
	if not outdir.endswith('/'):
		outdir = outdir + '/'
else:
	time_stamp = timestamp('time')
	print(time_stamp + ' FqC Error: Cannot find output directory.')
	sys.exit()

if os.listdir(outdir) != []:
	print('FqC Warning: Output directory not empty. Existing files may be overwritten.')

if os.path.isdir(tmpdir):
	if not tmpdir.endswith('/'):
		tmpdir = tmpdir + '/'
else:
	print('FqC Error: Cannot find temporary directory')
	os.system('echo $(date -u) "Exiting"')
	sys.exit()

# creating unique tmp-dir for this run
time_stamp = timestamp('full')
tmpdir = str(tmpdir + 'FqC_tmpdir_' + time_stamp + '_tmp/')
os.system('mkdir ' + tmpdir)

## reading config (returns a dictionary with key(type) and value(arguments/path)
if os.path.isfile(config):
	config_dict = config_reader(config)
elif os.path.isdir(config):
	config_dict = config_reader(str(config + 'utp.config'))
else:
	print('FqC Error: Cannot find config file at ' + config + '. Exiting...')
	sys.exit()

# development/debugging
if verbose:
	print('Read file(s) found:')
	for r in input_reads_files:
		print(r)



# perform read quality control on untrimmed reads (Fastqc)
## making necessary directories
tmpdir_fqc = str(tmpdir + 'fastqc')
if not os.path.isdir(tmpdir_fqc):
	os.system(str('mkdir ' + tmpdir_fqc))
outdir_fqc_untrimmed = str(outdir + 'fastqc_untrimmed/')
if not os.path.isdir(outdir_fqc_untrimmed):
	os.system(str('mkdir ' + outdir_fqc_untrimmed))

## running fastqc on the untrimmed files
reads_untrimmed = ' '.join(input_reads_files)
fastqc_in_command = str(config_dict['fastqc_executable'] + ' -o ' + outdir_fqc_untrimmed +
	' -d ' + tmpdir_fqc + ' -t ' + str(threads) + ' -noextract ' + 
	config_dict['fastqc_arguments'] + ' ' + reads_untrimmed)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' FastQC initiated on input reads')
	print(fastqc_in_command)
os.system(fastqc_in_command)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' FastQC completed')

## running multiqc on the fastqc results to summarize the fastqc results
report_in_name = str(outdir_fqc_untrimmed + '.fastqc_summary.untrimmed.txt')
multiqc_in_command = str(config_dict['multiqc_executable'] + ' ' + config_dict['multiqc_arguments'] + ' -n ' +
	report_in_name + ' ' + outdir_fqc_untrimmed)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' MultiQC summarizing of FastQC initiated')
	print(multiqc_in_command)
os.system(multiqc_in_command)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' MultiQC completed')

# trim adapters of reads (trimmomatic)
input_reads_trimmed = []
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Trimming reads with Trimmomatic')
if single_end:
	# single-end reads do not need to be paired and thus can just be run 1-by-1
	for rfile in input_reads_files:
		routname = os.path.basename(rfile)
		routname = routname.split('.')[0]
		routname = str(outdir + routname + '.trimmed.fq.gz')
		input_reads_trimmed.append(routname)
		trimmomatic_command = str(config_dict['trimmomatic_executable'] + ' SE -threads ' + str(threads) + ' ' +
			config_dict['trimmomatic_phred'] + ' -summary ' + routname + '.summary.txt ' + 
			rfile + ' ' + routname + ' ' + config_dict['trimmomatic_arguments'])
		if verbose:
			time_stamp = timestamp('time')
			print(time_stamp + ' Running Trimmomatic:')
			print(trimmomatic_command)
		os.system(trimmomatic_command)
else:
	# paired-end reads need to be paired properly before running the trimmer on it.
	r1files = []
	r2files = []
	for rfile in input_reads_files:
		if 'R1' in rfile:
			r1files.append(rfile)
		elif 'R2' in rfile:
			r2files.append(rfile)
		else:
			print('FqC Error: Cannot find "R1" or "R2" in file name ' + rfile + ', removing this file from analysis.')
	# sort the input lists because they are mismatching for some reason
	r1files.sort()
	r2files.sort()
	if len(r1files) == len(r2files):
		index = range(0, (len(r1files)))
		for i in index:
			# preparing the necessary names
			goutname = os.path.basename(r1files[i])
			goutname = goutname.split('.')[0]
			goutname = goutname.replace('R1','out')
			r1outname = os.path.basename(r1files[i])
			r1outname = r1outname.split('.')[0]
			u1outname = str(outdir + r1outname + '.trimmed.unpaired.fq.gz')
			r1outname = str(outdir + r1outname + '.trimmed.paired.fq.gz')
			r2outname = os.path.basename(r2files[i])
			r2outname = r2outname.split('.')[0]
			u2outname = str(outdir + r2outname + '.trimmed.unpaired.fq.gz')
			r2outname = str(outdir + r2outname + '.trimmed.paired.fq.gz')
			input_reads_trimmed.append(r1outname)
			input_reads_trimmed.append(r2outname)
			# setting trimming command
			trimmomatic_command = str(config_dict['trimmomatic_executable'] + ' PE -threads ' + str(threads) + ' ' +
				config_dict['trimmomatic_phred'] + ' -summary ' + tmpdir + goutname + '.summary.txt ' +
				r1files[i] + ' ' + r2files[i] + ' ' + 
				r1outname + ' ' + u1outname + ' ' +
				r2outname + ' ' + u2outname + ' ' +
				config_dict['trimmomatic_arguments'])
			if verbose:
				time_stamp = timestamp('time')
				print(time_stamp + ' Running Trimmomatic:')
				print(trimmomatic_command)
			os.system(trimmomatic_command)
	else:
		print('FqC Error: unequal number of paired read files.')
		time_stamp = timestamp('time')
		print(time_stamp + ' Exiting')
		sys.exit()

if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Finished read trimming')

# perform read quality control on trimmed reads (Fastqc)
## filepaths:
outdir_fqc_trimmed = str(outdir + 'fastqc_trimmed/')
if not os.path.isdir(outdir_fqc_trimmed):
	os.system(str('mkdir ' + outdir_fqc_trimmed))

## running fastqc on each file
reads_trimmed = ' '.join(input_reads_trimmed)
fastq_out_command = str(config_dict['fastqc_executable'] + ' -o ' + outdir_fqc_trimmed +
	' -d ' + tmpdir_fqc + ' -t ' + str(threads) + ' -noextract ' +
	config_dict['fastqc_arguments'] + ' ' + reads_trimmed)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' Running FastQC on trimmed reads')
	print(fastq_out_command)
os.system(fastq_out_command)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' FastQC complete')

## running multiqc on the fastqc results to summarize the fastqc results
report_out_name = str(outdir_fqc_trimmed + '.fastqc_summary.trimmed.txt')
multiqc_out_command = str(config_dict['multiqc_executable'] + ' ' + config_dict['multiqc_arguments'] + ' -n ' +
	report_out_name + ' ' + outdir_fqc_trimmed)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' MultiQC summarizing of FastQC initiated')
	print(multiqc_out_command)
os.system(multiqc_out_command)
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' MultiQC completed')

# cleaning up temporary directory
if not tmpcopy == 'no':
	os.system('mv ' + tmpdir + ' ' + tmpcopy)
elif not keep_tmp_data:
	time_stamp = timestamp('time')
	print(time_stamp + ' Cleaning temporary directory')
	os.system('rm -rf ' + tmpdir)

# end of script
if verbose:
	time_stamp = timestamp('time')
	print(time_stamp + ' FqC run completed')
