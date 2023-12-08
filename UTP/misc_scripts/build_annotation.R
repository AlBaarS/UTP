#!/usr/bin/env Rscript

# Script to build and import an annotation package for GO-term enrichment analysis
# Made by Alejandro Baars as part of his General Research Profile
# Utrecht University 2023

# checking and installing dependencies
if (!require('pacman')) {
	install.packages('pacman')
}

if (!require('BiocManager')) {
	install.packages('BiocManager')
}

# loading pacman and BiocManager, packages used to install the rest
library('pacman')
library('BiocManager')

# use pacman to check and install the remaining packages
pacman::p_load(tidyverse, Rsubread, optparse, AnnotationDbi, AnnotationForge, update = FALSE)

# packages that don't work with pacman
if (!require('GOstats', quietly = TRUE))
	BiocManager::install('GOstats')
if (!require('AnnotationHub', quietly = TRUE))
	BiocManager::install('AnnotationHub')

library('GOstats')
library('AnnotationHub')

# set script options
buildAnn_options <- list(
make_option(c('--input','-i'),
	type = 'character',
	default = '',
	help = 'Specify the input tables here in the order GO_table.csv,sym_table.csv,chrom_table.csv,metadata.txt (required)'),
make_option(c('--outdir','-o'),
	type = 'character',
	default = '',
	help = 'Specify the output directory here (required)'),
make_option('--skip_lines',
	type = 'numeric',
	default = 0,
	help = 'Specify the number of lines you want to skip (e.g. for a header). Default = 0'))

# parse input arguments
buildAnn_arguments <- OptionParser(option_list=buildAnn_options)
buildAnn_parsed <- parse_args(buildAnn_arguments)

# CHECK ESSENTIAL ARGUMENTS
if (buildAnn_parsed$input == '') {
	stop('Input parameter "input" is missing. Exiting...')
}

if (buildAnn_parsed$outdir == '') {
	stop('Input parameter "outdir" is missing. Exiting...')
}

# PROCESS ARGUMENTS
input_files <- str_split_1(buildAnn_parsed$input, ',')
outdir <- buildAnn_parsed$outdir
skip_lines <- buildAnn_parsed$skip_lines

## assign individual files to variables
input_file_GO <- input_files[[1]]
input_file_sym <- input_files[[2]]
input_file_chrom <- input_files[[3]]
input_file_meta <- input_files[[4]]



# IMPORT DATA
## GO-terms
input_GO <- as.data.frame(read_csv(input_file_GO, col_names = FALSE, skip = skip_lines))
colnames(input_GO) <- c('GID', 'GO', 'EVIDENCE')
input_GO <- dplyr::distinct(input_GO)
print('GO-term table:')
print(head(input_GO)) # check if everything looks OK

## Gene symbols
input_sym <- as.data.frame(read_csv(input_file_sym, col_names = FALSE, skip = skip_lines))
colnames(input_sym) <- c('GID','SYMBOL','GENENAME')
input_sym <- dplyr::distinct(input_sym)
print('Symbol table:')
print(head(input_sym))

## Chromosome table
input_chrom <- as.data.frame(read_csv(input_file_chrom, col_names = FALSE, skip = skip_lines))
colnames(input_chrom) <- c('GID', 'CHROMOSOME')
input_chrom <- dplyr::distinct(input_chrom)
print('Chromosome table:')
print(head(input_chrom))

## Metadata
source(input_file_meta)



# BUILD PACKAGE
## variables version, maintainer, author, tax_id, genus, and species come from the sourced metadata file
makeOrgPackage(
	gene_info = input_sym,
	chromosome = input_chrom,
	go = input_GO,
	version = version,
	maintainer = maintainer,
	author = author,
	outputDir = outdir,
	tax_id = tax_id,
	genus = genus,
	species = species,
	goTable = 'go')



# EPILOGUE
print(paste('Package created successfully for', genus, species, 'in', outdir))
