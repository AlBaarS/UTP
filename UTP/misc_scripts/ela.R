#!/usr/bin/env Rscript

# Script to analyze expression levels of genes over a dataset (Expression Level Analyzer)
# Made by Alejandro Baars as part of his General Research Profile
# Utrecht University 2023

# INITIALIZING
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
pacman::p_load(tidyverse, Rsubread, limma, edgeR, optparse, update = FALSE)



# SET ARGUMENTS
ela_options <- list(
	make_option(c('--input','-i'),
		type = 'character',
		default = '',
		help = '<Required> Input for your counts. For multiple files; separate with commas (path/to/file1,path/to/file2,path/to/fileN).'),
	make_option(c('--output','-o'),
		type = 'character',
		default = '',
		help = '<Required> Specify the output path and prefix.'),
	make_option('--conditions',
		type = 'character',
		default = '',
		help = '<Required> Specify the conditions as names in a list (cond1,cond2,condN). Number does not need to be equal to the number of files submitted.'),
	make_option('--replicates',
		type = 'character',
		default = '',
		help = '<Required> Specify the number of replicates per condition (4,3,4). The sum of this list needs to be equal to the total number of columns in all files.'),
	make_option('--mode',
		type = 'character',
		default = 'testing',
		help = paste0('Specify what the script should do. With "testing" (default) the script will show you various metrics of the submitted data. With',
			' "filter" the script will make a txt file with geneIDs of genes that fail to meet the submitted threshold.')),
	make_option('--threshold',
		type = 'numeric',
		default = 1,
		help = 'Specify the threshold of expression when mode == "filter". Default = 1'),
	make_option('--normalization',
		type = 'character',
		default = 'RPKM',
		help = 'Specify the normalization method of your data. Accepts either RPKM (default) or CPM.'),
	make_option('--samples',
		type = 'character',
		default = '',
		help = paste0('Specify the sample names in same order as they appear in the counts files and in the order you specified the ',
			'counts files, seperated by commas.'))
)

# PROCESS ARGUMENTS
# Parse argument list
ela_arguments <- OptionParser(option_list=ela_options)
ela_parsed <- parse_args(ela_arguments)

# Check if required arguments are present
if (ela_parsed$input == '') {
	stop('Input parameter missing. Exiting...')
}

if (ela_parsed$output == '') {
	stop('Output directory parameter missing. Exiting...')
}

if (ela_parsed$replicates == '') {
	stop('Replicates parameter missing. Exiting...')
}

if (ela_parsed$conditions == '') {
	stop('Conditions parameter missing. Exiting...')
}

# Store as variables
input_files <- str_split_1(ela_parsed$input, ',')
output <- ela_parsed$output
replicates <- as.numeric(str_split_1(ela_parsed$replicates, ','))
conditions <- str_split_1(ela_parsed$conditions, ',')
mode <- ela_parsed$mode
threshold <- ela_parsed$threshold
norm <- ela_parsed$normalization
samples <- str_split_1(ela_parsed$samples, ',')

basename <- tail(str_split_1(output, '/'), 1)



# FUNCTIONS
# function to read the counts data from one file into a table
readCounts <- function(counts, sample_names) {

	# import counts file
	input_counts <- read_tsv(counts, col_names = TRUE, comment = '#', show_col_types = FALSE)

	# prepare to reshape header (for shorter names)
	input_counts_header_old <- base::colnames(input_counts)
	input_counts_header_new <- vector('list', (length(input_counts_header_old) - 5))

	# store gene lengths in vector
	input_counts_gene_length <- pull(input_counts, Length)

	# filter the necessary column names
	for (i in 1:length(input_counts_header_old)) {
		if (i == 1) {
			input_counts_header_new[i] <- input_counts_header_old[i]
		} else if (i >= 7) {
			if (sample_names[[1]] != '') {	# check if sample names are specified and use those if available (can also be done by merger function)
				input_counts_header_new[(i - 5)] <- sample_names[(i - 6)]
			} else {	# use manual sample naming, which just uses the filename minus the path
				sample_split <- str_split_1(input_counts_header_old[i], '/')
				sample_name <- paste(sample_split[length(sample_split)])
				input_counts_header_new[(i - 5)] <- sample_name
			}
		}
	}

	# reshape and rename data
	input_counts_tidy <- dplyr::select(input_counts, 1, 7:length(input_counts))
	base::colnames(input_counts_tidy) <- input_counts_header_new

	return(list(input_counts_tidy, input_counts_gene_length))
}	

# function to merge multiple counts tables into a single table. Assumes that all submitted tables have the same genes/regions.
mergeCounts <- function(counts, sample_names) {

	# read the counts into tables and put these into a list
	input_counts_list <- list()
	for (i in 1:length(counts)) {
		input_counts_tmp <- readCounts(counts[i], '')	# use the earlier defined function to read the data
		input_counts <- input_counts_tmp[[1]]
		input_counts_gene_length <- c(input_counts_tmp[[2]])
		input_counts_list <- c(input_counts_list, list(input_counts))
	}
	# merge the obtained counts tables
	for (c in 1:length(counts)) {
		if (c == 1) {
			input_counts_merged <- input_counts_list[[c]]
		} else {
			tmp_df <- input_counts_list[[c]]
			if (all(input_counts_merged$Geneid == tmp_df$Geneid)) {
				input_counts_merged <- cbind(input_counts_merged, dplyr::select(tmp_df, -Geneid))
			}
		}
	}
	# use manual sample names
	if (sample_names[[1]] != '') {
		base::colnames(input_counts_merged) <- c('Geneid', sample_names)
	}

	return(list(input_counts_merged, input_counts_gene_length))
}

# Function to normalize counts
normalizeCounts <- function(input_counts, input_counts_gene_length, conditions_fct, replicates_fct, norm_library = 'TMMwsp', norm_sample = 'RPKM') {

	# creating normalized expression tibble
	## normalize using RPKM (which takes gene length into account)
	if (norm_sample == 'RPKM'){
		counts_norm <- input_counts %>% 
			DGEList(., genes=.$Geneid, group=conditions_fct) %>%
			normLibSizes(method = norm_library, ) %>%
			rpkm(gene.length = input_counts_gene_length) %>%
			as_tibble %>%
			cbind(input_counts$Geneid, .)
		base::colnames(counts_norm)[1] <- 'Geneid'
	## normalize using CPM
	} else if (norm_sample == 'CPM'){
		counts_norm <- input_counts %>% 
			DGEList(., genes=.$Geneid, group=conditions_fct) %>%
			normLibSizes(method = norm_library, ) %>%
			cpm() %>%
			as_tibble %>%
			cbind(input_counts$Geneid, .)
		base::colnames(counts_norm)[1] <- 'Geneid'
	} else {
		break(paste('Submitted normalization method "', norm_sample, '" not recognized, use "CPM" or "RPKM" (capitalized)', sep = ''))
	}

	return(counts_norm)
}

# Function to generate plots
plotRPKM <- function(counts, condition) {

	# Create array to filter the gene rowmeans by
	RPKM_array <- seq(0, 10, by = 0.05)

	# Infer column names from table
	cnames <- colnames(counts)[2:length(colnames(counts))]

	# Filter genes based on array
	expr_filtered <- tibble()
	colnames(expr_filtered) <- c('RPKM', cnames)
	for (a in RPKM_array) {
		no_genes_list <- vector('list', (length(cnames)))
		for (n in 1:length(cnames)) {
			no_genes <- nrow(filter(counts, counts[[cnames[[n]]]] >= a))
			no_genes_list[[n]] <- no_genes
		}
		no_genes_list <- append(no_genes_list, a, 0)
		names(no_genes_list) <- c('RPKM', cnames)
		expr_filtered <- rbind(expr_filtered, no_genes_list)
	}

	# Write table to file
	table_file <- paste0(output, '_', condition, '_table.csv')
	write_csv(expr_filtered, file = table_file)

	# Re-shape data for plotting
	expr_gathered <- gather(expr_filtered, 'sample', 'no_genes', 2:ncol(expr_filtered))

	# Plot results
	## Create plot
	array_plot <- ggplot(
		data = expr_gathered, 
		aes(
			x = RPKM,
			y = no_genes,
			color = sample
		)
	) +
	geom_point() +
	theme_classic() +
	labs(
		x = paste0(norm),
		y = paste0('Number of genes at ', norm, ' cutoff'),
		title = paste0(condition, ' ', basename, ' - Number of genes at ', norm, ' cutoff plotted over ', norm),
		subtitle = paste0('N = ', nrow(counts_normalized_df))
	)

	## Create filename
	plot_file <- paste0(output, '_', condition, '_plot.png')

	## Save plot
	ggsave(plot_file, plot = array_plot)

}


# DATA PREPERATION
## sort through replicates and conditions
conditions_fct <- c()
replicates_fct <- c()

for (r in 1:length(replicates)) {
	conditions_fct <- c(conditions_fct, rep(conditions[[r]], each=replicates[[r]]))
	replicates_fct <- c(replicates_fct, 1:replicates[[r]])
}

conditions_fct <- as.factor(conditions_fct)
replicates_fct <- as.factor(replicates_fct)

## select columns per condition
indices <- vector('list', length(replicates))

for (i in 1:length(replicates)) {
	index <- c((1:replicates[[i]]) + sum(replicates[(1:i - 1)]) + 1)
	indices[[i]] <- index
}

#print(indices) #tmpprint



# NORMALIZE INPUT AND STORE IN DATAFRAME
if (length(input_files) > 1) {
	# first we merge the count tables
	input_data <- mergeCounts(input_files, samples)
	input_table <- input_data[[1]]
	input_counts_gene_length <- input_data[[2]]
	# then we normalize everything
	counts_normalized_df <- normalizeCounts(input_table, input_counts_gene_length, 
		conditions_fct, replicates_fct, norm_sample = norm)
#	print('head(counts_normalized_df)') #tmpprint
#	print(head(counts_normalized_df)) #tmpprint
#	print(str(counts_normalized_df)) #tmpprint
} else {
	# no tables need to be merged so we normalize immediately
	input_data <- readCounts(input_files[[1]], samples)
	input_table <- input_data[[1]]
	input_counts_gene_length <- input_data[[2]]
	counts_normalized_df <- normalizeCounts(input_table, input_counts_gene_length, 
		conditions_fct, replicates_fct, norm_sample = norm)
#	print('head(counts_normalized_df)') #tmpprint
#	print(head(counts_normalized_df)) #tmpprint
#	print(str(counts_normalized_df)) #tmpprint
}

# Create table with row means per condition & row mean of all rows
counts_rowmean_full_df <- select(counts_normalized_df, Geneid)
for (ci in 1:length(conditions)) {
	counts_rowmean_full_df[conditions[[ci]]] <- rowMeans(select(counts_normalized_df, indices[[ci]]))
}

print(head(counts_rowmean_full_df)) #tmpprint

counts_rowmean <- rowMeans(select(counts_normalized_df, -Geneid))
counts_rowmean_df <- tibble(Geneid = counts_normalized_df$Geneid, rowmean = counts_rowmean)



# ANALYSIS
if (mode == 'testing') {	

	# Filter genes with no expression in any sample
	no_noExpr <- nrow(filter(counts_rowmean_df, rowmean == 0))
	print(paste0('Number of genes with 0 expression in all samples: ', as.character(no_noExpr)))

	for (cond in 1:length(conditions)) {									# create a plot for each condition
		table_tmp <- select(counts_normalized_df, Geneid, indices[[cond]])	# this way, samples can be compared to each-other
		#print(head(table_tmp)) #tmpprint									# and we get a good overview of the full data
		plotRPKM(table_tmp, conditions[[cond]])
	}
	plotRPKM(counts_rowmean_full_df, 'full')								# run the plotting function on the summarized data

} else if (mode == 'filter') {

	# select all genes that have mean RPKMs below the threshold in all conditions
	genes_rm_df <- counts_rowmean_full_df %>% filter(if_all(!contains('Geneid'), ~ . < threshold))
	print(paste0('Genes below threshold: ', nrow(genes_rm_df)))
	genes_rm <- select(genes_rm_df, Geneid)

	# write genes that need to be omitted from analysis to a txt file
	genes_rm_file <- paste0(output, '_genes_rm.txt')
	write_tsv(genes_rm, file = genes_rm_file, col_names = FALSE)

} else {

	stop(paste0('Unrecognized mode argument: "', mode, '". Use "testing" or "filter".'))

}
