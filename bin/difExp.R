#!/usr/bin/env Rscript

# Script to calculate differential expression from generated expression counts
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
pacman::p_load(tidyverse, limma, edgeR, statmod, corrplot, optparse, here, ggrepel, ggpointdensity, update = FALSE)

# version specification
difExp_version <- '0.1.7'
print(paste('difExp.R version:', difExp_version))

# SET ARGUMENTS
difExp_options <- list(
	make_option(c('--input_counts', '-i'),
		type = 'character', 
		default = '',
		help = 'Specify the input count files here, separated by commas (required).'),
	make_option(c('--output_dir', '-o'),
		type = 'character',
		default = '',
		help = 'Specify your output directory here (required).'),
	make_option(c('--name', '-n'),
		type = 'character',
		default = 'DE_analysis',
		help = 'Specify the name of your output files (default: DE_analysis)'),
	make_option(c('--replicates', '-r'),
		type = 'integer',
		default = '',
		help = 'Specify the number of replicates you have (required). Accepts a vector of numbers (e.g. 3,4,4,3).'),
	make_option(c('--conditions', '-c'),
		type = 'character',
		default = '',
		help = 'Specify the conditions of your samples in the order as they appear in the counts files, separated by commas (required).'),
	make_option(c('--samples','-s'),
		type = 'character',
		default = '',
		help = paste0('Specify the sample names in same order as they appear in the counts files and in the order you specified the ',
			'counts files, seperated by commas.')),
	make_option('--omit_genes',
		type = 'character',
		default = FALSE,
		help = 'Specify a txt file with gene names to be omitted from analysis (1 line per gene ID). Must match the gene IDs in the counts file.'),
	make_option('--test_method',
		type = 'character',
		default = 'QLFTest',
		help = 'Select the testing method for DE to be applied. Choose from "QLFTest" (default), "QLFTreat", or "glmLRT".'),
  make_option('--save_ggplot_obj',
    action = 'store_true',
    default = FALSE,
    help = 'Save the ggplot2 objects of the generated plots into an .rdata file.'),
	make_option('--config',
		type = 'character',
		default = '',
		help = 'Specify a custom config file with path.'),
	make_option(c('--verbose', '-V'),
		action = 'store_true',
		default = FALSE,
		help = 'Run in verbose mode (off by default).'),
	make_option('--debug',
		action = 'store_true',
		default = FALSE,
		help = 'Run in debug mode (off by default).')
)

# PARSING & PREPERATION
# set default config
script_path <- Sys.which('difExp.R')
script_path <- str_split_1(script_path, '"')
utp_path <- str_split_1(script_path, '/')
utp_path <- utp_path[1:(length(utp_path)-2)]
def_config <- str_c(c(utp_path, 'config/difExp.config'), collapse = '/')

# parse input arguments
difExp_arguments <- OptionParser(option_list=difExp_options)
difExp_parsed <- parse_args(difExp_arguments)

# break if required arguments are missing
if (difExp_parsed$input_counts == '') {
	stop('Input parameter missing. Exiting...')
}

if (difExp_parsed$output_dir == '') {
	stop('Output directory parameter missing. Exiting...')
}

if (difExp_parsed$replicates == '') {
	stop('Replicates parameter missing. Exiting...')
}

if (difExp_parsed$conditions == '') {
	stop('Conditions parameter missing. Exiting...')
}



# PROCESS INPUT ARGUMENTS
input_files <- str_split_1(difExp_parsed$input_counts, ',')
output_dir <- difExp_parsed$output_dir
output_name <- difExp_parsed$name
replicates <- as.numeric(str_split_1(difExp_parsed$replicates, ','))
#print(replicates) #tmpprint
#print(class(replicates)) #tmpprint
#print(str(replicates)) #tmpprint
conditions <- str_split_1(difExp_parsed$conditions, ',')
save_ggplot_obj <- difExp_parsed$save_ggplot_obj
config <- difExp_parsed$config
verbose <- difExp_parsed$verbose
debug <- difExp_parsed$debug
sample_names <- difExp_parsed$samples
omit_genes <- difExp_parsed$omit_genes
test_method <- difExp_parsed$test_method

# set base filename
# write data to this directory and name
if (endsWith(output_dir, '/') == FALSE) {
	output_dir <- paste(output_dir, '/', sep='')
}

output_base <- paste(output_dir, output_name, sep='')

# process input arguments
# set default config if no config is submitted
if (config == '') {
	config <- def_config
}

# read config file
source(config)

# check if outdir exists and create if not
if (dir.exists(output_dir) == FALSE) {
	dir.create(output_dir)
}

# check if ggplot objects are to be stored and create outdir if so
if (save_ggplot_obj == TRUE) {
  outdir_ggplot <- paste0(output_dir, 'ggplot/')
  dir.create(outdir_ggplot)
}

# check if input files exist and break if a submitted file does not exist
input_index <- length(input_files)
if (input_index == 1) {
	if (file.exists(input_files[[input_index]]) == FALSE) {
		stop(paste('Error: input counts file', input_files[[input_index]], 'does not exist'))
	}
} else {
	for (ii in 1:input_index) {
		if (file.exists(input_files[[ii]]) == FALSE) {
			stop(paste('Error: input counts file', input_files[[ii]], 'does not exist'))
		}
	}
}

# check if sample names have been submitted & create list if they were
if (sample_names != '') {
	sample_names <- str_split_1(sample_names, ',')
}

# if not verbose -> silence warnings
if (verbose == FALSE) {
	defaultW <- getOption("warn")
	options(warn = -1)
}

# close argument processing
if (verbose == TRUE) {
	print('Arguments succesfully processed')
}



# SET FUNCTIONS

# function to read the counts data from one file into a table
readCounts <- function(counts, sample_names, verbose, debug) {
	# set running mode
	if (verbose == TRUE) {
		print(paste('Started reading file',counts))
	}
	if (debug == TRUE) {
		browser()
	}

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

	# print the header of the input if the verbose tag is specified
	if (verbose == TRUE) {
		print(head(input_counts_tidy))
	}
	return(list(input_counts_tidy, input_counts_gene_length))
}	

# function to merge multiple counts tables into a single table. Assumes that all submitted tables have the same genes/regions.
mergeCounts <- function(counts, sample_names, verbose, debug) {
	# set running mode
	if (verbose == TRUE) {
		print('Start merging count files')
	}
	if (debug == TRUE) {
		browser()
	}

	# read the counts into tables and put these into a list
	input_counts_list <- list()
	for (i in 1:length(counts)) {
		input_counts_tmp <- readCounts(counts[i], '', verbose, debug)	# use the earlier defined function to read the data
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
	# print the header of the final table if the verbose tag is specified
	if (verbose == TRUE) {
		print(head(input_counts_merged))
	}
	return(list(input_counts_merged, input_counts_gene_length))
}

# function to normalize counts - re-organizes data, normalizes, returns tibble
normalizeCounts <- function(input_counts, input_counts_gene_length, conditions_fct, replicates_fct, verbose, debug, norm_library, norm_sample) {
	# set running mode
	if (verbose == TRUE) {
		print('Normalizing counts')
	}
	if (debug == TRUE) {
		browser()
	}

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

# if there are >2 conditions, conditions need to be paired appropriately so that e.g.
# A-B, A-C, A-D, B-C, B-D, C-D are compared, with the right number of replicates.
conditionsManyToPair <- function(conditions, replicates) {
	# make empty lists to pair conditions
	conditions_out <- vector('list', choose(length(conditions), 2)) # conditions re-organized into a list of lists
	cond_index_out <- vector('list', choose(length(conditions), 2)) # list of lists with the columns needed for the analysis
	replicates_out <- vector('list', choose(length(conditions), 2)) # list of vectors with replicate numbers
	cond_x_rep_out <- vector('list', choose(length(conditions), 2)) # list of conditions x number of replicates

	# organize conditions to pair
	h <- 1
	for (i in 1:(length(conditions) - 1)) {
		# get the first half of a pair
		pair_1 <- conditions[[i]]
		indx_1 <- c((1:replicates[[i]]) + sum(replicates[(1:i - 1)]) + 1)
		repl_1 <- c(1:replicates[[i]])
		cond_1 <- rep(conditions[[i]], times=replicates[[i]])
		for (j in (1 + i):length(conditions)) {
			# get the second half of a pair
			pair_2 <- conditions[[j]]
			indx_2 <- c((1:replicates[[j]]) + sum(replicates[(1:j - 1)]) + 1)
			repl_2 <- c(1:replicates[[j]])
			cond_2 <- rep(conditions[[j]], times=replicates[[j]])
			# pair the halves
			pairs <- c(pair_1, pair_2)
			indxs <- c(indx_1, indx_2)
			repls <- c(repl_1, repl_2)
			conds <- c(cond_1, cond_2)
			# append to list
			conditions_out[[h]] <- pairs
			cond_index_out[[h]] <- indxs
			replicates_out[[h]] <- repls
			cond_x_rep_out[[h]] <- conds
			h <- (h + 1)
		}
	}
	return(list(conditions_out, cond_index_out, replicates_out, cond_x_rep_out))
}

# set experimental design (i.e. tells the analyzing functions what the columns are)
makeDesign <- function(counts, conditions_fct, replicates_fct, conditions) {
	design_tibble <- as_tibble(data.frame(Sample=base::colnames(dplyr::select(counts, -Geneid)), conditions_fct, replicates_fct))
	design_matrix <- as.data.frame(model.matrix(~0+conditions_fct, data=design_tibble))

	colnames(design_matrix) <- sort(conditions)  #columns are placed in alphabetical order rather than the specified order

  design_matrix <- as.matrix(design_matrix[conditions])  # this re-orders the columns to the order in which they were specified

	if (verbose == TRUE) {
		print(design_matrix)
	}
	return(design_matrix)
}

# perform differential expression analysis
makeDGE <- function(counts, conditions_fct, replicates_fct, conditions, design_matrix, verbose, debug) {
	# set running mode
	if (verbose == TRUE) {
		print('Performing differential expression analysis')
	}
	if (debug == TRUE) {
		browser()
	}

	# transform & prepare data for differential expression analysis
	counts_DGE <- counts %>%
		dplyr::select(-Geneid) %>%
		DGEList(genes=counts$Geneid, group=conditions_fct) %>%
		estimateDisp(design_matrix)
	return(counts_DGE)
}

difExp <- function(counts_DGE, method, adjust_method, lfc, design_matrix, verbose, debug, mode = 'full', test_method = 'QLFTest') {
	# set running mode
	if (verbose == TRUE) {
		print('Performing differential expression analysis')
	}
	if (debug == TRUE) {
		browser()
	}

	# set contrast
	if (ncol(design_matrix) == 2) {
		contrast <- c(-1,1)
	} else {
		contrast <- NULL
	}

  print('')                   #tmpprint
  print('head(counts_DGE)')   #tmpprint
  print('')                   #tmpprint
  print(head(counts_DGE))     #tmpprint

	# apply model to data
	if (test_method == 'QLFTreat') {
		counts_DGELRT <- counts_DGE %>%
			glmQLFit(design=design_matrix, robust=robust, abundance.trend=abundance_trend, winsor.tail.p=winsor_p) %>%
			glmTreat(coef=ncol(.$design), contrast=contrast, lfc=lfc)
	} else if (test_method == 'QLFTest') {
		counts_DGELRT <- counts_DGE %>%
			glmQLFit(design=design_matrix, robust=robust, abundance.trend=abundance_trend, winsor.tail.p=winsor_p) %>%
			glmQLFTest(coef=ncol(.$design), contrast=contrast)
	} else if (test_method == 'glmLRT') {
		counts_DGELRT <- counts_DGE %>%
			glmFit(design=design_matrix) %>%
			glmLRT(coef=ncol(.$design), contrast=contrast)
	}

	# determine if full DGELRT table should be returned or the function should continue
	if (mode == 'half') {
		return(counts_DGELRT)
	} else if (mode == 'full') {

		# perform differential expression analysis & store results in a tibble
		if (test_method == 'QLFTreat') {
			counts_DE <- as_tibble(as.data.frame(decideTests(counts_DGELRT, adjust.method=adjust_method, p.value=0.05)))
		} else {
			counts_DE <- as_tibble(as.data.frame(decideTests(counts_DGELRT, adjust.method=adjust_method, p.value=0.05, lfc=lfc)))
		}

		counts_DE <- cbind(counts_DGELRT$genes, counts_DE)

		# re-name columns appropriately
		counts_name <- paste(c('-', colnames(design_matrix)[[1]], '_+', colnames(design_matrix)[[2]]), collapse='')
		base::colnames(counts_DE) <- c('Geneid',counts_name)

		# change numeric labels to character labels (1 to 'up', 0 to 'nonsig', -1 to 'down')
		counts_DE[[counts_name]] <- counts_DE %>% with(ifelse(
			.[[counts_name]] == 0, 'nonsig',
			ifelse(.[[counts_name]] == 1, 'up', 'down')))

		# create a copy of the DE results table with the P-values and Log2(FoldChange) in additional columns
		# in the case of >2 conditions, the DE tables are merged, but with the P-values and Log2(FC), this
		# doesn't work, so they are stored in a separate table.
		counts_DE_full <- counts_DE
		counts_DE_full[['log2FC']] <- counts_DGELRT$table[['logFC']]
		counts_DE_full[['PValue']] <- p.adjust(counts_DGELRT$table[['PValue']], 'BH')

		# return differential expression output
		return(list(counts_DE, counts_DE_full))
	}
}

# create volcano plot of DE data
ggplotVolcanoForDE <- function(data, lfc, name) {

	# prepare data
	de_colname <- colnames(data)[[2]]
	de_name <- paste(c(name, de_colname), collapse='_')
	my_colors <- c(up = 'blue', nonsig = 'black', down = 'red')

	## get gene IDs of genes that are differentially expressed to display on the plot
	data$de_label <- NA
	data$de_label[data[,2] != "nonsig"] <- data$Geneid[data[,2] != "nonsig"]
	data_grouped <- group_by_at(data, 2)
 
  ## prepare title
  comp <- str_replace(de_colname, '-', '')
  comp <- str_replace(comp, '\\+', '')
  comp <- str_split(comp, '_')[[1]]
  comp_title <- paste0(comp[[1]], ' versus ', comp[[2]])

	# create plot
	plot_volcano <- ggplot(
		data = data_grouped,
		aes(
			x = log2FC,
			y = -log10(PValue),
			color = data_grouped[[de_colname]],
			label = de_label)) +
		geom_point(size=3) +
		theme_classic() +
		scale_color_manual(values = my_colors) +
		labs(
			title = paste0('Differential expression in ', name, ', comparison ', comp_title),
			subtitle = paste0('N = ', nrow(data), ': ', 
				nrow(filter(data, data[[de_colname]] == 'up')), ' up, ', 
				nrow(filter(data, data[[de_colname]] == 'down')), ' down.'),
			color = 'Differential expression') +
		xlab('log2(FC)') +
		ylab('-log10(adjusted P-value)')

		if (nrow(filter(data, data[[de_colname]] != 'nonsig')) <= 10) {
			plot_volcano <- plot_volcano + geom_text_repel()
		}

	return(list(plot_volcano, de_name))

}

# create scatterplot for DE data
ggplotScatterForDE <- function(data, output_name, comparison) {

	# select the two columns
	colname_1 <- as.character(colnames(data)[2]) #1st column is always called 'Geneid'
	colname_2 <- as.character(colnames(data)[3])

	# create the output name
	name_out <- paste(output_name, '_', comparison, sep = '')

	# generate the plot
	scatterplot <- ggplot(data) +
		geom_pointdensity(aes(x = log2(data[[colname_1]]), y = log2(data[[colname_2]]))) +
		scale_color_gradient(name='Gene density') +
		theme_classic() +
		labs(
			x = paste0('Log2 Mean expression of ', colname_1),
			y = paste0('Log2 Mean expression of ', colname_2),
			title = paste0('Expression in ', output_name, ': ', colname_2, ' versus ', colname_1),
			subtitle = paste0('N = ', nrow(data))) +
		geom_abline(slope = 1) #this line (x = y) highlights any deviations from the diagonal

	return(list(scatterplot, name_out))

}



# DATA PREPERATION

if (verbose == TRUE) {
	print('Preparing data')
}

# create empty vectors to parse into
conditions_fct <- c()
replicates_fct <- c()

if (length(conditions) == 2) {

	for (r in 1:length(replicates)) {
		conditions_fct <- c(conditions_fct, rep(conditions[[r]], each=replicates[[r]]))
		replicates_fct <- c(replicates_fct, 1:replicates[[r]])
	}
 
} else {

	cond_replc_paired <- conditionsManyToPair(conditions, replicates)
	conditions_paired <- cond_replc_paired[[1]]     # vector with paired conditions (A-B, A-C, B-C)
	cond_index_paired <- cond_replc_paired[[2]]     # vector with appropriate columns to select for those pairings (123-456, 123-789, 456-789)
	replicates_vector <- cond_replc_paired[[3]]		# vector with replicates paired to each condition (abc-abc, abc-abc, abc-abc)
	conditions_vector <- cond_replc_paired[[4]]		# vector with paired conditions (AAA-BBB, AAA-CCC, BBB-CCC)

	for (r in 1:length(replicates)) {
		conditions_fct <- c(conditions_fct, rep(conditions[[r]], each=replicates[[r]]))
		replicates_fct <- c(replicates_fct, 1:replicates[[r]])
	}

}

conditions_fct <- as.factor(conditions_fct)
replicates_fct <- as.factor(replicates_fct)

#print(conditions_fct) #tmpprint
#print(replicates_fct) #tmpprint



# NORMALIZE INPUT AND STORE IN DATAFRAMES

# importing genes to be removed (if applicable)
if (omit_genes != FALSE) {
	rm_genes <- read_tsv(omit_genes, col_names = FALSE)
	colnames(rm_genes) <- 'Geneid'
}

# importing main data
if (length(input_files) > 1) {

	# first we merge the count tables
	input_data <- mergeCounts(input_files, sample_names, verbose, debug)
	input_table <- input_data[[1]]
	input_counts_gene_length <- input_data[[2]]

	# remove any genes if specified
	if (omit_genes != FALSE) {
		if (verbose == TRUE) {
			print(paste0('Removing ', nrow(rm_genes), ' genes from analysis.'))
		}
		input_table <- filter(input_table, !(Geneid %in% rm_genes$Geneid))
	}

	# then we normalize everything
	counts_normalized_df <- normalizeCounts(input_table, input_counts_gene_length, 
		conditions_fct, replicates_fct, verbose, debug,
		norm_library, norm_sample)

} else {

	# no tables need to be merged so we normalize immediately
	input_data <- readCounts(input_files[[1]], sample_names, verbose, debug)
	input_table <- input_data[[1]]
	input_counts_gene_length <- input_data[[2]]

	# remove any genes if specified
	if (omit_genes != FALSE) {
		if (verbose == TRUE) {
			print(paste0('Removing ', nrow(rm_genes), ' genes from analysis.'))
		}
		input_table <- filter(input_table, !(Geneid %in% rm_genes$Geneid))
	}

	# then we normalize everything
	counts_normalized_df <- normalizeCounts(input_table, input_counts_gene_length, 
		conditions_fct, replicates_fct, verbose, debug,
		norm_library, norm_sample)

}

counts_normalized_df_untidy <- log(dplyr::select(counts_normalized_df, -1))
rownames(counts_normalized_df_untidy) <- counts_normalized_df[['Geneid']]
#print(head(counts_normalized_df_untidy)) #tmpprint



# create table with the row average per condition
data_avr <- dplyr::select(counts_normalized_df, 'Geneid')
for (ist in 1:length(conditions)) {
	cond <- conditions[[ist]]
	cols <- 
	data_avr_in <- dplyr::select(counts_normalized_df, (1:replicates[[ist]]) + sum(replicates[(1:ist - 1)]) + 1)
	data_avr[[cond]] <- rowMeans(data_avr_in)
}



# PERFORM DIFFERENTIAL EXPRESSION ANALYSIS
## make design matrix
data_design <- makeDesign(counts_normalized_df, conditions_fct, replicates_fct, conditions)

## apply model to the data
data_DGE <- makeDGE(counts_normalized_df, conditions_fct, replicates_fct, conditions, data_design, verbose, debug)

## calculate differential expression by condition
if (length(conditions) == 2) {

	data_DE_list <- difExp(data_DGE,  method, adjust_method, lfc, data_design, verbose, debug, mode = 'full', test_method = test_method)
	data_DE <- data_DE_list[[1]]
	data_DE_full <- data_DE_list[[2]]

} else {

	data_DE <- dplyr::select(counts_normalized_df, 'Geneid')
	data_DE_full <- vector('list', (length(conditions_paired)))

	for (z in 1:length(conditions_paired)) {

		# select subset of data
		data_DGE_subset <- data_DGE[,((cond_index_paired[[z]]) - 1)]	# -1 because here there is no first column with GeneIDs
		data_counts_subset <- dplyr::select(counts_normalized_df, 1, cond_index_paired[[z]])

		# create design
		data_design_tmp <- makeDesign(data_counts_subset, conditions_vector[[z]], replicates_vector[[z]], conditions_paired[[z]])

		# perform analysis & format/merge results
		data_DE_list_tmp <- difExp(data_DGE_subset, method, adjust_method, lfc, data_design_tmp, verbose, debug, mode = 'full', test_method = test_method)
		data_DE_tmp <- data_DE_list_tmp[[1]]
		data_DE_full[z] <- data_DE_list_tmp[2]
		data_DE <- full_join(data_DE, data_DE_tmp, by='Geneid')

	}

}

comparisons <- colnames(dplyr::select(data_DE, -1))



# WRITE DATA TO OUTPUT FILES
output_table_all_name <- paste0(output_base, '_all.csv')
output_table_sig_name <- paste0(output_base, '_sig.csv')
output_table_avr_name <- paste0(output_base, '_avr.csv')
output_table_nrm_name <- paste0(output_base, '_nrm.csv')

if (verbose == TRUE) {
	print(paste('Writing full DE table to', output_table_all_name))
	print(paste('Writing DE genes to', output_table_sig_name))
}

# select significant genes
data_DE_sig <- filter(data_DE, if_any(-Geneid, ~ . != 'nonsig'))

#print('data_DE_sig')		#tmpprint
#print(head(data_DE_sig))	#tmpprint

write_csv(data_DE, output_table_all_name)
write_csv(data_DE_sig, output_table_sig_name)
write_csv(data_avr, output_table_avr_name)
write_csv(counts_normalized_df, output_table_nrm_name)

if (length(conditions) == 2) {
	output_table_full_name <- paste(output_base, '_full.csv', sep='')
	if (verbose == TRUE) {
		print(paste('Writing full DE table with Log2(foldchange) and P-values'))
	}
	data_DE_full_posthoc <- data_DE_full
	data_DE_full_posthoc$PValue <- p.adjust(data_DE_full[['PValue']], 'BH')
	write_csv(data_DE_full_posthoc, output_table_full_name)
} else {
	for (oti in 1:length(data_DE_full)) {
		data_DE_tmp_out <- data_DE_full[[oti]]
		conds_table_out <- conditions_paired[[oti]]
		output_table_full_name <- paste(output_base, '_-', conds_table_out[[1]], '_+', conds_table_out[[2]], '_full.csv', sep='')
		if (verbose == TRUE) {
			print(paste('Writing full DE table with Log2(foldchange) and P-values'))
		}
		data_DE_full_posthoc <- data_DE_tmp_out
		data_DE_full_posthoc$PValue <- p.adjust(data_DE_tmp_out[['PValue']], 'BH')
		write_csv(data_DE_full_posthoc, output_table_full_name)
	}
}



# DATA VISUALIZATION
## volcano plots
if (length(conditions) == 2) {
	plot_volcano_out <- ggplotVolcanoForDE(data_DE_full, lfc, output_name)
	plot_volcano <- plot_volcano_out[[1]]
	name_volcano <- plot_volcano_out[[2]]
	file_volcano <- paste0(output_dir, name_volcano, '_volcano.png')
	ggsave(file_volcano, plot = plot_volcano)
  if (save_ggplot_obj == TRUE) {
    ggfile_volcano <- paste0(outdir_ggplot, name_volcano, '_volcanoplot.rds')
    saveRDS(plot_volcano, file = ggfile_volcano)
  }
} else {
	for (ci in 1:length(data_DE_full)) {	# changed from length(conditions_vector)
		plot_volcano_out <- ggplotVolcanoForDE(data_DE_full[[ci]], lfc, output_name)
		plot_volcano <- plot_volcano_out[[1]]
		name_volcano <- plot_volcano_out[[2]]
		file_volcano <- paste(output_dir, name_volcano, '_volcano.png', sep='')
		ggsave(file_volcano, plot = plot_volcano)
    if (save_ggplot_obj == TRUE) {
      ggfile_volcano <- paste0(outdir_ggplot, name_volcano, '_volcanoplot.rds')
      saveRDS(plot_volcano, file = ggfile_volcano)
    }
	}
}

## correlation plots
data_correlation <- cor(dplyr::select(counts_normalized_df, -Geneid))
name_correlation <- paste(output_base, '_corrplot.png', sep='')
number_samples <- replicates * length(conditions)
cplt_size <- 200 + 60 * number_samples # dynamic plot size

png(file = name_correlation, width = cplt_size, height = cplt_size, units = 'px')
corrplot.mixed(
	data_correlation,
	upper = 'circle',
	lower = 'number')
dev.off()



## scatterplots
if (length(conditions) == 2) {
	plot_scatter_out <- ggplotScatterForDE(data_avr, output_name, comparisons)
	plot_scatter <- plot_scatter_out[[1]]
	name_scatter <- plot_scatter_out[[2]]
	file_scatter <- paste(output_dir, name_scatter, '_scatterplot.png', sep='')
	ggsave(file_scatter, plot = plot_scatter)
  if (save_ggplot_obj == TRUE) {
    ggfile_scatter <- paste0(outdir_ggplot, name_scatter, '_scatterplot.rds')
    saveRDS(plot_scatter, file = ggfile_scatter)
  }
} else {
	data_avr_colnames <- colnames(data_avr)
	for (si in 1:length(conditions_paired)) {
		scatter_conds <- conditions_paired[[si]]
		scatter_cond_1 <- which(data_avr_colnames == scatter_conds[[1]])
		scatter_cond_2 <- which(data_avr_colnames == scatter_conds[[2]])
		scatter_data_tmp <- select_(data_avr, 1, scatter_cond_2, scatter_cond_1) # necessary as they are also inverted in order
		scatter_comp <- as.character(comparisons[[si]])
#		print('scatter_comp') #tmpprint
#		print(scatter_comp) #tmpprint
		plot_scatter_out <- ggplotScatterForDE(scatter_data_tmp, output_name, scatter_comp)	 # in the case with 2 conditions
		plot_scatter <- plot_scatter_out[[1]]
		name_scatter <- plot_scatter_out[[2]]
		file_scatter <- paste(output_dir, name_scatter, '_scatterplot.png', sep='')
		ggsave(file_scatter, plot = plot_scatter)
    if (save_ggplot_obj == TRUE) {
      ggfile_scatter <- paste0(outdir_ggplot, name_scatter, '_scatterplot.rds')
      saveRDS(plot_scatter, file = ggfile_scatter)
    }
	}
}

## plotMDS
name_MDS <- paste(output_base, '_MDS.png', sep='')
png(file = name_MDS, width = 800, height = 800, units = 'px')
MDS_colours <- c('red', 'blue', 'green', 'magenta', 'orange', 'pink', 'yellow', 'blue4')
plot_MDS <- plotMDS(counts_normalized_df_untidy, col = MDS_colours[conditions_fct])
dev.off()



# EPILOGUE
# turn on warnings again in environment
if (verbose == FALSE) {
	options(warn = defaultW)
}
