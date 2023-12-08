#!/usr/bin/env Rscript

# Script to install any necessary packages for the difExp.R script
# Made by Alejandro Baars as part of his General Research Profile
# Utrecht University 2023

print('Checking and installing packages')

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

print('Checking and installing complete')
