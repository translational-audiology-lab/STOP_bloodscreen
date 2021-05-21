# -----------------------------------------------------------------------------#
# Run this in R as below to get all intermediate data files.
#    source("master_data_generation.R")
# All other R scripts that generate data files for data analysis are executed. 
# -----------------------------------------------------------------------------#
# initiated on 2021-05-07
# authors : Mun-Gwan Hong
# -----------------------------------------------------------------------------#
rm(list = ls())

# Read the table and create RData file of clinical info
source("gen-s1-clinical.v01.R")

# Parse Olink proteomic data
source("gen-s1-olink_proteomic.v01.R")

# Preprocessing of Olink Proteomic data
#   - details in `report-QC-Olink_proteomic_data.Rmd`
source("gen-s1-olink_proteomic.v02.R")

# Limit the samples to those with proteomic data
source("gen-s1-clinical.v02.R")
