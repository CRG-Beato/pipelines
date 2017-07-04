#==================================================================================================
# Created on: 2017-07-07
# Usage: ./diffex.R
# Author: Javier Quilez (GitHub: jaquol)
# Goal: Performs differential expression analysis between 2 conditions using sleuth
#==================================================================================================



#==================================================================================================
# Configuration variables and paths
#==================================================================================================

# installation --only needs to be done once
#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
#install.packages("devtools")
#devtools::install_github("pachterlab/sleuth")
#install.packages("RSQLite")
#biocLite("biomaRt")
#install.packages("DESeq2")
#install.packages("pheatmap")

# load R libraries
library("sleuth")
library("RSQLite")
library("DESeq2")
library("pheatmap")

# variables
samples <- c('rz_005_01_01_rnaseq', 'rz_005_02_01_rnaseq', 'rz_005_03_01_rnaseq',
           'rz_006_01_01_rnaseq', 'rz_006_02_01_rnaseq', 'rz_006_03_01_rnaseq',
           'rz_007_01_01_rnaseq', 'rz_007_02_01_rnaseq', 'rz_007_03_01_rnaseq',
           'rz_008_01_01_rnaseq', 'rz_008_02_01_rnaseq', 'rz_008_03_01_rnaseq')
data_type <- "rnaseq"
project <- "rzaurin"
analysis <- "2016-12-12_differential_expression_analysis"
species = "homo_sapiens"
assembly_version = "hg38_mmtv"
sequencing_type <- 'paired_end'

# paths
HOME <- "/users/GR/mb/jquilez"
SAMPLES <- paste0(HOME, "/data/", data_type, "/samples")
ANALYSIS <- paste0(HOME, "/projects/", project, "/analysis/", analysis)
ifile_metadata <- paste0(HOME, "/data/beato_lab_metadata.db")