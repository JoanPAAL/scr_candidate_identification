#!/usr/bin/env Rscript

# This script calculates FPKMs 

# It should be placed on drac: 
# $ scp RIBOSEQ/R_SCRIPTS/calculate_fpkm.R jpallares@161.116.70.109:SCRIPTS/R_SCRIPTS/

# Usage: calculate_fpkm.R path_to_working_directory featurecounts.out output_name
# Example: time SCRIPTS/R_SCRIPTS/calculate_fpkm.R DATA/INSECT_RIBOSEQ_DATA/CULEX_TARSALIS/ featurecounts_merged.out fpkm.tsv


# Libraries
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("DGEobj.utils", quietly = TRUE))
    install.packages("DGEobj.utils")

BiocManager::install(version = "3.16", quietly=TRUE)
BiocManager::install("edgeR", quietly=TRUE)

library(DGEobj.utils)


# Execution
args = commandArgs(trailingOnly=TRUE)

setwd(args[1]) 
df <- read.csv(args[2], skip=1, sep = "\t")

#colnames(df)
mtx <- data.matrix(df[7:10]) # Counts as a Numeric Matrix of the Sample columns
length <- c(df$Length) # As a vector
Geneid <- c(df$Geneid)

# Relate Counts and Lenght:
fpkm_data <- convertCounts(mtx, unit = "fpkm", geneLength = length, 
                          log = FALSE, normalize = "tmm")

df_out <-  as.data.frame(cbind(Geneid, fpkm_data))
#df_out

write.table(df_out, file=args[3], quote=FALSE, sep='\t', row.names=FALSE)
