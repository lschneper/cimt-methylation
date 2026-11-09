

# Author: Eve Schoeffler (adapted from scripts by Lisa Schneper and Yinan Zheng)

# If you have a large RGset and are running out of memory trying to run preprocessENmix, run it in batches.
# Make batches of > 800 samples.
# This script takes arguments so that it can run on different "batches" within a folder and should be run in a for-loop according to the number of batches. There are different sections which should be (un)commented such as QCinfo, preprocessENmix, quantile normalization...according to your needs. At the end, you will have preprocessed beta matrices (in batches -- they will need to be combined in order to represent the beta matrix for your entire sample set).

# Set lib paths and load libraries
.libPaths("path to R package library")

library(boot)
library(minfi)
library(magrittr)
library(ENmix)
library(IlluminaHumanMethylationEPICmanifest)
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

# Set working directory
setwd("path to working directory (contains qc files in Final_QC)")

# Take in argument
args <- commandArgs(TRUE)
batch <- as.numeric(args[1])

rds_files <- list.files("Preprocess/RGSets")
# Read in RDS file
print(paste("Reading RGSet", batch))
rgset <- readRDS(paste0("Preprocess/RGSets/",rds_files[batch]))

# Read in QC file
print(paste("Reading in QC info..."))
qc_info <- readRDS("Final_QC/qc_info.EPIC.saliva.RDS")

# Preprocess
pre <- preprocessENmix(rgset, QCinfo = qc_info, nCores = 1)
print(paste("Saving methyl set", batch))
saveRDS(
  pre,
  paste0(
    "Preprocess/MethylSets/MethylSet_ENmix_preprocessed.blood.",
    batch,
    ".RDS"
  )
)
print("Done!")
