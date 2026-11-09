

# Author: Eve Schoeffler (adapted from scripts from Yinan Zheng)

# This script is intended to split a large RGChannelSet into plate-sized batches.

# Load libraries
.libPaths("path to R package library")
library(minfi)

# Set working directory
setwd("path to working directory")

# Read in RDS file
big_rgset <- readRDS("path to extended RGset")

#### Create a list of object names to hold subsets of RDS file
# Initialize the starting and ending samples
start <- 1
end <- 96

# Define the last sample
last.sample <- big_rgset@colData@nrows

subsets <- vector(mode = "list", length = 28)

# Create a folder to hold the batches
dir.create("Plate-sized_batches")

for (elem in 1:length(subsets)) {
  # Create the subset
  print(paste0("Saving subset ", elem))
  saveRDS(big_rgset[, start:end], file = paste0("Plate-sized_batches/RGSet", elem))
  
  # Increment start and end
  start <- start + 96
  end <- ifelse(end + 96 > last.sample, last.sample, end + 96)
}

# Create folders for the next step (QC splitter)
dir.create("QCsplitter")
dir.create("QCsplitter/bisul/")
dir.create("QCsplitter/qcmat/")
dir.create("QCsplitter/mu/")
dir.create("QCsplitter/beta/")
dir.create("QCsplitter/detP/")
dir.create("QCsplitter/nbead/")
dir.create("QCsplitter/badsample/")

# Quit out of R
quit(save = "no")
