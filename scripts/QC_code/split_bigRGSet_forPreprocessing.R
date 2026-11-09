# Author: Eve Schoeffler (adapted from scripts from Yinan Zheng)

# Split a large RGSet before preprocessing to run on batches (avoid memory issues)

.libPaths("path to R package library")
library(minfi)
setwd("path to working directory")

# Make a preprocessing directory
dir.create("Preprocess")
dir.create("Preprocess/RGSets")
dir.create("Preprocess/MethylSets")

# Read RGSet
rgset <- readRDS("path to RGChannelSet file")

# Split and save RGSet
last_sample <- rgset@colData@nrows
start <- 1
end <- 800
subsets <- 4 # Adjust as needed!

for (elem in 1:subsets) {
  # Create the subset
  print(paste0("Saving subset ", elem, " from index ", start, " to index ", end, "..."))
  saveRDS(rgset[, start:end], file = paste0("Preprocess/RGSets/RGChannelSet.saliva.", elem, ".RDS"))
  print(paste0("Done saving subset ", elem, "!"))
  
  # Increment start and end
  start <- start + 800
  end <- ifelse(end + 800 > last_sample, last_sample, end + 800)
}
