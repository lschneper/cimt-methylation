# Author: Eve Schoeffler (adapted from scripts from Yinan Zheng)
# 
# This script is intended to run the function QCsplitter on a single RGSet.

# Load libraries
.libPaths("path to R package library")
library(minfi)

# Set working directory
setwd("path to working directory")

# Take in argument
args <- commandArgs(TRUE)
batch <- as.numeric(args[1])

# Make a list of all the batches
batches <- list.files()

# Read in the batch that was given as an argument
print(paste0("Reading in RGSet for batch ", batch))
rgset <- readRDS(batches[batch])

# Load QCsplitter function
source("path to splitter_functions.R")

# Apply QCsplitter
qc <- QCsplitter(rgset,detPthre=0.000001,detPtype="negative",nbthre=3,samplethre=0.05,CpGthre=0.05,
                            bisulthre=NULL,outlier=TRUE, chunk_num=10)

##### Save the outputs

# Set output folders
out_bisul_folder <- "../QCsplitter/bisul/"
out_qcmat_folder <- "../QCsplitter/qcmat/"
out_mu_folder <- "../QCsplitter/mu/"
out_beta_folder <- "../QCsplitter/beta/"
out_detP_folder <- "../QCsplitter/detP/"
out_nbead_folder <- "../QCsplitter/nbead/"
out_badsample_folder <- "../QCsplitter/badsample/"

# Output bisul chunk
saveRDS(qc[["bisul"]], paste0(out_bisul_folder, "/allEPIC", batch, "_bisul_chunk.rds"))
message("Exported ", batch, " batch bisul.")

# Output a list of qc chunks split across CpGs
for(i in 1:length(qc[["qcmatlist"]])){
    saveRDS(qc[["qcmatlist"]][[i]], paste0(out_qcmat_folder, "/allEPIC_", batch, "_qcmat_CpG_chunk", i, ".rds"))
  }

message("Exported ", batch, " batch ", length(qc[["qcmatlist"]]), " chunks of QC matrix, split across CpGs.")

# Output total intensity value chunk
saveRDS(qc[["mu"]], paste0(out_mu_folder, "/allEPIC_", batch, "_total_intensity_chunk.rds"))
message("Exported ", batch, " batch total intensity values.")

# Output beta value chunk
saveRDS(qc[["beta"]], paste0(out_beta_folder, "/allEPIC_", batch, "_beta_chunk.rds"))
message("Exported ", batch, " batch beta values.")

# Output detection P value chunk
saveRDS(qc[["detP"]], paste0(out_detP_folder, "/allEPIC_", batch, "_detP_chunk.rds"))
message("Exported ", batch, " batch detection P values.")

# Output nbead chunk
saveRDS(qc[["nbead"]], paste0(out_nbead_folder, "/allEPIC_", batch, "_nbead_chunk.rds"))
message("Exported ", batch, " batch nbead values.")

# Output badsample chunk
saveRDS(qc[["badsample"]], paste0(out_badsample_folder, "/allEPIC_", batch, "_badsample_chunk.rds"))
message("Exported ", batch, " batch badsample.")

message("DONE!")
message(Sys.time())
writeLines(capture.output(sessionInfo()),"sessionInfo.txt")
