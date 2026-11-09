# Set libPaths
.libPaths("path to R package library")

# Load libraries
library(illuminaio)
library(boot)
library(magrittr)
library(doParallel)
# uncomment for 450K data
#library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(ENmix)

setwd("path to /Preprocess/ directory")

# Combine the methyl set batches into one methylset
msets <- list.files("MethylSets")
msets_info <- file.info(msets)
sorted_msets <- msets[order(msets_info$mtime)]
mset_all <- NULL
for (methyl_set in sorted_msets) {
  mset_new <- readRDS(methyl_set)
  mset_all <- combineArrays(mset_all, mset_new)
}


# print(paste("Started processing at: ",Sys.time(),sep=""))

# Quantile normalize
qnorm <- norm.quantile(mset_all)

# Save
saveRDS(qnorm, file="MethylSet_ENmix_qnorm.combinedEPIC.saliva.n2690.RDS")
rm(pre)
gc()

# Probe type bias correction
beta<-rcp(qnorm)
rm(qnorm)
gc()

# Save
saveRDS(beta,file="BetaMatrix.combinedEPIC.saliva.n2690.RDS")

print(paste("Finished processing at ",Sys.time(),sep=""))

quit(save="no") 
