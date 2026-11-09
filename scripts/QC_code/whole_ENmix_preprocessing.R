library(boot)
library(minfi)
library(magrittr)
library(ENmix)
library(IlluminaHumanMethylationEPICmanifest)
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

.libPaths("path to R package library")
# Set working directory
setwd("path to working directory")

# Read in RDS file
rgset <- readRDS("path to extended RGChannelSet")

# Run QC
qc <- QCinfo(rgSet = rgset,distplot = F)

# Save QC info
saveRDS(qc, "qc_info.blood.RDS")
print("Done saving QC info object!")

qc_info <- readRDS("qc_info.blood.RDS")

# Preprocess and quantile normalize
pre <- preprocessENmix(rgset, QCinfo = qc_info, nCores = 1)
saveRDS(pre,"MethylSet_ENmix_prepreprocessed.CVH.blood.RDS")
print("Done saving first stage methyl set!")

# Read in the methyl set
pre <- readRDS("MethylSet_ENmix_prepreprocessed.CVH.blood.RDS")

# Quantile normalize
qnorm <- norm.quantile(pre)

# Save
saveRDS(qnorm, "MethylSet_quantile_normalized.CVH.blood.RDS")

# Probe type bias correction
beta <- rcp(qnorm)
rm(qnorm)
gc()

# Save
saveRDS(beta, "Beta_qnorm.CVH.blood.RDS")

# qcfilter
qcfilter <- qcfilter(beta)
saveRDS(qcfilter, "Beta.qcfilter.CVH.blood.RDS")
