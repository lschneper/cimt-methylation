# Author: Eve Schoeffler

######### Note ###############################
# When loading idats from the CVH dataset, we found that there were differing numbers of probes between plates. To load all of the data in the same way, we used force = TRUE and randomized the order of loading. This pared down the larger number of probes to the smaller.
##############################################

# Set libPaths
.libPaths("path to R package library")

# Load libraries
library(boot)
library(minfi)
library(magrittr)
library(ENmix)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Set working directory path
setwd("path to output file directory (for RDS_files)/")


#Accept input start and end arguments
args<-commandArgs(TRUE)
start<-as.numeric(as.character(args[1]))
end<-as.numeric(as.character(args[2]))

# Read in idat sample sheet - need to modify path accordingly
idats<-readRDS("path to idat_ss.RDS")
# note this file has two columns:
# e.g.
# filename  basename
#205841310180_R01C01_Grn.idat  idatpath/205841310180_R01C01

tmp<-idats$basename[start:end]
rg<-read.metharray(tmp,extended=T, verbose = T, force = T)

# Save
saveRDS(rg,file=paste0("idat_",start,"_",end,"rgsetext_ES_beadpool.RDS"))

# Quit
quit(save="no")
