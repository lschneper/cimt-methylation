.libPaths("path to R package library")
library(bacon)
library(meta)
library(dplyr)
library(tidyr)
library(ggplot2)

args<-commandArgs(trailingOnly=TRUE)
age<-as.integer(args[1])
model_num<-as.integer(args[2])

INPUT_FOLDER<-paste0("path to input file directory",age,"/sensitivity/")

EPIC_INPUT_FILE<-paste0(INPUT_FOLDER,"y",age,"_epic_sensitivity_model",model_num,".ewasoutput.txt")
HM450K_INPUT_FILE<-paste0(INPUT_FOLDER,"y",age,"_450K_sensitivity_model",model_num,".ewasoutput.txt")

OUTPUT_FOLDER="path to output directory"

dir.create(OUTPUT_FOLDER)
OUTPUT_FILESTEM=paste0("y",age,"_meta_sensitivity_model",model_num)
epicfile=EPIC_INPUT_FILE
hm450file=HM450K_INPUT_FILE

epic <- read.table(epicfile, sep = "\t", header = T)
epic<-epic[epic$Coefficient !="Coefficient",]
epic$Coefficient<-as.numeric(epic$Coefficient)
epic$StdError<-as.numeric(epic$StdError)

hm450 <- read.table(hm450file, sep = "\t", header = T)
hm450<-hm450[hm450$Coefficient !="Coefficient",]
hm450$Coefficient<-as.numeric(hm450$Coefficient)
hm450$StdError<-as.numeric(hm450$StdError)


# Find overlap of CpGs among all results
epic<-epic[epic$probe %in% hm450$probe,]
hm450 <- hm450[hm450$probe %in% epic$probe, ]
epico<-epic[match(hm450$probe,epic$probe),]
epico$study<-"epic"
hm450$study<-"450k"

both_arrays<-rbind(epico,hm450)
both_arrays<-both_arrays[order(both_arrays$probe,both_arrays$study),]

cpg_list<-split(both_arrays,both_arrays$probe)

# Run meta analysis
metalist<-list()
cat("Began processing chunk with", length(cpg_list), "CpGs at ",Sys.time(),"\n")
metalist <- lapply(cpg_list,function(x){
  metagen(
    TE=x$Coefficient,
    seTE=x$StdError,
    studlab = x$study,
    comb.fixed=TRUE,comb.random=FALSE)
})


names(metalist)<-names(cpg_list)

saveRDS(metalist,file=paste0(OUTPUT_FOLDER,OUTPUT_FILESTEM,"_metalist_output.RDS"))
gc()
cat("finished processing chunk and began extrating stats with", length(cpg_list), "CpGs at ",Sys.time(),"\n")
meta_stats<-lapply(metalist,function(m){
  w.fixed<-round(100*m$w.fixed/sum(m$w.fixed),digits=2)
  w.random<-round(100*m$w.random/sum(m$w.random),digits=2)
  names(w.fixed)<-m$studlab
  names(w.random)<-m$studlab
  data.frame(
    TE_fixed=try(m$TE.fixed),
    seTE_fixed=try(m$seTE.fixed),
    p_fixed=try(m$pval.fixed),
    TE_random=try(m$TE.random),
    seTE_random=try(m$seTE.random),
    p_random=try(m$pval.random),
    I2=try(m$I2),
    Tau2=try(m$tau2),
    CochranQ=try(m$Q),
    w450k_fixed=try(w.fixed["450k"]),
    wEPIC_fixed=try(w.fixed["epic"]),
    w450k_random=try(w.random["450k"]),
    wEPIC_random=try(w.random["epic"])
  )
})


metalistdf<-do.call(rbind,meta_stats)
# remove random effect
metalistdf<-metalistdf[,grep("random",names(metalistdf),invert=T)]
metalistdf$fdr_fixed<-p.adjust(metalistdf$p_fixed,method="fdr")
metalistdf$bonferroni_fixed<-p.adjust(metalistdf$p_fixed,method="bonferroni")

gc()

cat("finished extrating stats with", length(cpg_list), "CpGs at ",Sys.time(),"\n")
saveRDS(metalistdf,file=paste0(OUTPUT_FOLDER,OUTPUT_FILESTEM,"_metalistdf_output.RDS"))

  
