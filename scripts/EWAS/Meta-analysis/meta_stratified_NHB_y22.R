# title: meta_NHB_stratified.R
# author: Lisa Schneper
# date: 05-04-2026

.libPaths("/path to R package library")

library(metafor)
library(optparse)

extract_remlsummary<-function(res_list){
  do.call(rbind, lapply(res_list,function(x){
    if(is.null(x)|| all(is.na(x))) return(NULL)
    
    data.frame(
      beta_NHB = x$b[1],
      beta_race = x$b[2],
      
      se_NHB = x$se[1],
      se_race = x$se[2],
      
      p_NHB = x$pval[1],
      p_race = x$pval[2]
      
      
    )
  }))
}
cleandf_fun<-function(df){
  df<-df[!duplicated(df),]
  df<-df[!df$probe%in% "probe",]
  df$beta<-as.numeric(df$Coefficient)
  df$se<-as.numeric(df$StdError)
  return(df)
}
setwd("path to EWAS output")
args<-commandArgs(trailingOnly=TRUE)
# change for other ages
agenum=22
#agenum=as.integer(args[1])
age<-paste0("y",agenum)

file1<-read.delim("path to NHB_stratified_ewas_model1_results.txt")
file1$race<-"NHB"
file1$array="epic"

file2<-read.delim("path to notNHB_stratified_ewas_model1_results.txt")
file2$race<-"OTHER"
file2$array="epic"


file1<-cleandf_fun(file1)
file2<-cleandf_fun(file2)

file1_2_cpg<-file1$probe[file1$probe %in% file2$probe]
file1sub<-file1[file1$probe %in% file1_2_cpg,]
file2sub<-file2[file2$probe %in% file1_2_cpg,]

both_arrays<-rbind(file1sub,file2sub)
both_arrays$sex<-factor(both_arrays$race,level=c("NHB","OTHER"))

both_arrays<-both_arrays[order(both_arrays$probe,both_arrays$race),]

cpg_list<-split(both_arrays,both_arrays$probe)
#end_cpg[end_cpg > length(cpg_list)] <- length(cpg_list)
#cpg_list<-cpg_list[start_cpg:end_cpg]
# Run meta analysis
metalist<-list()
cat("Began processing chunk with", length(cpg_list), "CpGs at ",Sys.time(),"\n")
metalist <- lapply(cpg_list,function(x){
  rma(
    yi=x$beta,
    sei=x$se,
    mods = ~x$race,
    method="FE"
  )
})
cat("Finished processing ", length(cpg_list), "CpGs at ",Sys.time(),"\n")

reml_meta<-extract_remlsummary(metalist)
reml_meta$fdr<-p.adjust(reml_meta$p_NHB,method="fdr")
reml_meta$CpG<-names(cpg_list)
write.table(reml_meta,file=paste0("NHB_OTHER_stratified_",age,"_meta_FE_output.txt"),sep="\t",row.names=F,quote=F)
