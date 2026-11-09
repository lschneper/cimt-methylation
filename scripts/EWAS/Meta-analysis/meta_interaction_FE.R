# title: meta_interaction.R
# author: Lisa Schneper
# date: 05-04-2026

.libPaths("path to R package library")

library(metafor)
library(optparse)
readlist_fun<-function(filelist){
  all_lists<-lapply(filelist,readRDS)
  combined<-do.call(c,all_lists)
  return(combined)
}

extract_summary<-function(res_list){
  do.call(rbind, lapply(res_list,function(x){
    if(is.null(x)|| all(is.na(x))) return(NULL)
    
    data.frame(
      CpG = x$CpG,
      beta_main = x$beta_main,
      se_main=x$se_main,
      p_main=x$p_main,
      p_int_global=x$p_int_global,
      beta_int = x$beta_int,
      se_int=x$se_int,
      p_int=x$p_int,
      n=x$n
      
    )
  }))
}

extract_remlsummary<-function(res_list){
  do.call(rbind, lapply(res_list,function(x){
    if(is.null(x)|| all(is.na(x))) return(NULL)
    
    data.frame(
      beta_meta = x$b,
      se_meta=x$se,
      p_meta=x$pval,
      z_meta=x$zval,
      ci_lb = x$ci.lb,
      ci_ub = x$ci.ub,
      Q=x$QE,
      Q_p=x$QEp,
      I2=x$I2,
      tau2=x$tau2,
      k=x$k
    )
  }))
}
setwd("path to file directory")
args<-commandArgs(trailingOnly=TRUE)

modelnum=as.integer(args[1])
model=paste0("model",modelnum)
agenum=as.integer(args[2])
age<-paste0("y",agenum)

files<-list.files(pattern=model)

epicfilelist<-files[grep(paste0(age,"_epic"),files)]
hm450Kfilelist<-files[grep(paste0(age,"_450K"),files)]
epicfiles<-readlist_fun(epicfilelist)
hm450Kfiles<-readlist_fun(hm450Kfilelist)
epicdf<-extract_summary(epicfiles)
hm450Kdf<-extract_summary(hm450Kfiles)
epicdf_meta<-epicdf[epicdf$CpG %in% hm450Kdf$CpG,]
hm450Kdf_meta<-hm450Kdf[hm450Kdf$CpG %in% epicdf$CpG,]
hm450Kdf_metao<-hm450Kdf_meta[match(epicdf_meta$CpG,hm450Kdf_meta$CpG),]
epico<-epicdf_meta
hm450<-hm450Kdf_metao
##### 
epico$study<-"epic"
hm450$study<-"450k"

both_arrays<-rbind(epico,hm450)
both_arrays<-both_arrays[order(both_arrays$CpG,both_arrays$study),]
both_arrays$probe<-both_arrays$CpG
write.table(both_arrays,file=paste0(model,"_",age,"_interaction_results.txt"),sep="\t",quote=F,row.names=F)
cpg_list<-split(both_arrays,both_arrays$probe)
#end_cpg[end_cpg > length(cpg_list)] <- length(cpg_list)
#cpg_list<-cpg_list[start_cpg:end_cpg]
# Run meta analysis
metalist<-list()
cat("Began processing chunk with", length(cpg_list), "CpGs at ",Sys.time(),"\n")
metalist <- lapply(cpg_list,function(x){
  rma(
    yi=x$beta_int,
    sei=x$se_int,
    method="FE"
    )
})
cat("Finished processing ", length(cpg_list), "CpGs at ",Sys.time(),"\n")

reml_meta<-extract_remlsummary(metalist)
reml_meta$fdr<-p.adjust(reml_meta$p_meta,method="fdr")
reml_meta$CpG<-names(cpg_list)
write.table(reml_meta,file=paste0(model,"_",age,"_FE_meta_output.txt"),sep="\t",row.names=F,quote=F)
