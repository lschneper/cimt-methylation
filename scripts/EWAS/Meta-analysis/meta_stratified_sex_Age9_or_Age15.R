# title: meta_NHB_stratified.R
# author: Lisa Schneper
# date: 05-04-2026

.libPaths("path to R libary packages")

library(metafor)
library(optparse)

extract_remlsummary<-function(res_list){
  do.call(rbind, lapply(res_list,function(x){
    if(is.null(x)|| all(is.na(x))) return(NULL)
    
    data.frame(
      beta_NHB = x$b[1],
      beta_race = x$b[2],
      beta_array=x$b[3],
      
      se_NHB = x$se[1],
      se_race = x$se[2],
      se_array=x$se[3],
      
      p_NHB = x$pval[1],
      p_race = x$pval[2],
      p_array=x$pval[3]
      
    )
  }))
}
cleandf_fun<-function(df){
  df<-df[!duplicated(df),]
  df<-df[!df$probe%in% "probe",]
  df$beta<-as.numeric(df$beta)
  df$se<-as.numeric(df$se)
  return(df)
}
setwd("path to sex_stratify working directory")
args<-commandArgs(trailingOnly=TRUE)
#modelnum=9
#agenum=15

agenum=as.integer(args[1])
age<-paste0("y",agenum)

file1<-read.delim(paste0("sex_stratify_",age,"_450K_MALE.ewasoutput.txt"))
file1$sex<-"MALE"
file1$array="450K"

file2<-read.delim(paste0("sex_stratify_",age,"_epic_MALE.ewasoutput.txt"))
file2$sex<-"MALE"
file2$array="epic"

file3<-read.delim(paste0("sex_stratify_",age,"_450K_FEMALE.ewasoutput.txt"))
file3$sex<-"FEMALE"
file3$array="450K"

file4<-read.delim(paste0("sex_stratify_",age,"_epic_FEMALE.ewasoutput.txt"))
file4$sex<-"FEMALE"
file4$array="epic"


file1<-cleandf_fun(file1)
file2<-cleandf_fun(file2)
file3<-cleandf_fun(file3)
file4<-cleandf_fun(file4)
                                                      
file1_2_cpg<-file1$probe[file1$probe %in% file2$probe]
file1_2_3_cpg<-file1_2_cpg[file1_2_cpg %in% file3$probe]
file1_2_3_4_cpg<-file1_2_3_cpg[file1_2_3_cpg %in% file4$probe]
file1sub<-file1[file1$probe %in% file1_2_3_4_cpg,]
file2sub<-file2[file2$probe %in% file1_2_3_4_cpg,]
file3sub<-file3[file3$probe %in% file1_2_3_4_cpg,]
file4sub<-file4[file4$probe %in% file1_2_3_4_cpg,]

both_arrays<-rbind(file1sub,file2sub,file3sub,file4sub)
both_arrays$sex<-factor(both_arrays$sex,level=c("MALE","FEMALE"))
both_arrays$array<-factor(both_arrays$array,level=c("epic","450K"))
both_arrays<-both_arrays[order(both_arrays$probe,both_arrays$sex,both_arrays$array),]

cpg_list<-split(both_arrays,both_arrays$probe)

# Run meta analysis
metalist<-list()
cat("Began processing chunk with", length(cpg_list), "CpGs at ",Sys.time(),"\n")
metalist <- lapply(cpg_list,function(x){
  rma(
    yi=x$beta,
    sei=x$se,
    mods = ~x$sex + x$array,
    method="FE"
  )
})
cat("Finished processing ", length(cpg_list), "CpGs at ",Sys.time(),"\n")

reml_meta<-extract_remlsummary(metalist)
reml_meta$fdr<-p.adjust(reml_meta$p_NHB,method="fdr")
reml_meta$CpG<-names(cpg_list)
write.table(reml_meta,file=paste0("sexMale_Female_stratified_",age,"_meta_FE_output.txt"),sep="\t",row.names=F,quote=F)
