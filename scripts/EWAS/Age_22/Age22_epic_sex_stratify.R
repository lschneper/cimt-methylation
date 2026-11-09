# author: Lisa Schneper
# date: 05-01-26

Sys.setenv(OMP_NUM_THREADS=4)
Sys.setenv(MKL_NUM_THREADS=4)
Sys.setenv(BLAS_NUM_THREADS=4)
Sys.setenv(LAPACK_NUM_THREADS=4)

.libPaths("path to R package library")
library(optparse)

library(parallel)
library(qs2)

args<-commandArgs(trailingOnly=TRUE)
start_cpg<-as.integer(args[1])
end_cpg<-as.integer(args[2])
sex<-as.integer(args[3])

sex_desc<-ifelse(sex==0,"MALE",ifelse(sex==1,"FEMALE",NA))
model_num<-1
OUTPUT_FOLDER="path to output directory"
dir.create(OUTPUT_FOLDER)
OUTPUT_FILESTEM=paste0("y22_epic_sensitivity_model",model_num,"_",sex_desc,"_",start_cpg,"_",end_cpg)

# # model 1 = Basic model
model_1 <- "Bcell CD4T CD8T MIX Mono Neu NK bloodEPIC_SVPC1 bloodEPIC_SVPC2 bloodEPIC_SVPC3 bloodEPIC_SVPC4 ageYA_cvh"
designmodel_1<- "+Bcell+CD4T+CD8T+MIX+Mono+Neu+NK+bloodEPIC_SVPC1+bloodEPIC_SVPC2+bloodEPIC_SVPC3+bloodEPIC_SVPC4+ageYA_cvh"
# 
DESIGN<-model_1
DESIGN_MODEL<-designmodel_1

design_string <- DESIGN
DESIGN <- unlist(strsplit(design_string, " "))

survey <- readRDS("path to RDS file with Age 22 dependent variable, covariates and Array ID")
#survey$SexFemale<-ifelse(survey$cm1bsex %in% "2 Girl",1,0)
rownames(survey) <- survey$Array_ID
# subset survey to race
survey<-survey[survey$sexFemale %in% sex,]
# Read in beta matrix
beta<-qs_read(file="path to BetaMatrix.blood.n1271.qcfilter.020425.qs")

# Subset survey data
#covariates
survey_sub<-survey[row.names(survey)%in% colnames(beta),]
covar<-survey_sub[,names(survey_sub)%in% DESIGN]
covar<-covar[complete.cases(covar),]
# y variable
y<-survey_sub[which(row.names(survey_sub) %in% row.names(covar)),"log_Mean_Mean_CIMT"]
names(y)<-row.names(covar)

# Read in EPIC manifest
man <- (readRDS("path to EPIC.hg19.manifest.rds"))
man$ID<-names(man)
man<-as.data.frame(man)
# # # Exclude bad probes and sex chromosomes
manSub <- man[man$seqnames %in% c("chrX","chrY","chrM") | !substr(man$ID,1,2) %in% "cg","ID"]

#
# Subset betamatrix by samples and remove nonCpG probes, probes mapping to X,Y and M
beta<-beta[,colnames(beta) %in% row.names(covar)]
beta <- beta[!rownames(beta) %in% manSub,]
betao<-beta[,match(row.names(covar),colnames(beta))]
end_cpg[end_cpg > nrow(betao)] <- nrow(betao)
betao<-betao[start_cpg:end_cpg,]

rm(beta)
m <- log2(betao) - log2(1 - betao)
rm(betao)
gc()


m_t <- t(m)

rm(m)
gc()

cpg_lm <- function(cpg_index){
  # Extract the CpG name
  cpg_name<-colnames(m_t)[cpg_index]
  design_m <- cbind(m_t[,cpg_index], covar)
  colnames(design_m)[1]<-cpg_name
  form<-as.formula(paste0("y~",cpg_name,DESIGN_MODEL))
  temp <- try(lm(form,data=design_m), silent=TRUE)
  if (inherits(temp,"try-error")){
    return(rep(NA,10))
  }
  # Extract model_summary
  coef_table<-summary(temp)$coefficients # Coefficients for the CpG orig
  df<-data.frame(probe=cpg_name,
  beta=coef_table[cpg_name,"Estimate"],
  t=coef_table[cpg_name,"t value"],
  se=coef_table[cpg_name,"Std. Error"],
  p=coef_table[cpg_name,"Pr(>|t|)"],
  n = nobs(temp)
  )
  return(df)
}
  


cpg_indices<-1:ncol(m_t)

modelresults_list<-lapply(1:ncol(m_t), cpg_lm)
modelresultsdf<-do.call(rbind,modelresults_list)
write.table(modelresultsdf, file=paste0(OUTPUT_FOLDER, OUTPUT_FILESTEM,".ewasoutput.txt"),sep="\t",quote=F,row.names=F)
rm(modelresults_list)
rm(m_t)
gc()

quit(save="no")

