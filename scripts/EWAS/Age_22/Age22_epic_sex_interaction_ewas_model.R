# author: Lisa Schneper
# date: 04-10-26

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

model_num<-9
OUTPUT_FOLDER="path to output directory"
dir.create(OUTPUT_FOLDER)
OUTPUT_FILESTEM=paste0("y22_epic_sensitivity_model",model_num,"_",start_cpg,"_",end_cpg)

# # model 1 = Basic model
model_1 <- "sexFemale Bcell CD4T CD8T MIX Mono Neu NK bloodEPIC_SVPC1 bloodEPIC_SVPC2 bloodEPIC_SVPC3 bloodEPIC_SVPC4 ageYA_cvh"
designmodel_1<- "+sexFemale+Bcell+CD4T+CD8T+MIX+Mono+Neu+NK+bloodEPIC_SVPC1+bloodEPIC_SVPC2+bloodEPIC_SVPC3+bloodEPIC_SVPC4+ageYA_cvh"
# 
# #model 9 == Model 1 with self reported sex interaction
model_9<- "sexFemale sexFemale Bcell CD4T CD8T MIX Mono Neu NK bloodEPIC_SVPC1 bloodEPIC_SVPC2 bloodEPIC_SVPC3 bloodEPIC_SVPC4 ageYA_cvh"
designmodel_9<- "*sexFemale+Bcell+CD4T+CD8T+MIX+Mono+Neu+NK+bloodEPIC_SVPC1+bloodEPIC_SVPC2+bloodEPIC_SVPC3+bloodEPIC_SVPC4+ageYA_cvh"

DESIGN<-model_9
DESIGN_MODEL<-designmodel_9
DESIGN_MODEL1<-designmodel_1

design_string <- DESIGN
DESIGN <- unlist(strsplit(design_string, " "))

survey <- readRDS("path to RDS file with Age 22 dependent variable, covariates and Array ID")
#survey$SexFemale<-ifelse(survey$cm1bsex %in% "2 Girl",1,0)
rownames(survey) <- survey$Array_ID

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
  form2<-as.formula(paste0("y~",cpg_name,DESIGN_MODEL1))
   reduced_model<- try(lm(form2,data=design_m), silent=TRUE)

#   if (inherits(temp,"try-error")){
#     return(rep(NA,10))
#   }
     p_int_global<-anova(reduced_model,temp)$`Pr(>F)`[2]
#     test<-list(mod=temp,reducemod=reduced_model,anova=p_int_global)
#     return(test)

  # Extract model_summary
  coef_table<-summary(temp)$coefficients
  vc<-vcov(temp)
  main_row<-cpg_name
  int_rows<-grep(paste0("^",cpg_name, ":"),row.names(coef_table),value=TRUE)
  rows_to_keep<-c(main_row,int_rows)

  vc_sub<-vc[rows_to_keep,rows_to_keep]
  coef_Sub<-coef_table[rows_to_keep,,drop=FALSE]

  result<-list(
    CpG=cpg_name,
  beta_main=coef_Sub[main_row,"Estimate"],
  se_main = coef_Sub[main_row,"Std. Error"],
  p_main =coef_Sub[main_row,"Pr(>|t|)"],

  beta_int=coef_Sub[int_rows,"Estimate"],
  se_int = coef_Sub[int_rows,"Std. Error"],
  p_int = coef_Sub[int_rows,"Pr(>|t|)"],
  p_int_global = p_int_global,
  int_terms=int_rows,

  coef_sub=coef_Sub,
  vcov_sub=vc_sub,
  n=nobs(temp)
  )
  return(result)
}
  

cpg_indices<-1:ncol(m_t)

modelresults_list<-lapply(1:ncol(m_t), cpg_lm)
saveRDS(modelresults_list, file=paste0(OUTPUT_FOLDER, OUTPUT_FILESTEM,".ewasoutput.RDS"))
rm(modelresults_list)
rm(m_t)
gc()

quit(save="no")

