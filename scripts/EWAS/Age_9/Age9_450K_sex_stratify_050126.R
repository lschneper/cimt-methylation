# title: Compiled EWAS Code
# author: Eve Schoeffler
# date: 02-28-2025
# Modified by Lisa Schneper
# date: 04-10-26

Sys.setenv(OMP_NUM_THREADS=4)
Sys.setenv(MKL_NUM_THREADS=4)
Sys.setenv(BLAS_NUM_THREADS=4)
Sys.setenv(LAPACK_NUM_THREADS=4)

.libPaths("/drives/drive1/home/lweis/R/x86_64-redhat-linux-gnu-library/4.4/")
library(optparse)

library(parallel)
library(qs2,lib.loc = "~/R/x86_64-redhat-linux-gnu-library/4.5")

args<-commandArgs(trailingOnly=TRUE)
start_cpg<-as.integer(args[1])
end_cpg<-as.integer(args[2])
# model_num<-as.integer(args[3])
sex<-as.integer(args[3])
#start_cpg<-as.integer(1)
#end_cpg<-as.integer(20)
#race<-as.integer(2)
#for (i in c(1,3,4,5,8,9)){
sex_desc<-ifelse(sex==0,"MALE",ifelse(sex==1,"FEMALE",NA))
model_num<-1
OUTPUT_FOLDER="/drives/drive1/ff_chaya/Lisa/April_2026/ewas_results/stratify/sex_stratify"
dir.create(OUTPUT_FOLDER)
OUTPUT_FILESTEM=paste0("y9_450K_sensitivity_model",model_num,"_",sex_desc,"_",start_cpg,"_",end_cpg)

# # model 1 = Basic model
model_1 <- "Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age"
designmodel_1<- "+Leukocytes+SVPC1+SVPC2+SVPC3+SVPC4+SVPC5+SVPC6+SVPC7+age"
# 
# # model 3 = Basic model + self reported race
#model_3<- "raceEthY22 sexFemale Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age"
#designmodel_3<-"+raceEthY22+sexFemale+Leukocytes+SVPC1+SVPC2+SVPC3+SVPC4+SVPC5+SVPC6+SVPC7+age"
# 
# #model 4 = Basic model + ancestry PCs
# model_4<- "sexFemale Bcell CD4T CD8T MIX Mono Neu NK bloodEPIC_SVPC1 bloodEPIC_SVPC2 bloodEPIC_SVPC3 bloodEPIC_SVPC4 ageYA_cvh GPC1 GPC2 GPC3 GPC4 GPC5"
# designmodel_4<- "+sexFemale+Bcell+CD4T+CD8T+MIX+Mono+Neu+NK+bloodEPIC_SVPC1+bloodEPIC_SVPC2+bloodEPIC_SVPC3+bloodEPIC_SVPC4+ageYA_cvh + GPC1 + GPC2 + GPC3 + GPC4 + GPC5"

#model 8 = Model 1 with self-reported race interaction
#model_8<- "raceEthY22 sexFemale Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age"
#designmodel_8<-"*raceEthY22+sexFemale+Leukocytes+SVPC1+SVPC2+SVPC3+SVPC4+SVPC5+SVPC6+SVPC7+age"

# #model 9 == Model 1 with self reported sex interaction
# model_9<- "sexFemale Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age"
# designmodel_9<- "*sexFemale+Leukocytes+SVPC1+SVPC2+SVPC3+SVPC4+SVPC5+SVPC6+SVPC7+age"

# DESIGN<-if(model_num==1){
#   model_1} else if (model_num==2){
#     model_2} else if (model_num==3){
#       model_3} else if(model_num==4){
#         model_4} else if (model_num==5){
#           model_5} else if (model_num==6){
#             model_6} else if(model_num==7) {
#               model_7} else if (model_num==8) {
#                 model_8} else if (model_num==9) {
#                   model_9} else {NA}
# DESIGN_MODEL<-if(model_num==1){
#   designmodel_1} else if (model_num==2){
#     designmodel_2} else if (model_num==3){
#       designmodel_3} else if(model_num==4){
#         designmodel_4} else if (model_num==5){
#           designmodel_5} else if (model_num==6){
#             designmodel_6} else if(model_num==7) {
#               designmodel_7} else if (model_num==8) {
#                 designmodel_8} else if (model_num==9) {
#                   designmodel_9} else {NA}
DESIGN<-model_1
DESIGN_MODEL<-designmodel_1
#DESIGN="sexFemale Bcell CD4T CD8T MIX Mono Neu NK bloodEPIC_SVPC1 bloodEPIC_SVPC2 bloodEPIC_SVPC3 bloodEPIC_SVPC4 ageYA_cvh"
design_string <- DESIGN
DESIGN <- unlist(strsplit(design_string, " "))

#survey <- readRDS("/drives/drive1/ff_chaya/Lisa/nc_samplesheets/y15_saliva_EPIC_techPC_Middleton_cimt.RDS") # y15
#survey <- readRDS("/drives/drive1/ff_chaya/Lisa/nc_samplesheets/y22_blood_031326.RDS") # y22
survey <- readRDS("/drives/drive1/ff_chaya/Lisa/nc_samplesheets/y9_450K_ses_toast_ancPC_032726.RDS")
#survey$SexFemale<-ifelse(survey$cm1bsex %in% "2 Girl",1,0)
rownames(survey) <- survey$Array_ID
# subset survey to race
survey<-survey[survey$sexFemale %in% sex,]
# Read in beta matrix
beta<-qs_read(file="/drives/drive1/ff_chaya/Lisa/nc_betafiles/BetaMatrix.k450.n1724.qcfilter.020425.qs")

# CpGs - significant at any time point
#sigcpgs<-read.delim("/drives/drive1/ff_chaya/Lisa/ewas_results/meta_output/significantCpGs_041026.txt")
# Subset beta for sigcpgs
#beta<-beta[row.names(beta) %in% sigcpgs$x,]

# Subset survey data
#covariates
survey_sub<-survey[row.names(survey)%in% colnames(beta),]
covar<-survey_sub[,names(survey_sub)%in% DESIGN]
covar<-covar[complete.cases(covar),]
# y variable
y<-survey_sub[which(row.names(survey_sub) %in% row.names(covar)),"log_Mean_Mean_CIMT"]
names(y)<-row.names(covar)

# Read in 450K manifest
man <- (readRDS("/drives/drive1/geneticData/referenceDataSets/Methylation/Zhou/hm450.hg19.manifest.dataframe.rds"))
man$ID<-row.names(man)
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

