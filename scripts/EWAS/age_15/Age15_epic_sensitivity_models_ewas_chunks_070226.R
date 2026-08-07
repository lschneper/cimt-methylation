# title: Compiled EWAS Code
# author: Eve Schoeffler
# date: 02-28-2025
# Modified by Lisa Schneper
# date: 03-28-26

Sys.setenv(OMP_NUM_THREADS=4)
Sys.setenv(MKL_NUM_THREADS=4)
Sys.setenv(BLAS_NUM_THREADS=4)
Sys.setenv(LAPACK_NUM_THREADS=4)

.libPaths("~/R/x86_64-redhat-linux-gnu-library/4.4/")
library(optparse)
library(bacon)
library(qqman)
library(parallel)
library(qs2,lib.loc = "~/R/x86_64-redhat-linux-gnu-library/4.5")

args<-commandArgs(trailingOnly=TRUE)
start_cpg<-as.integer(args[1])
end_cpg<-as.integer(args[2])
model_num<-as.integer(args[3])

#start_cpg<-as.integer(1)
#end_cpg<-as.integer(20)
#model_num<-2

OUTPUT_FOLDER="PATH_TO_OUTPUT_DIRECTORY"
dir.create(OUTPUT_FOLDER)
OUTPUT_FILESTEM=paste0("y15_epic_sensitivity_",start_cpg,"_",end_cpg,"_model",model_num)

# model 1 = Basic model
model_1<- "sexFemale Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age"

# model2 = Toast model
model_2<- "sexFemale toast_1 toast_2 toast_3 toast_4 toast_5 SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7  age"

# model 3 = Basic model + self reported race
model_3<- "sexFemale Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age raceEthY22"

# model 4 = Basic model + ancestry PCs
model_4<- "sexFemale Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age GPC1 GPC2 GPC3 GPC4 GPC5"

# model 5 = Basic model + smoking
model_5<- "sexFemale Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age y15_eversmoke"

# model 6 = Basic model + mother education
model_6<- "sexFemale Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age eduMom_yrs"

# model 7 = Basic model + mother education + smoking
model_7<- "sexFemale Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age eduMom_yrs y15_eversmoke"

# model 10 = Basic model + fetal exposure
model_10<- "sexFemale Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age y15_eversmoke_fetalexp"

# model 11 = Basic model + Clinic city
model_11<- "sexFemale Leukocytes SVPC1 SVPC2 SVPC3 SVPC4 SVPC5 SVPC6 SVPC7 age City"

DESIGN<-if(model_num==1){
  model_1} else if (model_num==2){
    model_2} else if (model_num==3){
      model_3} else if(model_num==4){
        model_4} else if (model_num==5){
          model_5} else if (model_num==6){
            model_6} else if(model_num==7) {
              model_7} else if (model_num==10) {
                model_10} else if (model_num==11) {
                  model_11} else {NA}

design_string <- DESIGN
DESIGN <- unlist(strsplit(design_string, " "))

survey <- readRDS("PATH_TO_VARIABLE_FILE")

rownames(survey) <- survey$Array_ID

# Read in beta matrix
beta<-qs_read(file="PATH_TO_BETA_MATRIX")

# Subset survey data
#covariates
survey_sub<-survey[row.names(survey)%in% colnames(beta),]
covar<-survey_sub[,names(survey_sub)%in% DESIGN]
covar<-covar[complete.cases(covar),]
# y varialbe
y<-survey_sub[which(row.names(survey_sub) %in% row.names(covar)),"log_Mean_Mean_CIMT"]
names(y)<-row.names(covar)

# Read in EPIC manifest
man <- readRDS("PATH_TO_MANIFEST")
man$ID<-row.names(man)
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
  temp <- lm(y ~.,data=design_m)
  
  # Extract model_summary
  coef_summary<-summary(temp)$coefficients[2,] # Coefficients for the CpG
  n_obs <- format(length(temp$fitted.values),scientific=F)
  # Combine and name the results
  result<-c(coef_summary,n=n_obs)
  return(result)
}

# library(parallel)
chunk_size<-100000
cpg_indices<-1:ncol(m_t)
#chunk <- 1:ncol(m_t) # original
# chunk<-1:1000
chunks<-split(cpg_indices,ceiling(seq_along(cpg_indices)/chunk_size))
ncores<-1

modelresults_list<-mclapply(chunks, function(chunk){
  cat("Began processing chunk with", length(chunk), "CpGs at ",Sys.time(),"\n")
  chunk_results<-lapply(chunk, cpg_lm)
  cat("Finished processing chunk with", length(chunk), "CpGs at ",Sys.time(),"\n")
  gc()
  do.call(rbind,chunk_results)
}, mc.cores=ncores)

modelresults_df <- as.data.frame(do.call(rbind, modelresults_list))
names(modelresults_df) <- c("Coefficient", "StdError", "t", "P.Value","n")
for(i in 1:ncol(modelresults_df)){
  modelresults_df[,i]<-as.numeric(modelresults_df[,i])
}
modelresults_df$probe <- colnames(m_t)

write.table(modelresults_df, file=paste0(OUTPUT_FOLDER, OUTPUT_FILESTEM,"_model",model_num,".ewasoutput.txt"), sep = "\t", quote = F, row.names = F)
writeLines(capture.output(sessionInfo()),paste0(OUTPUT_FOLDER, OUTPUT_FILESTEM,"_model",model_num,".ewasoutput_sessionInfo.txt"))
quit(save="no")
