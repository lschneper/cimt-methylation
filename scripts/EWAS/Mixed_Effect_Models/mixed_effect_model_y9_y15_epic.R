# title: Mixed Effect EWAS Model Code - random intercepts
# author: Lisa Schneper - adapting code from Eve Schoeffler and Muhammad Al Amin
# date: 01-16-2026
# 
 
Sys.setenv(OMP_NUM_THREADS=4)
Sys.setenv(MKL_NUM_THREADS=4)
Sys.setenv(BLAS_NUM_THREADS=4)
Sys.setenv(LAPACK_NUM_THREADS=4)

.libPaths("path to R package library")
library(optparse)
library(bacon)
library(qqman)
library(lme4)
library(magrittr)
library(dplyr)
library(qs2)

OUTPUT_FOLDER="path to output directory/"
OUTPUT_FILESTEM="output filestem"


y9survey <- readRDS("path to RDS file with Age 9 dependent variable, covariates, Array ID and subject ID (CVID)") # y9
y9survey$SexFemale<-ifelse(y9survey$cm1bsex %in% "2 Girl",1,0)
y9survey$wave_age<- 9

rownames(y9survey) <- y9survey$meth_id
# subset survey to only include longitudinal samples
y15survey<-readRDS("path to RDS file with Age 15 dependent variable, covariates, Array ID and subject ID (CVID)")
y15survey$SexFemale<-ifelse(y15survey$cm1bsex %in% "2 Girl",1,0)
y15survey<-y15survey[y15survey$CVID %in% y9survey$CVID,]
y15survey$wave_age<-15
row.names(y15survey)<-y15survey$meth_id
y9survey<-y9survey[y9survey$CVID %in% y15survey$CVID,]
y9surveyo<-y9survey[match(y15survey$CVID,y9survey$CVID),]
survey<-rbind(y9surveyo,y15survey)

rm(y9survey)
rm(y15survey)

# Read in EPIC manifest
man <- readRDS("path to EPIC.hg19.manifest.dataframe.rds")
man$ID<-row.names(man)
# # # Exclude bad probes and sex chromosomes
manSub <- man[man$seqnames %in% c("chrX","chrY","chrM") | !substr(man$ID,1,2) %in% "cg","ID"]

# Read in beta matrix for y9 & 15
beta <- qs_read("path to BetaMatrix.saliva.n2627.qcfilter.020425.qs") 
# Subset betamatrix by samples
beta<-beta[,colnames(beta) %in% survey$meth_id]
beta <- beta[!rownames(beta) %in% manSub,]
betao<-beta[,match(survey$meth_id,colnames(beta))]
rm(beta)
m <- log2(betao) - log2(1 - betao)
rm(betao)
gc()

m_t <- t(m)
rm(m)
gc()
##############################
cpg_lm <- function(cpg_index){
  
  cpg_name<-colnames(m_t)[cpg_index]
  
  tryCatch({
    
  model<-lmer(
    m_t[,cpg_name]~age+SexFemale+SVPC1+SVPC2+SVPC3+SVPC4+SVPC5+SVPC6+SVPC7+Leukocytes+(1|CVID),data=survey
    )
  if (isSingular(model)) {
    return(data.frame(
      probe=cpg_name,
      Estimate=NA,
      `Std. Error`=NA,
      `t value`=NA,
      `Pr(>|t|)`= NA,
      n_obs=NA,
      status="singular",
      check.names=F
    ))
  }
  ranef_df<-ranef(model)$CVID
  ranef_df$CVID<-row.names(ranef_df)
  dnAmfeatures<-ranef_df %>%
    rename(
      dnAm_intercept=`(Intercept)`
    )
  
  temp<-lm(logcimt ~ dnAm_intercept,data=left_join(y9surveyo,dnAmfeatures,by="CVID"))
  
  # Extract model_summary
  coef_summary<-summary(temp)$coefficients[2,] # Coefficients for the CpG
  data.frame(
    probe=cpg_name,
    Estimate=coef_summary["Estimate"],
    `Std. Error`= coef_summary["Std. Error"],
    `t value`= coef_summary["t value"],
      `Pr(>|t|)`= coef_summary["Pr(>|t|)"],
    n_obs=format(length(temp$fitted.values),scientific=F),
    status="ok",
    check.names=F
  )
  }, error = function(e) {
    data.frame(
      probe=cpg_name,
      Estimate=NA,
      `Std. Error`=NA,
      `t value`=NA,
      `Pr(>|t|)`= NA,
      n_obs=NA,
      status="error",
      check.names=F
    )
  })
}
      

chunk <- 1:ncol(m_t)
#chunk<-1:4
modelresults <- list()
modelresults <- lapply(chunk, cpg_lm)

# Need to rename the sixth column because some are n and some are n_obs
library(data.table)
DT<-rbindlist(
  lapply(modelresults, function(x) {
    setnames(x, 6, "n_obs")
    x
  }),
  use.names=TRUE,
  fill=TRUE
)
 
modelresults_df<-DT[substr(DT$probe,1,2) %in% "cg" & DT$status=="ok",]
gc()
######
modelresults_df$fdr<-p.adjust(modelresults_df$P.Value,method="fdr")

bc <- bacon(modelresults_df$t,effectsizes=modelresults_df$Estimate,standarderrors=modelresults_df$`Pr(>|t|)`)
inflation(bc)

modelresults_df$bacon.P.Value <- as.numeric(pval(bc)[,1])
modelresults_df$bacon.est <- as.numeric(es(bc)[,1])
modelresults_df$bacon.se <- as.numeric(se(bc)[,1])
modelresults_df<-as.data.frame(modelresults_df)
row.names(modelresults_df)<-modelresults_df$probe
man_hits <- man[rownames(modelresults_df), ]
man_hitso<-man_hits[match(modelresults_df$probe,man_hits$ID),]
modelresults_df$Gene <- man_hitso$gene
modelresults_df$Chr <- as.character(man_hitso$chrm_A)
modelresults_df$Start <- man_hitso$start
modelresults_df$End <- man_hitso$end
modelresults_df$fdr<-p.adjust(modelresults_df$`Pr(>|t|)`,method="fdr")
modelresults_df$bacon.fdr<-p.adjust(modelresults_df$bacon.P.Value,method="fdr")
head(modelresults_df)
table(modelresults_df$fdr<0.05)
table(modelresults_df$bacon.fdr<0.05)
write.table(modelresults_df, file=paste0(OUTPUT_FOLDER, OUTPUT_FILESTEM, ".ewasoutput.txt"), sep = "\t", quote = F, row.names = F)

bed <- modelresults_df[, c("Chr", "Start", "End", "bacon.P.Value")]
bed <- bed[!is.na(bed$Chr), ]
bed <- bed[order(bed$Chr, bed$Start), ]

names(bed) <- c("chrom", "start", "end", "P.Value")
bed <- bed[order(bed$chrom, bed$start), ]
write.table(bed, paste0(OUTPUT_FOLDER, OUTPUT_FILESTEM, ".baconadj.combpinput.bed"), quote = F, row.names = F, sep = "\t")
lambda_fun <- function(pval) { 
  z <- qnorm(pval/2, lower.tail = F)
  gc <- median(z^2, na.rm = T)/qchisq(0.5,1)
  return(gc)
}
#####
bac_eif<-lambda_fun(modelresults_df$bacon.P.Value)
bac_eif
eif<-lambda_fun(modelresults_df$`Pr(>|t|)`)
eif

write.table(as.data.frame(DT[DT$status %in% "singular","probe"]),file=paste0(OUTPUT_FOLDER, OUTPUT_FILESTEM, ".singularprobes.txt"), quote = F, row.names = F, sep = "\t")
table(DT$status)
