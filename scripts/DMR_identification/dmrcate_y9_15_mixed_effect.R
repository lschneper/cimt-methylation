.libPaths("path to R package library")
library(DMRcate)
library(bacon)
library(GenomicRanges)

bacon_fun<-function(df){
  bc<-bacon(teststatistics = df$t,effectsizes = df$Coefficient,standarderrors = df$StdError)
  df$bacon.P.Value <- as.numeric(pval(bc)[,1])
  df$bacon.est <- as.numeric(es(bc)[,1])
  df$bacon.se <- as.numeric(se(bc)[,1])
  df$bacon.t<-as.numeric(tstat(bc)[,1])
  return(df)
}


mod1<-readRDS("path to y9y15_mixedeffect_model1_meta_dfoutput.RDS")
mod1$z_fixed<-mod1$TE_fixed/mod1$seTE_fixed
bc<-bacon(teststatistics = mod1$z_fixed,effectsizes = mod1$TE_fixed,standarderrors = mod1$seTE_fixed)
mod1$bacon.P.Value_fixed <- as.numeric(pval(bc)[,1])
mod1$bacon.est_fixed <- as.numeric(es(bc)[,1])
mod1$bacon.se_fixed <- as.numeric(se(bc)[,1])
mod1$bacon.t_fixed<-as.numeric(tstat(bc)[,1])

man <- readRDS("path to EPIC.hg19.manifest.dataframe.rds")
man$ID<-row.names(man)
man<-man[,1:3]

mod1b<-merge(man,mod1,by.x="row.names",by.y="probe")
mod1b$bacon.fdr_fixed<-p.adjust(mod1b$bacon.P.Value_fixed,method="fdr")
mod1b$is.sig<-mod1b$bacon.fdr_fixed<0.05
mod1b$seqnames<-as.character(mod1b$seqnames)

probenames<-as.character(mod1b$Row.names)
dmrcateinput<-GRanges(
  seqnames=mod1b$seqnames,
  ranges=IRanges(start=mod1b$start,end=mod1b$start),
  diff=mod1b$bacon.est_fixed,
  stat=mod1b$bacon.t_fixed,
  ind.fdr=mod1b$bacon.fdr_fixed,
  is.sig=mod1b$is.sig
)


cleandiff<-as.numeric(mod1b$bacon.est_fixed)
cleant<-as.numeric(mod1b$bacon.t_fixed)
cleanfdr<-as.numeric(mod1b$bacon.fdr_fixed)

dmrcateinput<-GRanges(
  seqnames=mod1b$seqnames,
  ranges=IRanges(start=mod1b$start,width=1),
  ID=probenames,
  stat=cleant,  
  diff=cleandiff,
  rawpval=mod1b$bacon.P.Value_fixed,
  ind.fdr=cleanfdr,
  is.sig=cleanfdr<0.05
    )
names(dmrcateinput)<-mcols(dmrcateinput)$ID
library(qs2)
betas<-qs_read("path to BetaMatrix.saliva.n2627.qcfilter.020425.qs")

betas<-betas[row.names(betas) %in% probenames,]
betaso<-betas[match(probenames,row.names(betas)),]

annot_custom_beta<-new("CpGannotated",ranges=dmrcateinput,betas=betaso)
dmrcate_results<-dmrcate(annot_custom,lambda = 1000,C=2)
dmrcateresultdf<-data.frame(
  coord=dmrcate_results@coord,
  no.cpgs=dmrcate_results@no.cpgs,
  min_smoothed_fdr=dmrcate_results@min_smoothed_fdr,
  Stouffer=dmrcate_results@Stouffer,
  HMFDR=dmrcate_results@HMFDR,
  Fisher=dmrcate_results@Fisher,
  maxdiff=dmrcate_results@maxdiff,
  meandiff=dmrcate_results@meandiff
)


dir.create("output directory path")
saveRDS(dmrcate_results,"output file path")
saveRDS(dmrcateresultdf,"ouitput file path (should differ from above)")
