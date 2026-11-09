.libPaths("path to R  package library")
library(DMRcate)
library(bacon)
library(GenomicRanges)
library(qs2)

bacon_fun<-function(df){
  bc<-bacon(teststatistics = df$z_fixed,effectsizes = df$TE_fixed,standarderrors = df$seTE_fixed)
  df$bacon.P.Value_fixed <- as.numeric(pval(bc)[,1])
  df$bacon.est_fixed <- as.numeric(es(bc)[,1])
  df$bacon.se_fixed <- as.numeric(se(bc)[,1])
  df$bacon.t_fixed<-as.numeric(tstat(bc)[,1])
  return(df)
}


mod1b<-readRDS("path to y15_meta_sensitivity_model1_meta_analysis_output_bacon_chipseeker.RDS")
mod1b<-bacon_fun(mod1b)

mod1b$is.sig<-mod1b$bacon.fdr_fixed<0.05
mod1b$seqnames<-as.character(mod1b$seqnames)

probenames<-as.character(mod1b$Name)
dmrcateinput<-GRanges(
  seqnames=mod1b$seqnames,
  ranges=IRanges(start=mod1b$start,end=mod1b$start),
  ID=probenames,
  stat=mod1b$bacon.t_fixed,
  diff=mod1b$bacon.est_fixed,
  rawpval=mod1b$bacon.P.Value_fixed,
  ind.fdr=mod1b$bacon.fdr_fixed,
  is.sig=mod1b$is.sig
)

names(dmrcateinput)<-mcols(dmrcateinput)$ID

betas<-qs_read("path to BetaMatrix.saliva.n2627.qcfilter.020425.qs")

betas<-betas[row.names(betas) %in% probenames,]
betaso<-betas[match(probenames,row.names(betas)),]

annot_custom<-new("CpGannotated",ranges=dmrcateinput,betas=betaso)

dmrcate_results<-dmrcate(annot_custom,lambda = 1000,C=2)
dir.create("output directory")
saveRDS(dmrcate_results,"output file")

#### y9###
mod1b<-readRDS("path to y9_meta_sensitivity_model1_output.RDS")
mod1b<-bacon_fun(mod1b)

mod1b$is.sig<-mod1b$bacon.fdr_fixed<0.05
mod1b$seqnames<-as.character(mod1b$seqnames)

probenames<-as.character(mod1b$Name)
dmrcateinput<-GRanges(
  seqnames=mod1b$seqnames,
  ranges=IRanges(start=mod1b$start,end=mod1b$start),
  ID=probenames,
  stat=mod1b$bacon.t_fixed,
  diff=mod1b$bacon.est_fixed,
  rawpval=mod1b$bacon.P.Value_fixed,
  ind.fdr=mod1b$bacon.fdr_fixed,
  is.sig=mod1b$is.sig
)

names(dmrcateinput)<-mcols(dmrcateinput)$ID
betas<-qs_read("path to BetaMatrix.saliva.n2627.qcfilter.020425.qs")

betas<-betas[row.names(betas) %in% probenames,]
betaso<-betas[match(probenames,row.names(betas)),]

annot_custom<-new("CpGannotated",ranges=dmrcateinput,betas=betaso)

dmrcate_results<-dmrcate(annot_custom,lambda = 1000,C=2)
dir.create("/output file path")
saveRDS(dmrcate_results,"output file")
