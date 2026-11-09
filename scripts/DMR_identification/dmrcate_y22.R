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


mod1<-read.delim("path to y22_epic_model1.ewasoutput.txt")
mod1a<-bacon_fun(mod1)

man <- readRDS("path to EPIC.hg19.manifest.dataframe.rds")
man$ID<-row.names(man)
man<-man[,1:3]

mod1b<-merge(man,mod1a,by.x="row.names",by.y="probe")
mod1b$bacon.fdr<-p.adjust(mod1b$bacon.P.Value,method="fdr")
mod1b$is.sig<-mod1b$bacon.fdr<0.05
mod1b$seqnames<-as.character(mod1b$seqnames)

probenames<-as.character(mod1b$Row.names)
y22dmrcateinput<-GRanges(
  seqnames=mod1b$seqnames,
  ranges=IRanges(start=mod1b$start,end=mod1b$start),
  diff=mod1b$bacon.est,
  stat=mod1b$bacon.t,
  ind.fdr=mod1b$bacon.fdr,
  is.sig=mod1b$is.sig
)


cleandiff<-as.numeric(mod1b$bacon.est)
cleant<-as.numeric(mod1b$bacon.t)
cleanfdr<-as.numeric(mod1b$bacon.fdr)

y22dmrcateinput<-GRanges(
  seqnames=mod1b$seqnames,
  ranges=IRanges(start=mod1b$start,width=1),
  ID=probenames,
  stat=cleant,  
  diff=cleandiff,
  rawpval=mod1b$bacon.P.Value,
  ind.fdr=cleanfdr,
  is.sig=cleanfdr<0.05
    )
names(y22dmrcateinput)<-mcols(y22dmrcateinput)$ID
library(qs2)
betas<-qs_read("path to BetaMatrix.blood.n1271.qcfilter.020425.qs")

betas<-betas[row.names(betas) %in% probenames,]
betaso<-betas[match(probenames,row.names(betas)),]

annot_custom<-new("CpGannotated",ranges=y22dmrcateinput,betas=betaso)

y22_dmrcate_results<-dmrcate(annot_custom,lambda = 1000,C=2)
dir.create("output file directory path")
saveRDS(y22_dmrcate_results,"output file path and name")

