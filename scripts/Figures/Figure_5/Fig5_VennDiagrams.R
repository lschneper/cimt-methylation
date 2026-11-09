library(eulerr)
library(ggVennDiagram)
library(ggplot2)
library(Cairo)
library(GenomicRanges)

setwd("path to working directory")
OUTDIR=getwd()
long_dmrs<-as.data.frame(readRDS("mixed effect model DMRs"))
long_dmrs_gr<-makeGRangesFromDataFrame(long_dmrs,keep.extra.columns = T)
y22_dmrs<-readRDS("Age 22 DMRs")
y22_dmrs_gr<-makeGRangesFromDataFrame(y22_dmrs,keep.extra.columns = T)

# dmr1 dmr2
create_df_for_Venn_fun<-function(dmr1,dmr2,setnameDMR1,setnameDMR2){
dmr1_gr<-makeGRangesFromDataFrame(dmr1)
dmr2_gr<-makeGRangesFromDataFrame(dmr2)

hits<-findOverlaps(dmr1_gr,dmr2_gr)

dmr1_idx<-unique(queryHits(hits))
dmr2_idx<-unique(subjectHits(hits))
# Counts
n1_only<-length(dmr1_gr)-length(dmr1_idx)
n2_only<-length(dmr2_gr)-length(dmr2_idx)
n_overlap1<-length(dmr1_idx)
n_overlap2<-length(dmr2_idx)

df<-data.frame(
  Set=c(setnameDMR1,setnameDMR2),
  Total=c(length(dmr1_gr),length(dmr2_gr)),
  Overlapping=c(n_overlap1,n_overlap2),
  Unique=c(n1_only,n2_only)
)
return(df)
}

y22_long_df<-create_df_for_Venn_fun(dmr1=y22_dmrs,dmr2=long_dmrs,setnameDMR1 = "Age 22",setnameDMR2 = "Age9 - Age 15")

create_vennlist<-function(dmr1,dmr2,setDMRname1,setDMRname2){
  dmr1_gr<-makeGRangesFromDataFrame(dmr1)
  dmr2_gr<-makeGRangesFromDataFrame(dmr2)
  combined<-GenomicRanges::reduce(c(dmr1_gr,dmr2_gr),ignore.strand=TRUE,min.gapwidth=0L)
  
  inr1<-overlapsAny(combined,dmr1_gr,ignore.strand=TRUE)
  inr2<-overlapsAny(combined,dmr2_gr,ignore.strand=TRUE)
  
  venn_list<-list(
    setname1=which(inr1),
    setname2=which(inr2)
    )
  names(venn_list)<-c(setDMRname1,setDMRname2)
  return(venn_list)

}

y22_long_vennlist<-create_vennlist(long_dmrs,y22_dmrs,"Ages\n9 to 15","Age\n22")

CairoPDF(file=paste0(OUTDIR,"/Fig5_DMRvenn_y22_long_paneld.pdf"),family="Arial",height=8,width=8)
plot(venn(y22_long_vennlist), quantities = list(cex=1.75), labels = list(cex=2),fills = c("white","#CCC","lightblue"))
dev.off()

persistent<-read.delim("persistent_dmrs_from_combp_dmrcate_for_venn.txt")
persistent_long_vennlist<-create_vennlist(long_dmrs,persistent,"Ages\n9 to 15","Persistent\nDMRs")
CairoPDF(file=paste0(OUTDIR,"/Fig5_DMRvenn_persistent_long_panele.pdf"),family="Arial",height=8,width=8)
plot(venn(persistent_long_vennlist), quantities = list(cex=1.75), labels = list(cex=2),fills = c("white","#CCC","lightblue"))
dev.off()

### DMPs
y22_epic_dmps=readRDS("Age 22 epic results")
y22_epic_dmps<-y22_epic_dmps[y22_epic_dmps$bacon.fdr<0.05,]
long_meta<-readRDS("Mixed effect model ewas results")
names(long_meta)
library(bacon)
bacon_fun<-function(df){
  df$z_fixed<-df$TE_fixed/df$seTE_fixed
  bc<-bacon(teststatistics = df$z,effectsizes = df$TE_fixed,standarderrors = df$seTE_fixed)
  df$bacon.P.Value_fixed <- as.numeric(pval(bc)[,1])
  df$bacon.est_fixed <- as.numeric(es(bc)[,1])
  df$bacon.se_fixed <- as.numeric(se(bc)[,1])
  return(df)
}
long_meta<-bacon_fun(long_meta)
long_meta$bacon.fdr_fixed<-p.adjust(long_meta$bacon.P.Value_fixed,method="BH")
long_meta<-long_meta[long_meta$bacon.fdr_fixed<0.05,]
common<-read.table("persistent_ewas_hits.txt",sep="\t")

dmp_long_y22<-list("Ages\n9 to 15"=long_meta$probe,"Age\n22"=y22_epic_dmps$probe)
CairoPDF(file=paste0(OUTDIR,"/Fig5_DMPvenn_Age22_long_panelb.pdf"),family="Arial",height=8,width=8)
plot(venn(dmp_long_y22), quantities = list(cex=1.75), labels = list(cex=2),fills = c("white","#CCC","lightblue"))
dev.off()

dmp_long_common<-list("Ages\n9 to 15"=long_meta$probe,"Persistent\nDMPs"=common$probe)
CairoPDF(file=paste0(OUTDIR,"/Fig5_DMPvenn_common_long_panelc.pdf"),family="Arial",height=8,width=8)
plot(venn(dmp_long_common), quantities = list(cex=1.75), labels = list(cex=2),fills = c("white","#CCC","lightblue"))
dev.off()
writeLines(capture.output(sessionInfo()),"/drives/drive1/ff_chaya/Lisa/github_code/Figures/Figure_5/VENN.sessionInfo.txt")
