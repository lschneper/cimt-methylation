library(eulerr)
library(ggVennDiagram)

setwd("/path to working directory")

y9<-readRDS("path to age 9 ewas results")
y9$probe<-y9$Name
y15<-readRDS("path to age 15 ewas results")
y15$probe<-y15$Name
y22<-readRDS("path to age 22 ewas results")
threshold=0.05
top9<-y9[y9$bacon.fdr_fixed<threshold,]
top15<-y15[y15$bacon.fdr_fixed< threshold,]
top22<-y22[y22$bacon.fdr< threshold,]

top9up<-y9[y9$bacon.fdr_fixed<threshold & y9$TE_fixed>0,]
top15up<-y15[y15$bacon.fdr_fixed< threshold & y15$TE_fixed>0,]
top22up<-y22[y22$bacon.fdr< threshold & y22$bacon.est>0,]

top9down<-y9[y9$bacon.fdr_fixed<threshold & y9$TE_fixed<0,]
top15down<-y15[y15$bacon.fdr_fixed< threshold & y15$TE_fixed<0,]
top22down<-y22[y22$bacon.fdr< threshold & y22$bacon.est<0,]

top9up_prom<-y9[!is.na(y9$overlapping_promoterSymbols) & y9$bacon.fdr_fixed<threshold & y9$TE_fixed>0,]
top15up_prom<-y15[!is.na(y15$overlapping_promoterSymbols) & y15$bacon.fdr_fixed< threshold & y15$TE_fixed>0,]
top22up_prom<-y22[!is.na(y22$overlapping_promoterSymbols) & y22$bacon.fdr< threshold & y22$bacon.est>0,]

top9down_prom<-y9[!is.na(y9$overlapping_promoterSymbols) & y9$bacon.fdr_fixed<threshold & y9$TE_fixed<0,]
top15dowm_prom<-y15[!is.na(y15$overlapping_promoterSymbols) & y15$bacon.fdr_fixed< threshold & y15$TE_fixed<0,]
top22down_prom<-y22[!is.na(y22$overlapping_promoterSymbols) & y22$bacon.fdr< threshold & y22$bacon.est<0,]


eulerprep_fun<-function(top9,top15,top22){
# How many overlap across all three?
overlap_915 <- intersect(top9$probe, top15$probe)
overlap_922 <- intersect(top9$probe, top22$probe)
overlap_1522 <- intersect(top15$probe, top22$probe)
overlap_all <- intersect(overlap_915, top22$probe) # 42 probes

simple_table <- top9[top9$probe %in% overlap_all, c("seqnames", "start", "end")]
names(simple_table) <- c("V1", "V2", "V3")

#################

all_probes <- unique(c(top9$probe, top15$probe, top22$probe))
venn_df <- data.frame(probe = all_probes,
                      Y9 = ifelse(all_probes %in% top9$probe, 1, 0),
                      Y15 = ifelse(all_probes %in% top15$probe, 1, 0),
                      Y22 = ifelse(all_probes %in% top22$probe, 1, 0))
venn_data <- list(
  "Age 9" = which(venn_df$Y9 == 1),
  "Age 15" = which(venn_df$Y15 == 1),
  "Age 22" = which(venn_df$Y22 == 1)
)
}

############################
all<-eulerprep_fun(top9,top15,top22)
library(Cairo)
library(ggVennDiagram)
library(eulerr)
CairoPDF("path to save Fig3_y9m_y15m_y22e.all.cpg_overlap.pdf",family="Arial",width=8,height=8)
plot(venn(all), quantities = list(cex=1.5), labels = list(cex=1.75), fills = c("white","#CCC","orange"))
dev.off()

###########################
simpletable_fun<-function(top9,top15,top22){
    # How many overlap across all three?
    overlap_915 <- intersect(top9$probe, top15$probe)
    overlap_922 <- intersect(top9$probe, top22$probe)
    overlap_1522 <- intersect(top15$probe, top22$probe)
    overlap_all <- intersect(overlap_915, top22$probe) 
    }# 53 probes

alloverlap<-simpletable_fun(top9,top15,top22)    
y9alloverlap<-y9[y9$probe %in% alloverlap,]
write.table(y9alloverlap[,c(6,1:9,37:53)],file="persistent_DMPs_annotation.txt",sep="\t",quote=F,row.names=F)



