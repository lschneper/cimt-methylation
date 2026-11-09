library(eulerr)
library(ggVennDiagram)
library(ggplot2)
library(Cairo)
DIR<-"path to working directory"
setwd(DIR)

# Read in file with overlap information
all_multiinter<-read.delim("dmrs_from_combp_dmrcate_for_venn.txt")

all_multiinter$dmr<-paste(all_multiinter$seqnames,all_multiinter$start,all_multiinter$end,sep="_")
row.names(all_multiinter)<-all_multiinter$dmr
# Venn diagrams

all_multiinter$probe<-row.names(all_multiinter)
top9dmr<-all_multiinter[all_multiinter$y9==TRUE,]
top15dmr<-all_multiinter[all_multiinter$y15==TRUE,]


top22dmr<-all_multiinter[all_multiinter$y22==TRUE,]


eulerprep_fun<-function(top9,top15,top22){
  # How many overlap across all three?
  overlap_915 <- intersect(top9$probe, top15$probe)
  overlap_922 <- intersect(top9$probe, top22$probe)
  overlap_1522 <- intersect(top15$probe, top22$probe)
  overlap_all <- intersect(overlap_915, top22$probe) # 42 probes

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

all<-eulerprep_fun(top9dmr,top15dmr,top22dmr)
#CairoPDF(file=paste0(DIR,"/Fig3_DMRvenn_091726.pdf"),family="Arial",height=8,width=8)
CairoPDF(file=paste0(DIR,"/Fig3_DMRvenn_dmrcate_combp_092626.pdf"),family="Arial",height=8,width=8)
plot(venn(all), quantities = list(cex=1.5), labels = list(cex=1.75),fills = c("white",
                                                      "#CCC",
                                                      "lightblue"))
dev.off()


#check<-read.delim("dmrs_from_2methods_with_coordinates.txt")
#table(check$category)

capture.output(
  print(sessionInfo()),
  file = paste0(DIR, "/Fig3_LS_dmr_overlap_091726.R_sessionInfo.txt"))
