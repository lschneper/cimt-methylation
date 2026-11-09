library(magrittr)
library(dplyr)
setwd("path to working directory")
OUTDIR="path to output directory"
# Need to use all for Venn Diagrams
all_multiinter<-read.table(paste0(OUTDIR,"dmrs_from_dmrcate_combp_with_coordinates.txt"),sep="\t",header=T,stringsAsFactors=F)

y9y15dmrs<-readRDS(paste0(OUTDIR,"y9y15_DMRs_DMRcate_combp_chipseeker.RDS"))
y9dmr<-readRDS("path to y9_DMRs_DMRcate_combp_chipseeker file (adjust mode to read in as necessary)")
y15dmr<-readRDS("path to y15_DMRs_DMRcate_combp_chipseeker (adjust mode to read in as necessary)")
y22dmr<-readRDS("path to y22_DMRs_DMRcate_combp_chipseeker (adjust mode to read in as necessary)")
all_multiinter$dmr<-paste(all_multiinter$seqnames,all_multiinter$start,all_multiinter$end,sep="_")

persistent_multiinter<-all_multiinter[all_multiinter$y9_present==TRUE & all_multiinter$y15_present==TRUE & all_multiinter$y22_present==TRUE,]
y9ewas<-readRDS("path to y9_meta_sensitivity_model1_ewas_output_chipseeker.RDS")
y15ewas<-readRDS("path to y15_meta_sensitivity_model1_ewas_output_chipseeker.RDS")
y22ewas<-readRDS("path to Age 22 epic_sensitivity_model1_ewas_output_chipseeker.RDS")
y22ewas$MAPINFO=y22ewas$start


y9y15y22sig<-read.table(file="path to persistent_ewas_hits.txt",sep="\t",stringsAsFactors=F)
names(y9y15y22sig) <- y9y15y22sig[1,]
y9y15y22sig<-y9y15y22sig[2:nrow(y9y15y22sig),]
y9y15ewas<-readRDS("path to y9_y15_mixed effect model ewas meta-analysis fixed_effects results file")
add_dmrID_fun<-function(df){
  
  df$dmr<-paste(as.character(df$seqnames),df$start,df$end,sep="_")
  return(df)
}

y9y15dmr<-add_dmrID_fun(y9y15dmrs)


dmr_annotation_fun<-function(ewasdf,dmrdf,persistent){
  library(data.table)
  library(dplyr)
  ewasdf<-data.table::copy(ewasdf)
  dmrdf<-data.table::copy(dmrdf)
  persistent<-data.table::copy(persistent)
  multiinter<-persistent[,c("seqnames","start","end","dmr")]
  multiinter$seqnames<-as.character(multiinter$seqnames)
  names(multiinter)<-c("seqnames","start","end","common_test_DMR")
  
  ewasdf$seqnames<-as.character(ewasdf$seqnames)
  dmrdf$seqnames<-as.character(dmrdf$seqnames)
  
  
  setDT(ewasdf)
  setDT(multiinter)
  setkey(multiinter,seqnames,start,end)
  res<-foverlaps(
    ewasdf, multiinter,by.x = c("seqnames","start","end"),by.y = c("seqnames","start","end"),type="within",nomatch=NA)
  dmrpriorities<-c("dmr","dmr_id")
  dmrcol<-names(dmrdf)[names(dmrdf) %in% dmrpriorities][1]
  dmrdf2<-dmrdf[,c("seqnames","start","end",dmrcol)]
  names(dmrdf2)<-c("seqnames","i.start","i.end","wave_test_DMR")
  setDT(dmrdf2)
  setkey(dmrdf2,seqnames,i.start,i.end)
  res2<-foverlaps(
   res,dmrdf2, by.x=c("seqnames","i.start","i.end"),by.y=c("seqnames","i.start","i.end"),type="within",nomatch=NA)
  res2<-as.data.frame(res2)
  fdrpriority<-c("bacon.fdr_fixed","bacon.fdr")
  fdrcol<-names(res2)[names(res2)%in%fdrpriority][1]
  pvalpriority<-c("bacon.P.Value_fixed","bacon.P.Value")
  pvalcol<-names(res2)[names(res2)%in%pvalpriority][1]
  probepriority<-c("probe","Name","ID")
  probecol<-names(res2)[names(res2)%in%probepriority][1]
  
  res2$ewas_hit<-ifelse(res2[,probecol]%in%y9y15y22sig$probe,1,NA)
  res2$sigtype<-NA
  res2$sigtype[!is.na(res2$common_test_DMR)] <- "commonDMR"
  res2$sigtype[is.na(res2$sigtype) & !is.na(res2$wave_test_DMR)] <- "waveDMR"
  res2$sigtype[is.na(res2$sigtype) & !is.na(res2$ewas_hit)] <- "ewas_hit"
  
  res3<-res2 %>%
    group_by(wave_test_DMR) %>%
    mutate(prom_min=if_else(.data[[pvalcol]]==min(.data[[pvalcol]]),overlapping_promoterSymbols,NA_character_)) %>%
    ungroup()
  
  res3<-res3 %>%
    mutate(Chr=as.numeric(gsub("chr","",seqnames))) %>%
    mutate(MAPINFO=as.numeric(res3$i.i.start)) %>%
    arrange(Chr,MAPINFO) %>%
    group_by(Chr) %>%
    mutate(chr_len=max(MAPINFO)) %>%
    ungroup() %>%
    mutate(tot=as.numeric(MAPINFO) + cumsum(as.numeric(lag(chr_len,default=0))))
  
  res3 <- res3 %>%
    mutate(chr_group=Chr %% 2,
           logP=-log10(.data[[pvalcol]])
           )
  return(res3)
  }


y9ewas2<-dmr_annotation_fun(ewasdf=y9ewas,dmrdf=y9dmr,persistent = persistent_multiinter)
y15ewas2<-dmr_annotation_fun(ewasdf=y15ewas,dmrdf=y15dmr,persistent = persistent_multiinter)
y22ewas2<-dmr_annotation_fun(ewasdf=y22ewas,dmrdf=y22dmr,persistent = persistent_multiinter)

y9y15ewasdf2=y9y15ewas;dmrdf=y9y15dmr; persistent=persistent_multiinter
y9y15<-dmr_annotation_fun(ewasdf=y9y15ewas,dmrdf=y9y15dmr,persistent=persistent_multiinter)

saveRDS(y9ewas2,file=paste0(OUTDIR,"/y9_ewas_manhattaninput_dmrcatecombp_092626.RDS"))
saveRDS(y15ewas2,file=paste0(OUTDIR,"/y15_ewas_manhattaninput_dmrcatecombp_092626.RDS")) 
saveRDS(y22ewas2,file=paste0(OUTDIR,"/y22_ewas_manhattaninput_dmrcatecombp_092626.RDS"))
saveRDS(y9y15,file=paste0(OUTDIR,"/y9y15longitudinal_ewas_manhattaninput_dmrcatecombp_092626.RDS"))
