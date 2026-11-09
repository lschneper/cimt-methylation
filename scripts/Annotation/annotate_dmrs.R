.libPaths("path to R package library")
library(txdbmaker)

library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomeInfoDb)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(ggplot2)
library(ChIPseeker)

inputfile="path to DMR file"

ewasfile<-"path to corresponding EWAS file"
filepre<-gsub(".rds$|.RDS$|.txt$|.bed","",inputfile)

manifest_annotation<-ifelse(grepl("450k",inputfile,ignore.case = TRUE),"IlluminaHumanMethylation450kanno.ilmn12.hg19","IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

if (grepl("\\.rds$", inputfile, ignore.case = TRUE)) {
  meta <-readRDS(inputfile)
} else if (grepl("\\.txt$",inputfile, ignore.case=TRUE)){
  meta<-read.delim(inputfile,header=T)
} else {meta<-read.delim(inputfile,header=F)}

if (grepl("\\.rds$", ewasfile, ignore.case = TRUE)) {
  ewas <-readRDS(ewasfile)
} else if (grepl("\\.txt$",ewasfile, ignore.case=TRUE)){
  ewas<-read.delim(ewasfile,header=T)
}
if(!"probe" %in% colnames(ewas)){
  ewas$probe <-ewas$Name
}
# Load Gencode v37 annotation
txdb <- loadDb(
  "path to gencode.v37lift37.annotation.txdb"
)

gr <- makeGRangesFromDataFrame(meta,
  seqnames.field =  "V1",
  start.field = "V2", end.field = "V3",
  keep.extra.columns=TRUE)


# ---- Coordinate matching ----
seqlevelsStyle(gr) <- seqlevelsStyle(txdb)

# ---- 1. annotatePeak (nearest + flanking) ----
peak_anno <- annotatePeak(
  gr,
  TxDb = txdb,
  level = "gene",
  tssRegion = c(-3000, 3000),
  sameStrand = FALSE,
  addFlankGeneInfo = TRUE
)

anno_df <- as.data.frame(peak_anno)

# ---- Clean and map gene IDs ----
anno_df$clean_geneId <- gsub("\\..*", "", anno_df$geneId)
anno_df$geneSymbol <- mapIds(
  org.Hs.eg.db,
  keys = anno_df$clean_geneId,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
anno_df$ENTREZID <- mapIds(
  org.Hs.eg.db,
  keys = anno_df$clean_geneId,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
anno_df$geneSymbol <- ifelse(!is.na(anno_df$geneSymbol),
                             anno_df$geneSymbol,
                             anno_df$clean_geneId)

# ---- 2. findOverlaps: all genes that physically overlap each DMR ----
gene_ranges <- genes(txdb)
hits <- findOverlaps(gr, gene_ranges)

dmr2genes <- split(names(gene_ranges)[subjectHits(hits)], queryHits(hits))

overlapping_geneIDs <-  rep(NA_character_,length(gr))
idx <- as.integer(names(dmr2genes))
overlapping_geneIDs[idx] <-vapply(
  dmr2genes,
  function(x) paste(unique(x),collapse=","),
  character(1)
)
anno_df$overlapping_geneIDs<-overlapping_geneIDs

all_ids<-unique(unlist(strsplit(
  na.omit(anno_df$overlapping_geneIDs),","
)))

all_ids<-gsub("\\..*","",all_ids)
gene_symbol_map<-mapIds(
  org.Hs.eg.db,
  keys=all_ids,
  column="SYMBOL",
  keytype = "ENSEMBL",
  multiVals="first"
)

anno_df$overlapping_geneSymbols<-vapply(
  anno_df$overlapping_geneIDs,
  function(ids){
    if(is.na(ids))return(NA_character_)
    ens_ids<-gsub("\\..*","",unlist(strsplit(ids,",")))
    syms<-gene_symbol_map[ens_ids]
    paste(unique(na.omit(syms)),collapse=",")
  },
  character(1)
)

# Annotate promoter overlaps - use 2kb up, 0.5 kb down
promoters_gr<-promoters(genes(txdb),upstream=2000, downstream=500)

# Overlap DMRs with promoters
prom_hits<-findOverlaps(gr, promoters_gr)
# split overlaps by DMR
dmr2prom <- split(
  mcols(promoters_gr)$gene_id[subjectHits(prom_hits)],
  queryHits(prom_hits)
)
# Build promoter gene list per DMR
overlapping_promoterIDs <-  rep(NA_character_,length(gr))
idx <- as.integer(names(dmr2prom))
overlapping_promoterIDs[idx] <-vapply(
  dmr2prom,
  function(x) paste(unique(x),collapse=","),
  character(1)
)
anno_df$overlapping_promoterIDs<-overlapping_promoterIDs
# Map promoter genes to symbols
all_ids<-unique(unlist(strsplit(
  na.omit(anno_df$overlapping_promoterIDs),","
)))

all_ids<-gsub("\\..*","",all_ids)
gene_symbol_map<-mapIds(
  org.Hs.eg.db,
  keys=all_ids,
  column="SYMBOL",
  keytype = "ENSEMBL",
  multiVals="first"
)
# Add promoter symbols to annotation table
anno_df$overlapping_promoterSymbols<-vapply(
  anno_df$overlapping_promoterIDs,
  function(ids){
    if(is.na(ids))return(NA_character_)
    ens_ids<-gsub("\\..*","",unlist(strsplit(ids,",")))
    syms<-gene_symbol_map[ens_ids]
    paste(unique(na.omit(syms)),collapse=",")
  },
  character(1)
)

# Annotate relation to CpG islands
anno_df$dmr_id<-seq_len(nrow(anno_df))
gr <- makeGRangesFromDataFrame(anno_df,
                                  seqnames.field =  "seqnames",
                                  start.field = "start", end.field = "end",
                                  keep.extra.columns=TRUE)
# load manifest annotation
anno_epic <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
# Subset manifest by ewas data

anno_epic<-anno_epic[anno_epic$Name %in% ewas$probe,]

anno_epic<-anno_epic[,c("chr","pos","Relation_to_Island")]
anno_epic$end<-anno_epic$pos + 1
anno_epic<-anno_epic[,c("chr","pos","end","Relation_to_Island")]
anno_gr<- makeGRangesFromDataFrame(anno_epic,
                                   seqnames.field =  "chr",
                                   start.field = "pos", end.field = "end",
                                   keep.extra.columns=TRUE)

anno_dmr_overlaps<- findOverlaps(anno_gr,gr)
qh<-queryHits(anno_dmr_overlaps)
sh<-subjectHits(anno_dmr_overlaps)
Relation_to_Island<-as.character(mcols(anno_gr)$Relation_to_Island)
anno_dmr_overlapsdf<-data.frame(dmr=sh,
                                    Relation_to_Island=Relation_to_Island[qh])

# Collapse North and South for shores and shelves to prevent counting ties twice
anno_dmr_overlapsdf$Relation_to_Island<-ifelse(
  anno_dmr_overlapsdf$Relation_to_Island=="Island","Island", ifelse(
    anno_dmr_overlapsdf$Relation_to_Island=="OpenSea","OpenSea", ifelse(
      anno_dmr_overlapsdf$Relation_to_Island %in% c("N_Shelf","S_Shelf"),"Shelf", ifelse(
      anno_dmr_overlapsdf$Relation_to_Island %in% c("N_Shore","S_Shore"),"Shore", NA))))

dmr_props_all<-anno_dmr_overlapsdf %>%
  group_by(dmr,Relation_to_Island) %>%
  summarise(n=n(),.groups="drop") %>%
  group_by(dmr) %>%
  mutate(prop=n/sum(n))

dominant<-dmr_props_all %>%
  group_by(dmr) %>%
  slice_max(n,n=1, with_ties=FALSE)


anno_df2<-merge(anno_df,dominant,by.x="dmr_id",by.y="dmr")

# Save data
saveRDS(anno_df2,file=paste0(filepre,"_chipseeker.RDS"))
print(paste0("Finished processing annotate_dmr for ",filepre))

# 
# ---- Keep unfiltered copy ----
anno_df_all <- anno_df2


#Test enrichment of features
# ---- Filter significant DMRs ----
anno_df <- anno_df_all[anno_df_all$V7< 0.05 & anno_df_all$V5>=5, ]

# Classify hits
hits<-anno_df
hitdf<-data.frame(promoter_1kb=nrow(hits[hits$annotation == "Promoter (<=1kb)",]),
                  promoter_2kb=nrow(hits[hits$annotation=="Promoter (1-2kb)",]),
                  promoter_3kb=nrow(hits[hits$annotation=="Promoter (2-3kb)",]),
                  intron=nrow(hits[grepl("Intron", hits$annotation), ]),
                  exon=nrow(hits[grepl("Exon", hits$annotation), ]),
                  `5UTR`=nrow(hits[hits$annotation == "5' UTR",]),
                  `3UTR`=nrow(hits[hits$annotation == "3' UTR",]),
                  distal_intergenic=nrow(hits[hits$annotation == "Distal Intergenic",]),
                  check.names = F
)

# Classify all
alldf<-data.frame(promoter_1kb=nrow(anno_df_all[anno_df_all$annotation == "Promoter (<=1kb)",]),
                  promoter_2kb=nrow(anno_df_all[anno_df_all$annotation=="Promoter (1-2kb)",]),
                  promoter_3kb=nrow(anno_df_all[anno_df_all$annotation=="Promoter (2-3kb)",]),
                  intron=nrow(anno_df_all[grepl("Intron", anno_df_all$annotation), ]),
                  exon=nrow(anno_df_all[grepl("Exon", anno_df_all$annotation), ]),
                  `5UTR`=nrow(anno_df_all[anno_df_all$annotation == "5' UTR",]),
                  `3UTR`=nrow(anno_df_all[anno_df_all$annotation == "3' UTR",]),
                  distal_intergenic=nrow(anno_df_all[anno_df_all$annotation == "Distal Intergenic",]),check.names=F
)
features<-rbind(hitdf,alldf)
row.names(features)<-c("hits","all")
tfeatures<-as.data.frame(t(features))
#####
fishertest_fun<-function(df){
  allsum<-colSums(df)[2]
  sigsum<-colSums(df)[1]
  df$hitsminuscat<-sigsum-df$hits
  df$allminuscat<-allsum-df$all
  
  countsmatrixlist<-list()
  fishertestdf<-as.data.frame(matrix(ncol=5,nrow=nrow(df)))
  names(fishertestdf)<-c("element","pval","CI_low","CI_high","OR")
  fishertestdf$element<-row.names(df)
  for(i in 1:nrow(df)){
    countsmatrixlist[[i]]<-matrix(c(df[i,"hits"],df[i,"hitsminuscat"],df[i,"all"],df[i,"allminuscat"]),
                                  nrow=2, byrow = 2)
    ft<-fisher.test(countsmatrixlist[[i]],alternative = "two.sided")
    fishertestdf[i,"pval"]<-ft$p.value
    fishertestdf[i,"OR"]<-ft$estimate
    fishertestdf[i,"CI_low"]<-ft$conf.int[1]
    fishertestdf[i,"CI_high"]<-ft$conf.int[2]
  }
  fishertestdf$fdr<-p.adjust(fishertestdf$pval,method="fdr")
  fishertestdf$sig<-ifelse(fishertestdf$fdr<0.001,"***",ifelse(fishertestdf$fdr<0.01,"**",ifelse(fishertestdf$fdr<0.05,"*","NS")))
  fishertestdf<-merge(df,fishertestdf,by.x="row.names",by.y="element")
  #
  fishertestdf$sig[fishertestdf$sig=="NS"]<-NA
  fishertestdf$label_y_pos<-fishertestdf$OR + 0.1
  return(fishertestdf)
}
fishertestdf<-fishertest_fun(tfeatures)
#####
fishertestdf$Row.names<-factor(fishertestdf$Row.names,
                               levels=c("promoter_3kb","promoter_2kb","promoter_1kb","5UTR","exon","intron","3UTR","distal_intergenic"),
                               labels=c("Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5' UTR","Exon","Intron","3' UTR","Distal Intergenic"))

write.table(fishertestdf,file=paste0(filepre,"_promoter_fishertest.table.txt"),row.names=F,quote=F,sep="\t")


##### Relation to CpG Island Fisher test

relation_to_island_df<-rbind(hits=table(anno_df$Relation_to_Island),all=table(anno_df_all$Relation_to_Island))

trelationtoisland<-as.data.frame(t(relation_to_island_df))
if(ncol(trelationtoisland)<2){
  cpg<-(paste("no significant hits in ",inputfile))} else {
    allsum<-sum(trelationtoisland$all)
    sigsum<-sum(trelationtoisland$hits)
    trelationtoisland$hitsminuscat<-sigsum-trelationtoisland$hits
    trelationtoisland$allminuscat<-allsum-trelationtoisland$all
    
    
    cpg<-fishertest_fun(trelationtoisland)
    cpg$Row.names<-factor(cpg$Row.names,
                          levels=c("Island","Shore","Shelf","OpenSea"),
                          labels=c("CpG Island","Shore","Shelf","Open Sea"))
  }


write.table(cpg,file=paste0(filepre,"_relationtoCpGIsland_fishertest.table.txt"),row.names=F,quote=F,sep="\t")
