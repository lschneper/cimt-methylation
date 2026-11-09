TXDB="path to gencode.v37lift37.annotation.txdb"
MANIFEST_FILE="path to EPIC.hg19.manifest.original.tsv.gz.tsv"
OUTDIR<-"output directory path"

.libPaths("path to R package library")
library(GenomicRanges)
library(GenomeInfoDb)


library(AnnotationDbi)
library(org.Hs.eg.db)
library(ChIPseeker)


txdb<-loadDb(TXDB)

args<-commandArgs(trailingOnly=TRUE)

inputfile<-"path to EWAS output RDS file"

man<-read.delim("path to EPIC.hg19.manifest.original.tsv.gz.tsv")
man<-man[,c("Name","CHR","MAPINFO","Relation_to_UCSC_CpG_Island","Regulatory_Feature_Group")]
man$Chr<-paste0("chr",man$CHR)
man$Start<-man$MAPINFO
man$End<-man$MAPINFO+1
man<-man[,-2]
man$Relation_to_UCSC_CpG_Island[man$Relation_to_UCSC_CpG_Island==""]<-"Open_Sea"
names(man)[names(man) %in% "Relation_to_UCSC_CpG_Island"]<-"Relation_to_Island"
filepre<-gsub(".rds$|.RDS$|.txt$","",inputfile)


if (grepl("\\.rds$", inputfile, ignore.case = TRUE)) {
  meta <-readRDS(inputfile)
} else if (grepl("\\.txt$",inputfile, ignore.case=TRUE)){
  meta<-read.delim(inputfile,header=T)
}

if("ID" %in% names(meta)){print("ID column exists")} else {meta$ID<-row.names(meta)}

if("distanceToTSS" %in% names(meta) & "Relation_to_Island" %in% names(meta)) {meta2<-meta} else {meta2<-merge(man,meta,by.x="Name",by.y="ID",all.x=F,all.y=T)}
meta<-meta2

if("Chr" %in% names(meta)& "seqnames"%in% names(meta)){meta<-meta[,!names(meta) %in% "Chr"]}else {print("meta only has one chr column")}
if("Start" %in% names(meta)& "start"%in% names(meta)){meta<-meta[,!names(meta) %in% "Start"]}else {print("meta only has one start column")}
if("End" %in% names(meta)& "end"%in% names(meta)){meta<-meta[,!names(meta) %in% "End"]} else {print("meta only has one end column")}


seqnamespriority<-c("seqnames","Chr","chrom","chr")
seqnamescol<-names(meta)[names(meta)%in% seqnamespriority][1]
startpriority<-c("Start","start","MAPINFO")
startcol<-names(meta)[names(meta)%in% startpriority][1]

# Create GRanges object from input file
gr <- makeGRangesFromDataFrame(meta,
  seqnames.field =  seqnamescol,
  start.field = startcol, end.field = "end",
  keep.extra.columns=TRUE)

# ---- Coordinate matching ----
#seqlevelsStyle(gr) <- seqlevelsStyle(txdb)

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

## Annotate relation to CpG islands
# load manifest annotation
#anno_man <- as.data.frame(getAnnotation(manifest_annotation))
#anno_man<-anno_man[anno_man$Name %in% anno_df$probe,c("Name","Relation_to_Island")]
#anno_df2<-merge(anno_df,anno_man,by.x="probe",by.y="Name")
saveRDS(anno_df,file=paste0(OUTDIR,basename(filepre),"_chipseeker.RDS"))

filepre<-paste0(OUTDIR,basename(filepre))
filepre
##################################
anno_df$Relation_to_Island[anno_df$Relation_to_Island=="OpenSea"]<-"Open_Sea"  
# ---- Keep unfiltered copy ----

anno_df_all <- anno_df

# ---- Filter significant DMRs ----
# change from bacon.fdr_random to fdr_random for cross-sectional analysis

priority<-c("bacon.fdr_fixed","fdr_fixed","bacon.fdr","fdr")
fdr_col <- priority[priority %in% names(anno_df_all)][1]
if(is.na(fdr_col)) {
  stop("No FDR column found.")
}

anno_df<-anno_df_all[anno_df_all[[fdr_col]]<0.05,]

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
# allsum<-colSums(tfeatures)[2]
# sigsum<-colSums(tfeatures)[1]
# tfeatures$hitsminuscat<-sigsum-tfeatures$hits
# tfeatures$allminuscat<-allsum-tfeatures$all

########################

fishertest_fun<-function(df){
  allsum<-colSums(df)[2]
  sigsum<-colSums(df)[1]
  df$hitsminuscat<-sigsum-df$hits
  df$allminuscat<-allsum-df$all
  # for fisher test
  df$a<-df$hits
  df$b<-df$hitsminuscat
  df$c<-df$all-df$a
  df$d<-df$allminuscat-df$b
  
  countsmatrixlist<-list()
  
  fishertestdf<-as.data.frame(matrix(ncol=5,nrow=nrow(df)))
  names(fishertestdf)<-c("element","pval","CI_low","CI_high","OR")
  fishertestdf$element<-row.names(df)
  for(i in 1:nrow(df)){
    countsmatrixlist[[i]]<-matrix(
      c(df$a[i],df$b[i],df$c[i],df$d[i]),
      nrow=2,
      byrow=TRUE,
      dimnames=list(
        c("Significant","Nonsignificant"),
        c("Category","Other")
      )
    )
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

########################
promoterfishertest<-fishertest_fun(tfeatures)
promoterfishertest$Row.names<-factor(promoterfishertest$Row.names,
                               levels=c("promoter_3kb","promoter_2kb","promoter_1kb","5UTR","exon","intron","3UTR","distal_intergenic"),
                               labels=c("Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5' UTR","Exon","Intron","3' UTR","Distal Intergenic"))

write.table(promoterfishertest,file=paste0(filepre,"_promoter_fishertest.table.txt"),row.names=F,quote=F,sep="\t")

###########################

categories<-c("Island","Open_Sea","N_Shelf","S_Shelf","N_Shore","S_Shore")
relation_to_island_df<-rbind(hits=table(factor(hits$Relation_to_Island,levels=categories)),all=table(factor(anno_df_all$Relation_to_Island,levels=categories)))

trelationtoisland<-as.data.frame(t(relation_to_island_df))
if(ncol(trelationtoisland)<2){
  cpg<-(paste("no significant hits in ",inputfile))} else {
allsum<-sum(trelationtoisland$all)
sigsum<-sum(trelationtoisland$hits)
trelationtoisland$hitsminuscat<-sigsum-trelationtoisland$hits
trelationtoisland$allminuscat<-allsum-trelationtoisland$all


cpg<-fishertest_fun(trelationtoisland)
cpg$Row.names<-factor(cpg$Row.names,
                               levels=c("N_Shelf","N_Shore","Island","S_Shore","S_Shelf","Open_Sea"),
                               labels=c("North Shelf","North Shore","CpG Island","South Shore","South Shelf","Open Sea"))
}
write.table(cpg,file=paste0(filepre,"_relationtoCpGIsland_fishertest.table.txt"),row.names=F,quote=F,sep="\t")

#anno_df_all2<-merge(anno_df_all,annotation_hits,by.x="probe",by.y="Name")
#saveRDS(anno_df_all2,file="y15_meta_age_sex_7svs_newmiddleton_techpcs_031426_1_500000_metalistdf2.annotated.RDS")

quit(save="no")
