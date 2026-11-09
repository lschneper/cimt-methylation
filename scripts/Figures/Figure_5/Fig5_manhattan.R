#Authors: Eve Schoeffler and Lisa Schneper

# Load libraries
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(optparse)
library(showtext)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(Cairo)

font_add("Arial", 
         regular = "~/.fonts/arial.ttf",
         bold = "~/.fonts/ARIALBD.TTF")
showtext_auto()
showtext_opts(dpi = 300)



ewasfile<-"/path to y9y15longitudinal_ewas_manhattaninput_092526.RDS"
outpre<-"output_fileprefix"

OUTDIR="output directory path"
outdir=OUTDIR

  #####
# Read in files

if (grepl("\\.rds$", ewasfile, ignore.case = TRUE)) {
  ewas_res <- readRDS(ewasfile)
} else if (grepl("\\.txt$",ewasfile, ignore.case=TRUE)){
  ewas_res <- read.delim(ewasfile,header=T)
}

####
# Extract columns of interest
{
# Chr column
chr_priorities <-c("Chr","seqnames","SeqNames","Chr","chr")
chr_col<-chr_priorities[chr_priorities %in% names(ewas_res)][1]

# # Start column
start_priorities <-c("i.i.start","i.start","start","Start","Pos","pos")
start_col<-start_priorities[start_priorities %in% names(ewas_res)][1]
# 
# # End column
end_priorities <-c("i.i.end","i.end","end","End","END")
end_col<-end_priorities[end_priorities %in% names(ewas_res)][1]

# Pval column
pval_priorities<-c("bacon.P.Value_fixed","bacon.P.Value")
pval_col <- pval_priorities[pval_priorities %in% names(ewas_res)][1]
if(is.na(pval_col)) {
  stop("No PVal column found.")
}

# -log10 Pval column
logpval_priorities<-c("logP")
logpval_col <- logpval_priorities[logpval_priorities %in% names(ewas_res)][1]
if(is.na(logpval_col)) {
  stop("No log10PVal column found.")
}

# FDR column
fdr_priorities<-c("bacon.fdr_fixed","bacon.fdr")
fdr_col<-fdr_priorities[fdr_priorities %in% names(ewas_res)][1]

snp_col_priorities<-c("probe","Name")
snp_col<-snp_col_priorities[snp_col_priorities %in% names(ewas_res)][1]
}

wave_DMR_col_priorities<-c("wave_test_DMR","wave_DMR")
wave_DMR_col<-wave_DMR_col_priorities[wave_DMR_col_priorities %in% names(ewas_res)][1]

common_DMR_col_priorities<-c("common_test_DMR","common_DMR")
common_DMR_col<-common_DMR_col_priorities[common_DMR_col_priorities %in% names(ewas_res)][1]


ewas_res_bed<-ewas_res[,c(chr_col,start_col,end_col,logpval_col,pval_col,fdr_col,snp_col,"chr_len","tot","ewas_hit",wave_DMR_col,common_DMR_col,"sigtype","prom_min")]

# Add ewas_sig column
ewas_res_bed$commonCpG<-ifelse(ewas_res_bed$sigtype=="ewas_hit",1,0)


ewas_res_bed$waveDMR_cpg<-ifelse(ewas_res_bed$sigtype=="waveDMR",1,0)

ewas_res_bed$commonDMR_cpg<-ifelse(ewas_res_bed$sigtype=="commonDMR",1,0)
ewas_res_bed$priority<-ifelse(ewas_res_bed$commonDMR_cpg==1,"commonDMR",ifelse(ewas_res_bed$waveDMR_cpg==1,"waveDMR",ifelse(ewas_res_bed$commonCpG==1,"commonCpG",ifelse(ewas_res_bed$commonCpG==1,"waveCpG","NS"))))

axis_df<-ewas_res_bed %>%
  group_by(Chr) %>%
  summarise(center=mean(tot))
ewas_res_bed2 <- ewas_res_bed %>%
  mutate(chr_group=Chr %% 2,
         logP=logP
  )

ewas_res_bed2$priority2<-ifelse(ewas_res_bed2$priority %in% c("NS"),NA,ewas_res_bed2$priority)

fdr_cutoff<-max(ewas_res_bed2[[pval_col]][ewas_res_bed2[[fdr_col]]<=0.05],na.rm=T)
log10fdrcutoff<- -log10(fdr_cutoff)


CairoPDF(paste0(outdir,outpre,"_manhattan_resized.pdf"),width=12,height=8,family = "Arial")
p<-ggplot(ewas_res_bed2,aes(x=tot, y=logP)) +
  geom_point(aes(color=factor(chr_group)),size=0.8,alpha=0.6,
             color=ifelse(ewas_res_bed2$chr_group==0,"grey75","grey85")) +
 geom_hline(yintercept = log10fdrcutoff,color="blue",alpha=0.4) +
  
  geom_point(
    data=ewas_res_bed2[!is.na(ewas_res_bed2$priority2),],
    aes(x=tot,y=logP,color=priority2),
    size=0.8)+
  scale_color_manual(values=c(
    "commonDMR"="#DF536B",
    "commonCpG" ="#E69F00",
    "waveDMR" = "#2297e6"
  )) +
  scale_y_continuous(breaks=seq(0,ceiling(max(ewas_res_bed2$logP)),by=2)) +

  scale_x_continuous(labels=axis_df$Chr,breaks=axis_df$center) +
  geom_text_repel(
    data=ewas_res_bed2[!is.na(ewas_res_bed2$prom_min)& !is.na(ewas_res_bed2$priority2)&ewas_res_bed2$logP>=log10fdrcutoff,],
    aes(x=tot,y=logP,label=prom_min),
    size=5,
    family="Arial",
    fontface="bold",
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color="black") +
  theme_bw() +
  theme(legend.position = "none",axis.text.y=element_text(family="Arial",size=9),axis.text.x=element_text(family="Arial",size=9,angle=90,hjust=1, vjust=0.5),
        axis.title=element_text(size=12,family="Arial"),panel.grid=element_blank(),
        plot.title = element_text(family="Arial",hjust=0.5,size=14,face="bold"))+xlab("Chromosome")+ylab(expression(-log[10](italic(p)-value)))  + ggtitle(paste0("Age ",gsub("y","",outpre)))

p
dev.off()
