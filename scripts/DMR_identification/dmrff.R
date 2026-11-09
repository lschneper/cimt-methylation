.libPaths("path to R package library")
library(dmrff)
library(qs2)
library(bacon)
source("path to R_functions/bacon_fun")
setwd("path to working directory")

dmrff_fun<-function(statsdf,betas,outdir,outpre){
  statsdf$chr<-as.numeric(gsub("chr","",as.character(statsdf$seqnames)))
  # Remove NA in chromosomes and position
  statsdf<-statsdf[!is.na(statsdf$seqnames),]
  statsdf<-statsdf[!is.na(statsdf$start),]
  statsdfo<-statsdf[order(statsdf$chr,statsdf$start),]
  
  #betas <- readRDS(betas)
  betas<-betas[row.names(betas) %in% statsdfo$probe,]
  betaso<-betas[match(statsdfo$probe,row.names(betas)),]
  M<- log2(betaso) - log2(1 - betaso)
  
  rm(betas)
  rm(betaso)
  rm(statsdf)
  gc()
  
  dmrs<-dmrff(estimate = statsdfo$bacon.est,se=statsdfo$bacon.se,p.value = statsdfo$bacon.P.Value,methylation = M,chr=statsdfo$chr,pos=statsdfo$start)
  dmrsfilt<-dmrs[which(dmrs$p.adjust<0.1 & dmrs$n>1),]
  if(nrow(dmrsfilt)>0) {
  sites<-dmrff.sites(dmrsfilt,statsdfo$chr,statsdfo$start)
  sites<-cbind(sites,statsdfo[sites$site,c("overlapping_geneSymbols","overlapping_promoterSymbols","annotation","distanceToTSS","probe","start","end","bacon.fdr","bacon.P.Value","bacon.est")])
  sites<-cbind(sites,dmr=dmrsfilt[sites$region,c("chr","start","end","n","B","estimate","se","z","p.value","number.tests","p.adjust")])
write.table(sites,file=paste0(outdir,outpre,"_sigdmrffoutput.txt"),sep="\t",row.names=F,quote=F)
write.table(dmrs,file=paste0(outdir,outpre,"_alldmrffoutput.txt"),sep="\t",row.names=F,quote=F)


pdf(paste0(outdir,outpre,"_sigdmrffoutput_plots.pdf"))
  for(i in 1:nrow(dmrsfilt)){
    print(dmrff.plot(region.chr = dmrsfilt[i,"chr"],region.start = dmrsfilt[i,"start"],region.end=dmrsfilt[i,"end"],estimate = statsdfo$bacon.est,statsdfo$bacon.se,statsdfo$chr,statsdfo$start,ci=T,expand=1))
  }
dev.off()
} else {
  message("No significant DMRs found - skipping output and plots")
}
dmrpre<-dmrff.pre(estimate = statsdfo$bacon.est,se=statsdfo$bacon.se,methylation = M,chr=statsdfo$chr,pos=statsdfo$start)
return(dmrpre)
}
  
## Note - need to subset betas using samplesheets  
######################
y22<-readRDS("path to Age 22 EWAS output")
beta <- qs_read("path to BetaMatrix.blood.n1271.qcfilter.020425.qs") # y22
survey<-readRDS("path to file containing dependent variable and covariate data with array ID")
beta<-beta[,colnames(beta) %in% survey$Array_ID]

outdir<-"path to output directory"
outpre="output file prefix"
y22dmrff<-dmrff_fun(statsdf = y22,betas=beta,outdir=outdir,outpre=outpre)

gc()
################
# y15 epic
###################
# Read in EPIC manifest
man <- readRDS("/path to EPIC.hg19.manifest.dataframe.rds")
man$ID<-row.names(man)
# # # Exclude bad probes and sex chromosomes
manSub <- man[man$seqnames %in% c("chrX","chrY","chrM") | !substr(man$ID,1,2) %in% "cg","ID"]
man<-man[!man$ID %in% manSub,c("ID","seqnames","start","end")]
man$CHR<-gsub("chr","",as.character(man$seqnames))
man$MAPINFO<-man$start

survey<-readRDS("path to file containing dependent variable and covariate data with array ID")
# Read in beta matrix for y9 & 15
beta <- qs_read("path to BetaMatrix.saliva.n2627.qcfilter.020425.qs") 

# Subset beta matrix by y15 samples
beta<-beta[,colnames(beta) %in% survey$Array_ID]

y15_epic_stats<-read.delim("path to Age 15_epic_sensitivity_model1.ewasoutput.txt")

y15_epic_stats_m<-merge(man,y15_epic_stats,by.x="ID",by.y="probe")
y15_epic_stats_m$CHR<-as.numeric(as.character(y15_epic_stats_m$CHR))

beta<-beta[row.names(beta)%in% y15_epic_stats_m$ID,]

y15_epic_stats<-y15_epic_stats_m[y15_epic_stats_m$ID %in% row.names(beta),]
statsdf<-y15_epic_stats
betas<-beta
outdir ="path to output directory"
outpre="output prefix"

statsdf$chr<-as.numeric(gsub("chr","",as.character(statsdf$seqnames)))
# Remove NA in chromosomes and position
statsdf<-statsdf[!is.na(statsdf$CHR),]
statsdf<-statsdf[!is.na(statsdf$MAPINFO),]
statsdfo<-statsdf[order(statsdf$CHR,statsdf$MAPINFO),]
statsdfo$Coefficient<-as.numeric(statsdfo$Coefficient)
statsdfo$t<-as.numeric(statsdfo$t)
statsdfo$StdError<-as.numeric(statsdfo$StdError)
statsdfo<-bacon_fun(statsdfo)
statsdfo$bacon.fdr<-p.adjust(statsdfo$bacon.P.Value,method="fdr")


#betas <- readRDS(betas)
betas<-betas[row.names(betas) %in% statsdfo$ID,]
betaso<-betas[match(statsdfo$ID,row.names(betas)),]
M<- log2(betaso) - log2(1 - betaso)

rm(betas)
rm(betaso)
rm(statsdf)
gc()
dmrs<-dmrff(estimate = statsdfo$bacon.est,se=statsdfo$bacon.se,p.value = statsdfo$bacon.P.Value,methylation = M,chr=statsdfo$CHR,pos=statsdfo$MAPINFO)
dmrsfilt<-dmrs[which(dmrs$p.adjust<0.1 & dmrs$n>1),]
  sites<-dmrff.sites(dmrsfilt,statsdfo$CHR,statsdfo$MAPINFO)
  #  sites<-cbind(sites,statsdfo[sites$site,c("Gene","probe","start","end","bacon.P.Value","bacon.est")])
  sites<-cbind(sites,statsdfo[sites$site,c("ID","CHR","MAPINFO","bacon.fdr","bacon.P.Value","bacon.est")])
  sites<-cbind(sites,dmr=dmrsfilt[sites$region,c("chr","start","end","n","B","estimate","se","z","p.value","number.tests","p.adjust")])
  write.table(sites,file=paste0(outdir,outpre,"_sigdmrffoutput.txt"),sep="\t",row.names=F,quote=F)
  write.table(dmrs,file=paste0(outdir,outpre,"_alldmrffoutput.txt"),sep="\t",row.names=F,quote=F)
  
  
  pdf(paste0(outdir,outpre,"_sigdmrffoutput_plots.pdf"))
  for(i in 1:nrow(dmrsfilt)){
    print(dmrff.plot(region.chr = dmrsfilt[i,"chr"],region.start = dmrsfilt[i,"start"],region.end=dmrsfilt[i,"end"],estimate = statsdfo$bacon.est,statsdfo$bacon.se,statsdfo$CHR,statsdfo$MAPINFO,ci=T,expand=1))
  }
  dev.off()

y15epic_betas<-M
y15epic_statsdfo<-statsdfo

########
# y9 epic
#########

# y9 epic
# Read in EPIC manifest
man <- readRDS("path to EPIC.hg19.manifest.dataframe.rds")
man$ID<-row.names(man)
# # # Exclude bad probes and sex chromosomes
manSub <- man[man$seqnames %in% c("chrX","chrY","chrM") | !substr(man$ID,1,2) %in% "cg","ID"]
man<-man[!man$ID %in% manSub,c("ID","seqnames","start","end")]
man$CHR<-gsub("chr","",man$seqnames)
man$MAPINFO<-man$start

survey<-readRDS("path to file containing dependent variable and covariate data with array ID")

# Read in beta matrix for y9 & 15
beta <- qs_read("path to BetaMatrix.saliva.n2627.qcfilter.020425.qs") 
# Subset beta matrix by samples
beta<-beta[,colnames(beta)%in% survey$Array_ID]

y9_epic_stats<-read.delim("path to y9_epic_sensitivity_model1.ewasoutput.txt")

y9_epic_stats_m<-merge(man,y9_epic_stats,by.x="ID",by.y="probe")
y9_epic_stats_m$CHR<-as.numeric(as.character(y9_epic_stats_m$CHR))
y9_epic_stats_m$start<-as.numeric(as.character(y9_epic_stats_m$start))
y9_epic_stats_m$MAPINFO<-as.numeric(as.character(y9_epic_stats_m$MAPINFO))
y9_epic_stats_m$t<-as.numeric(as.character(y9_epic_stats_m$t))
y9_epic_stats_m$StdError<-as.numeric(as.character(y9_epic_stats_m$StdError))
y9_epic_stats_m$Coefficient<-as.numeric(as.character(y9_epic_stats_m$Coefficient))
y9_epic_stats_m<-bacon_fun(y9_epic_stats_m)

beta<-beta[row.names(beta)%in% y9_epic_stats_m$ID,]

y9_epic_stats<-y9_epic_stats_m[y9_epic_stats_m$ID %in% row.names(beta),]
statsdf<-y9_epic_stats
betas<-beta
outdir ="output directory file path"
outpre="output prefix"

# Remove NA in chromosomes and position
statsdf<-statsdf[!is.na(statsdf$CHR),]
statsdf<-statsdf[!is.na(statsdf$MAPINFO),]
statsdfo<-statsdf[order(statsdf$CHR,statsdf$MAPINFO),]
statsdfo$bacon.fdr<-p.adjust(statsdfo$bacon.P.Value,method="fdr")
#betas <- readRDS(betas)
betas<-betas[row.names(betas) %in% statsdfo$ID,]
betaso<-betas[match(statsdfo$ID,row.names(betas)),]
M<- log2(betaso) - log2(1 - betaso)
rm(betas)
rm(betaso)
rm(statsdf)
gc()

dmrs<-dmrff(estimate = statsdfo$bacon.est,se=statsdfo$bacon.se,p.value = statsdfo$bacon.P.Value,methylation = M,chr=statsdfo$CHR,pos=statsdfo$MAPINFO)
dmrsfilt<-dmrs[which(dmrs$p.adjust<0.1 & dmrs$n>1),]
  sites<-dmrff.sites(dmrsfilt,statsdfo$CHR,statsdfo$MAPINFO)
  #  sites<-cbind(sites,statsdfo[sites$site,c("Gene","probe","start","end","bacon.P.Value","bacon.est")])
  sites<-cbind(sites,statsdfo[sites$site,c("ID","CHR","MAPINFO","bacon.fdr","bacon.P.Value","bacon.est")])
  sites<-cbind(sites,dmr=dmrsfilt[sites$region,c("chr","start","end","n","B","estimate","se","z","p.value","number.tests","p.adjust")])
  write.table(sites,file=paste0(outdir,outpre,"_sigdmrffoutput.txt"),sep="\t",row.names=F,quote=F)
  write.table(dmrs,file=paste0(outdir,outpre,"_alldmrffoutput.txt"),sep="\t",row.names=F,quote=F)
  
  
  pdf(paste0(outdir,outpre,"_sigdmrffoutput_plots.pdf"))
  for(i in 1:nrow(dmrsfilt)){
    print(dmrff.plot(region.chr = dmrsfilt[i,"chr"],region.start = dmrsfilt[i,"start"],region.end=dmrsfilt[i,"end"],estimate = statsdfo$bacon.est,statsdfo$bacon.se,statsdfo$CHR,statsdfo$MAPINFO,ci=T,expand=1))
  }
  dev.off()

y9epic_betas<-M
y9epic_statsdfo<-statsdfo

####################
# y15 450K
###################
man <- readRDS("/path to hm450.hg19.manifest.dataframe.rds")
man$ID<-row.names(man)
# # # Exclude bad probes and sex chromosomes
manSub <- man[man$seqnames %in% c("chrX","chrY","ID","chrM"),"ID"]
man<-man[!man$ID %in% manSub,c("ID","seqnames","start","end")]
man$CHR<-gsub("chr","",man$seqnames)
man$MAPINFO<-man$start

survey<-readRDS("path to file containing dependent variable and covariate data with array ID")
# Read in stats
y15_450K_stats<-read.delim("path to y15_450K_sensitivity_model1.ewasoutput.txt")

# Read in beta matrix for y9 & 15
beta <- qs_read("path to BetaMatrix.k450.n1724.qcfilter.020425.qs")
beta<-beta[,colnames(beta)%in% survey$Array_ID]

y15_450K_stats_m<-merge(man,y15_450K_stats,by.x="ID",by.y="probe")
y15_450K_stats_m$CHR<-as.numeric(as.character(y15_450K_stats_m$CHR))
y15_450K_stats_m$start<-as.numeric(as.character(y15_450K_stats_m$start))
y15_450K_stats_m$MAPINFO<-as.numeric(as.character(y15_450K_stats_m$MAPINFO))
y15_450K_stats_m$t<-as.numeric(as.character(y15_450K_stats_m$t))
y15_450K_stats_m$StdError<-as.numeric(as.character(y15_450K_stats_m$StdError))
y15_450K_stats_m$Coefficient<-as.numeric(as.character(y15_450K_stats_m$Coefficient))
y15_450K_stats_m<-bacon_fun(y15_450K_stats_m)

beta<-beta[row.names(beta)%in% y15_450K_stats_m$ID,]

y15_450K_stats<-y15_450K_stats_m[y15_450K_stats_m$ID %in% row.names(beta),]
statsdf<-y15_450K_stats
betas<-beta
outdir ="/output directory file path"
outpre="output file prefix"

# Remove NA in chromosomes and position
statsdf<-statsdf[!is.na(statsdf$CHR),]
statsdf<-statsdf[!is.na(statsdf$MAPINFO),]
statsdfo<-statsdf[order(statsdf$CHR,statsdf$MAPINFO),]
statsdfo$bacon.fdr<-p.adjust(statsdfo$bacon.P.Value,method="fdr")

#betas <- readRDS(betas)
betas<-betas[row.names(betas) %in% statsdfo$ID,]
betaso<-betas[match(statsdfo$ID,row.names(betas)),]
M<- log2(betaso) - log2(1 - betaso)
rm(betas)
rm(betaso)
rm(statsdf)
gc()

dmrs<-dmrff(estimate = statsdfo$bacon.est,se=statsdfo$bacon.se,p.value = statsdfo$bacon.P.Value,methylation = M,chr=statsdfo$CHR,pos=statsdfo$MAPINFO)
dmrsfilt<-dmrs[which(dmrs$p.adjust<0.1 & dmrs$n>1),]
if(nrow(dmrsfilt)>0) {
  sites<-dmrff.sites(dmrsfilt,statsdfo$CHR,statsdfo$MAPINFO)
  sites<-cbind(sites,statsdfo[sites$site,c("ID","CHR","MAPINFO","bacon.fdr","bacon.P.Value","bacon.est")])
  sites<-cbind(sites,dmr=dmrsfilt[sites$region,c("chr","start","end","n","B","estimate","se","z","p.value","number.tests","p.adjust")])
  write.table(sites,file=paste0(outdir,outpre,"_sigdmrffoutput.txt"),sep="\t",row.names=F,quote=F)
  write.table(dmrs,file=paste0(outdir,outpre,"_alldmrffoutput.txt"),sep="\t",row.names=F,quote=F)
  
  
  pdf(paste0(outdir,outpre,"_sigdmrffoutput_plots.pdf"))
  for(i in 1:nrow(dmrsfilt)){
    print(dmrff.plot(region.chr = dmrsfilt[i,"chr"],region.start = dmrsfilt[i,"start"],region.end=dmrsfilt[i,"end"],estimate = statsdfo$bacon.est,statsdfo$bacon.se,statsdfo$CHR,statsdfo$MAPINFO,ci=T,expand=1))
  }
  dev.off()
} else {
  message("No significant DMRs found - skipping output and plots")
}
y15450K_dmrpre<-dmrff.pre(estimate = statsdfo$bacon.est,se=statsdfo$bacon.se,methylation = M,chr=statsdfo$CHR,pos=statsdfo$MAPINFO)
y15_450K_betas<-M
y15_450K_statsdfo<-statsdfo
##############
#y9 450K
################
man <- readRDS("path to hm450.hg19.manifest.dataframe.rds")
man$ID<-row.names(man)
# # # Exclude bad probes and sex chromosomes
manSub <- man[man$seqnames %in% c("chrX","chrY","ID","chrM"),"ID"]
man<-man[!man$ID %in% manSub,c("ID","seqnames","start","end")]
man$CHR<-gsub("chr","",man$seqnames)
man$MAPINFO<-man$start

survey<-readRDS("path to file containing dependent variable and covariate data with array ID")
# Read in stats
y9_450K_stats<-read.delim("path to y9_450K_sensitivity_model1.ewasoutput.txt")

# Read in beta matrix for y9 & 15
beta <- qs_read("/drives/drive1/ff_chaya/Lisa/nc_betafiles/BetaMatrix.k450.n1724.qcfilter.020425.qs")
beta<-beta[,colnames(beta)%in% survey$Array_ID]

y9_450K_stats_m<-merge(man,y9_450K_stats,by.x="ID",by.y="probe")
y9_450K_stats_m$CHR<-as.numeric(as.character(y9_450K_stats_m$CHR))
y9_450K_stats_m$start<-as.numeric(as.character(y9_450K_stats_m$start))
y9_450K_stats_m$MAPINFO<-as.numeric(as.character(y9_450K_stats_m$MAPINFO))
y9_450K_stats_m$t<-as.numeric(as.character(y9_450K_stats_m$t))
y9_450K_stats_m$StdError<-as.numeric(as.character(y9_450K_stats_m$StdError))
y9_450K_stats_m$Coefficient<-as.numeric(as.character(y9_450K_stats_m$Coefficient))
y9_450K_stats_m<-bacon_fun(y9_450K_stats_m)

beta<-beta[row.names(beta)%in% y9_450K_stats_m$ID,]

y9_450K_stats<-y9_450K_stats_m[y9_450K_stats_m$ID %in% row.names(beta),]
statsdf<-y9_450K_stats
betas<-beta
outdir ="output directory file path"
outpre="output file prefix"

# Remove NA in chromosomes and position
statsdf<-statsdf[!is.na(statsdf$CHR),]
statsdf<-statsdf[!is.na(statsdf$MAPINFO),]
statsdfo<-statsdf[order(statsdf$CHR,statsdf$MAPINFO),]
statsdfo$bacon.fdr<-p.adjust(statsdfo$bacon.P.Value,method="fdr")

#betas <- readRDS(betas)
betas<-betas[row.names(betas) %in% statsdfo$ID,]
betaso<-betas[match(statsdfo$ID,row.names(betas)),]
M<- log2(betaso) - log2(1 - betaso)
rm(betas)
rm(betaso)
rm(statsdf)
gc()

dmrs<-dmrff(estimate = statsdfo$bacon.est,se=statsdfo$bacon.se,p.value = statsdfo$bacon.P.Value,methylation = M,chr=statsdfo$CHR,pos=statsdfo$MAPINFO)
dmrsfilt<-dmrs[which(dmrs$p.adjust<0.1 & dmrs$n>1),]
if(nrow(dmrsfilt)>0) {
  sites<-dmrff.sites(dmrsfilt,statsdfo$CHR,statsdfo$MAPINFO)
  #  sites<-cbind(sites,statsdfo[sites$site,c("Gene","probe","start","end","bacon.P.Value","bacon.est")])
  sites<-cbind(sites,statsdfo[sites$site,c("ID","CHR","MAPINFO","bacon.fdr","bacon.P.Value","bacon.est")])
  sites<-cbind(sites,dmr=dmrsfilt[sites$region,c("chr","start","end","n","B","estimate","se","z","p.value","number.tests","p.adjust")])
  write.table(sites,file=paste0(outdir,outpre,"_sigdmrffoutput.txt"),sep="\t",row.names=F,quote=F)
  write.table(dmrs,file=paste0(outdir,outpre,"_alldmrffoutput.txt"),sep="\t",row.names=F,quote=F)
  
  
  pdf(paste0(outdir,outpre,"_sigdmrffoutput_plots.pdf"))
  for(i in 1:nrow(dmrsfilt)){
    print(dmrff.plot(region.chr = dmrsfilt[i,"chr"],region.start = dmrsfilt[i,"start"],region.end=dmrsfilt[i,"end"],estimate = statsdfo$bacon.est,statsdfo$bacon.se,statsdfo$CHR,statsdfo$MAPINFO,ci=T,expand=1))
  }
  dev.off()
} else {
  message("No significant DMRs found - skipping output and plots")
}

y9_450K_betas<-M
y9_statsdfo<-statsdfo
y9_450K_statsdfo<-y9_statsdfo
####### 
#harmonize datasets
###################
y9CpGs<-row.names(y9epic_betas) [row.names(y9epic_betas) %in% row.names(y9_450K_betas)]
y15CpGs<-row.names(y15epic_betas) [row.names(y15epic_betas) %in% row.names(y15_450K_betas)]

# y9
y9_epic_M<-y9epic_betas[row.names(y9epic_betas)%in% y9CpGs,]
y9_epic_Mord<-y9_epic_M[match(y9CpGs,row.names(y9_epic_M)),]

y9_450K_M<-y9_450K_betas[row.names(y9_450K_betas)%in% y9CpGs,]
y9_450K_Mord<-y9_450K_M[match(y9CpGs,row.names(y9_450K_M)),]

y9epic_statsdfo_sub<-y9epic_statsdfo[y9epic_statsdfo$ID %in% y9CpGs,]
y9epic_statsdfo_sub_ord<-y9epic_statsdfo_sub[match(y9CpGs,y9epic_statsdfo_sub$ID),]

y9_450K_statsdfo_sub<-y9_450K_statsdfo[y9_450K_statsdfo$ID %in% y9CpGs,]
y9_450K_statsdfo_sub_ord<-y9_450K_statsdfo_sub[match(y9CpGs,y9_450K_statsdfo_sub$ID),]

# y15
y15_epic_M<-y15epic_betas[row.names(y15epic_betas)%in% y15CpGs,]
y15_epic_Mord<-y15_epic_M[match(y15CpGs,row.names(y15_epic_M)),]

y15_450K_M<-y15_450K_betas[row.names(y15_450K_betas)%in% y15CpGs,]
y15_450K_Mord<-y15_450K_M[match(y15CpGs,row.names(y15_450K_M)),]

y15epic_statsdfo_sub<-y15epic_statsdfo[y15epic_statsdfo$ID %in% y15CpGs,]
y15epic_statsdfo_sub_ord<-y15epic_statsdfo_sub[match(y15CpGs,y15epic_statsdfo_sub$ID),]

y15_450K_statsdfo_sub<-y15_450K_statsdfo[y15_450K_statsdfo$ID %in% y15CpGs,]
y15_450K_statsdfo_sub_ord<-y15_450K_statsdfo_sub[match(y15CpGs,y15_450K_statsdfo_sub$ID),]


# Create dmrffpre objects
dmrffpre_fun<-function(df,M){
  df_dmrpre<-dmrff.pre(estimate = df$bacon.est,se=df$bacon.se,methylation = M,chr=df$CHR,pos=df$MAPINFO)
  return(df_dmrpre)
}

y9_450K_dmrpre<-dmrffpre_fun(df=y9_450K_statsdfo_sub_ord,M=y9_450K_Mord)
y9_epic_dmrpre<-dmrffpre_fun(df=y9epic_statsdfo_sub_ord,M=y9_epic_Mord)

y15_450K_dmrpre<-dmrffpre_fun(df=y15_450K_statsdfo_sub_ord,M=y15_450K_Mord)
y15_epic_dmrpre<-dmrffpre_fun(df=y15epic_statsdfo_sub_ord,M=y15_epic_Mord)



saveRDS(y9_450K_dmrpre,file="output file")

saveRDS(y15_450K_dmrpre,file="output file")

saveRDS(y15_epic_dmrpre,file="output file")

saveRDS(y9_epic_dmrpre,file="/output file")

saveRDS(y22dmrff,file="output file")

####################################
# y9 meta
######################
metalist<-list(epic=y9_epic_dmrpre,hm450=y9_450K_dmrpre)
meta<-dmrff.meta(metalist)
saveRDS(meta,file="output file")
meta_sigdmrs<-meta[[2]][meta[[2]]$n>1 & meta[[2]]$p.adjust<0.1,]
meta_sigdmrssites<-dmrff.sites(meta_sigdmrs,y9epic_statsdfo_sub_ord$CHR,y9epic_statsdfo_sub_ord$MAPINFO)
meta_sigdmrssites<-cbind(meta_sigdmrssites,y9epic_statsdfo_sub_ord[meta_sigdmrssites$site,])
meta_sigdmrssites<-cbind(meta_sigdmrssites,meta_sigdmrs[meta_sigdmrssites$region,c("chr","start","end","n","B","estimate","se","z","p.value","number.tests","p.adjust")])
write.table(meta_sigdmrssites,file="output file",sep="\t",row.names=F,quote=F)
y9meta_sigdmrsites<-meta_sigdmrssites
##############
# y15 meta
###################
metalist<-list(epic=y15_epic_dmrpre,hm450=y15_450K_dmrpre)
meta<-dmrff.meta(metalist)
saveRDS(meta,file="output file")
meta_sigdmrs<-meta[[2]][meta[[2]]$n>1 & meta[[2]]$p.adjust<0.1,]
meta_sigdmrssites<-dmrff.sites(meta_sigdmrs,y9epic_statsdfo_sub_ord$CHR,y9epic_statsdfo_sub_ord$MAPINFO)
meta_sigdmrssites<-cbind(meta_sigdmrssites,y9epic_statsdfo_sub_ord[meta_sigdmrssites$site,])
meta_sigdmrssites<-cbind(meta_sigdmrssites,meta_sigdmrs[meta_sigdmrssites$region,c("chr","start","end","n","B","estimate","se","z","p.value","number.tests","p.adjust")])
write.table(meta_sigdmrssites,file="output file",sep="\t",row.names=F,quote=F)
y15meta_sigdmrsites<-meta_sigdmrssites
############# 
### y22 with and without subsetting meta CpGs
man <- readRDS("path to EPIC.hg19.manifest.dataframe.rds")
man$ID<-row.names(man)
# # # Exclude bad probes and sex chromosomes
manSub <- man[man$seqnames %in% c("chrX","chrY","chrM") | !substr(man$ID,1,2) %in% "cg","ID"]
man<-man[!man$ID %in% manSub,c("ID","seqnames","start","end")]
man$CHR<-gsub("chr","",as.character(man$seqnames))
man$MAPINFO<-man$start

y22<-readRDS("path to y22_epic.ewasoutput_man_annotated_chipseeker.RDS")
y22beta <- qs_read("path to BetaMatrix.blood.n1271.qcfilter.020425.qs") # y22
y22survey<-readRDS("path to file containing dependent variable and covariate data with array ID")
y22beta<-y22beta[,colnames(y22beta) %in% survey$Array_ID]
################
y22_epic_stats_m<-y22
y22_epic_stats_m$CHR<-as.numeric(as.character(gsub("chr","",y22_epic_stats_m$seqnames)))

y22beta<-y22beta[row.names(y22beta)%in% y22_epic_stats_m$probe,]

y22_epic_stats<-y22_epic_stats_m[y22_epic_stats_m$probe %in% row.names(y22beta),]
statsdf<-y22_epic_stats
betas<-y22beta
outdir ="path to output directory"
outpre="output file prefix"

statsdf$chr<-as.numeric(gsub("chr","",as.character(statsdf$seqnames)))
# Remove NA in chromosomes and position
statsdf<-statsdf[!is.na(statsdf$CHR),]
statsdf$MAPINFO<-statsdf$start
statsdf<-statsdf[!is.na(statsdf$MAPINFO),]
statsdfo<-statsdf[order(statsdf$CHR,statsdf$MAPINFO),]
statsdfo$Coefficient<-as.numeric(statsdfo$Coefficient)
statsdfo$t<-as.numeric(statsdfo$t)
statsdfo$StdError<-as.numeric(statsdfo$StdError)
statsdfo<-bacon_fun(statsdfo)
statsdfo$bacon.fdr<-p.adjust(statsdfo$bacon.P.Value,method="fdr")


#betas <- readRDS(betas)
betas<-betas[row.names(betas) %in% statsdfo$probe,]
betaso<-betas[match(statsdfo$probe,row.names(betas)),]
M<- log2(betaso) - log2(1 - betaso)

rm(betas)
rm(betaso)
rm(statsdf)
gc()
dmrs<-dmrff(estimate = statsdfo$bacon.est,se=statsdfo$bacon.se,p.value = statsdfo$bacon.P.Value,methylation = M,chr=statsdfo$CHR,pos=statsdfo$MAPINFO)

dmrsfilt<-dmrs[which(dmrs$p.adjust<0.1 & dmrs$n>1),]
sites<-dmrff.sites(dmrsfilt,statsdfo$CHR,statsdfo$MAPINFO)
#  sites<-cbind(sites,statsdfo[sites$site,c("Gene","probe","start","end","bacon.P.Value","bacon.est")])
sites<-cbind(sites,statsdfo[sites$site,c("probe","CHR","MAPINFO","bacon.fdr","bacon.P.Value","bacon.est")])
sites<-cbind(sites,dmr=dmrsfilt[sites$region,c("chr","start","end","n","B","estimate","se","z","p.value","number.tests","p.adjust")])
write.table(sites,file=paste0(outdir,outpre,"_sigdmrffoutput.txt"),sep="\t",row.names=F,quote=F)
write.table(dmrs,file=paste0(outdir,outpre,"_alldmrffoutput.txt"),sep="\t",row.names=F,quote=F)


pdf(paste0(outdir,outpre,"_sigdmrffoutput_plots.pdf"))
for(i in 1:nrow(dmrsfilt)){
  print(dmrff.plot(region.chr = dmrsfilt[i,"chr"],region.start = dmrsfilt[i,"start"],region.end=dmrsfilt[i,"end"],estimate = statsdfo$bacon.est,statsdfo$bacon.se,statsdfo$CHR,statsdfo$MAPINFO,ci=T,expand=1))
}
dev.off()

# subset to meta cpgs
y22meta_betas<-M[row.names(M)%in% y15CpGs,]
y22meta_statsdfo<-statsdfo[statsdfo$probe %in% y15CpGs,]
outpre="output file prefix"


dmrs<-dmrff(estimate = y22meta_statsdfo$bacon.est,se=y22meta_statsdfo$bacon.se,p.value = y22meta_statsdfo$bacon.P.Value,methylation = y22meta_betas,chr=y22meta_statsdfo$CHR,pos=y22meta_statsdfo$MAPINFO)

dmrsfilt<-dmrs[which(dmrs$p.adjust<0.1 & dmrs$n>1),]
sites<-dmrff.sites(dmrsfilt,statsdfo$CHR,statsdfo$MAPINFO)
sites<-cbind(sites,statsdfo[sites$site,c("probe","CHR","MAPINFO","bacon.fdr","bacon.P.Value","bacon.est")])
sites<-cbind(sites,dmr=dmrsfilt[sites$region,c("chr","start","end","n","B","estimate","se","z","p.value","number.tests","p.adjust")])
write.table(sites,file=paste0(outdir,outpre,"_sigdmrffoutput.txt"),sep="\t",row.names=F,quote=F)
write.table(dmrs,file=paste0(outdir,outpre,"_alldmrffoutput.txt"),sep="\t",row.names=F,quote=F)


pdf(paste0(outdir,outpre,"_sigdmrffoutput_plots.pdf"))
for(i in 1:nrow(dmrsfilt)){
  print(dmrff.plot(region.chr = dmrsfilt[i,"chr"],region.start = dmrsfilt[i,"start"],region.end=dmrsfilt[i,"end"],estimate = statsdfo$bacon.est,statsdfo$bacon.se,statsdfo$CHR,statsdfo$MAPINFO,ci=T,expand=1))
}
dev.off()
saveRDS(session,"output file path /sessioninfo.RDS")
