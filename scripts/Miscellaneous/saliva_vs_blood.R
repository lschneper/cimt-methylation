
# Author: Ayanava Ganguly
# Edited by: Eve Schoeffler
# Edited by: Lisa Schneper
.libPaths("path to R library directory")

library(magrittr)
library(dplyr)
library(tidyr)
library(gtools)
library(readxl)
library(stringr)
library(qs2)
library(GenomicRanges)
library(irr)
library(limma)
library(missMDA)

outdir="output directory path"
# read in Zhou annotation manifest
reference <- as.data.frame(readRDS("path to EPIC.hg19.manifest.rds"))
 
bloodids<-readRDS("path to RDS containing Age 22 blood methylation sample sheet")

salivaids<-readRDS("path to Age 22 saliva sample sheet RDS")
blood<-qs_read("path to Age 22 blood beta matrix")
blood<-as.data.frame(blood)
blood<-blood[,names(blood) %in% bloodids$Array_ID]
bloodids<-bloodids[bloodids$Array_ID %in% names(blood),]
# Create dataframe with IDs matching order of methyl array
bloodidso<-bloodids[match(names(blood),bloodids$Array_ID),]

saliva<-qs_read("path to Age 22 saliva beta matrix")
saliva<-as.data.frame(saliva)
saliva<-saliva[,names(saliva) %in% salivaids$Array_ID]
#saliva<-saliva[,!names(saliva) %in% badsamples$Array_ID]
salivaids<-salivaids[salivaids$Array_ID %in% names(saliva),]
salivaidso<-salivaids[match(names(saliva),salivaids$Array_ID),]
salivacommon<-salivaidso[salivaidso$CVID %in% bloodidso$CVID,]
bloodcommon<-bloodidso[bloodidso$CVID %in% salivaidso$CVID,]

bloodcommon<-bloodcommon[,c("CVID","Array_ID")]
bloodidso<-bloodidso[bloodidso$Array_ID %in% bloodcommon$Array_ID,]

blood<-blood[,names(blood) %in% bloodidso$Array_ID]

salivaidso<-salivaidso[salivaidso$Array_ID %in% salivacommon$Array_ID,]

saliva<-saliva[,names(saliva) %in% salivaidso$Array_ID]
salivacommon<-salivacommon[,c("CVID","Array_ID")]
commonkey<-merge(bloodcommon,salivacommon,by="CVID")
names(commonkey)[2:3]<-c("blood","saliva")
bloodo<-blood[,match(commonkey$blood,names(blood))]
salivao<-saliva[,match(commonkey$saliva,names(saliva))]

interleaved_cols<-as.vector(rbind(names(bloodo),names(salivao)))
bloodofilt<-bloodo[row.names(bloodo) %in% row.names(salivao),]
salivaofilt<-salivao[row.names(salivao) %in% row.names(bloodo),]
# methyl_w7_blood <- log2((methyl_y22) / (1 - methyl_y22))
bloodsaliva<-cbind(bloodofilt,salivaofilt)[,interleaved_cols]
bloodsaliva_inorder <- cbind(bloodofilt,salivaofilt)


# remove all CpGs that are masked or on non-autosomal chromosomes
# Keep masked, remove sex
reference$ID <- rownames(reference)

refSub <- c(reference[reference$seqnames == "chrX","ID"],
            reference[reference$seqnames == "chrY","ID"],
            reference[reference$seqnames == "chrM","ID"])

bloodsaliva <- bloodsaliva[!rownames(bloodsaliva) %in% refSub,]
bloodsalivaM<-log2((bloodsaliva) / (1 - bloodsaliva))

bloodsaliva_inorder <- bloodsaliva_inorder[!rownames(bloodsaliva_inorder)%in% refSub,]
bloodsalivaM_order<-log2((bloodsaliva_inorder) / (1 - bloodsaliva_inorder))

# Run EWAS for tissue
design_pre_blood <- commonkey[, c(1,2)]
design_pre_blood$Blood <- 1
design_pre_saliva <- commonkey[, c(1,3)]
design_pre_saliva$Blood <- 0
names(design_pre_blood)[2] <- "Array_ID"
names(design_pre_saliva)[2] <- "Array_ID"
design_pre <- rbind(design_pre_blood, design_pre_saliva)
rownames(design_pre) <- design_pre$Array_ID
design_pre$Blood <- as.factor(design_pre$Blood)
design_tissue <- model.matrix(~design_pre$Blood)
block <- factor(design_pre$CVID) 
bloodsalivaM_order2<-as.matrix(bloodsalivaM_order)
storage.mode(bloodsalivaM_order2)<-"double"
corFit <- duplicateCorrelation(bloodsalivaM_order2, design_tissue, block = block)

dir.create(outdir)
saveRDS(rownames(bloodsalivaM_order), paste0(outdir,"/saliva_blood_corFit.CpG_order.RDS"))
saveRDS(corFit, paste0(outdir,"/saliva_blood_corFit.RDS"))
fit <- lmFit(bloodsalivaM_order, design_tissue, block = block, correlation = corFit$consensus.correlation)
efit <- eBayes(fit)
saveRDS(efit, paste0(outdir,"/saliva_blood_eFit.RDS"))
top <- topTable(efit, number = Inf)
top$FC <- 2^top$logFC

#####################
# Start from salivaofilt & bloodofilt
bloodofiltM<-log2((bloodofilt) / (1 - bloodofilt))
salivaofiltM<-log2((salivaofilt) / (1 - salivaofilt))

pearsoncortest_list<-list()
pearsoncorrelationestvector<-vector()
pearsoncorrelationpval<-vector()
for (i in 1:nrow(bloodofiltM)){
  pearsoncortest_list[[i]]<-try(cor.test(as.numeric(bloodofiltM[i,]),as.numeric(salivaofiltM[i,]),na.rm=T))
  names(pearsoncortest_list)[i]<-row.names(bloodofiltM)[i]}
saveRDS(pearsoncortest_list,file=paste0(outdir,"/CpGlevel_pearson_correlation.rds"))

for (i in 1:length(pearsoncortest_list)){
  pearsoncorrelationestvector[i]<-try(pearsoncortest_list[[i]]$estimate)
  pearsoncorrelationpval[i]<-try(pearsoncortest_list[[i]]$p.value)}

pearsoncorrelationdf<-data.frame(cpg=names(pearsoncortest_list),estimate=as.numeric(pearsoncorrelationestvector),pval=as.numeric(pearsoncorrelationpval))
saveRDS(pearsoncorrelationdf,file=paste0(outdir,"/CpGlevel_pearsoncorrelation_df.rds"))
######

pearsoncorrelationestvectordf<-as.data.frame(pearsoncorrelationestvector)
names(pearsoncorrelationestvectordf)<-"estimate"
ggplot(pearsoncorrelationestvectordf,aes(x=as.numeric(estimate)))+geom_histogram(bins=20,fill="lightgrey",color="black")+theme_bw()

spearmancortest_list<-list()
spearmancorrelationestvector<-vector()
spearmancorrelationpval<-vector()
for (i in 1:ncol(bloodofiltM)){
  spearmancortest_list[[i]]<-cor.test(bloodofiltM[,i],salivaofiltM[,i],na.rm=T,method = "spearman")
  spearmancorrelationestvector[i]<-pearsoncortest_list[[i]]$estimate
  spearmancorrelationpval[i]<-cortest_list[[i]]$p.value}
spearmancorrelationestvectordf<-as.data.frame(spearmancorrelationestvector)
names(spearmancorrelationestvectordf)<-"estimate"
ggplot(spearmancorrelationestvectordf,aes(x=estimate))+geom_histogram(bins=20,fill="lightgrey",color="black")+theme_bw()


bloodofiltM$tissue<-"blood"
salivaofiltM$tissue<-"saliva"

bloodsalivaM<-cbind(bloodofiltM[,1],salivaofiltM[,1])
colnames(bloodsalivaM)<-c("blood","saliva")
row.names(bloodsalivaM)<-row.names(bloodofilt)

ggplot(bloodsalivaM,aes(x=blood,y=saliva))+geom_point()+geom_smooth(method="lm")+theme_bw()+xlim(-15,15)+ylim(-15,15)

#######################
# Compute ICC
library(irr)
icc_results<-numeric(nrow(bloodsalivaM))

for(i in 1:nrow(bloodsalivaM)){
  cpg_values <- bloodsalivaM[i,]
  cpg_pairs<-as.data.frame(matrix(nrow=78, ncol=2))
  for (j in 1:78){
    sample1<-commonkey$blood[j]
    sample2<-commonkey$saliva[j]
    cpg_pairs[j,1]<-cpg_values[sample1]
    cpg_pairs[j,2]<-cpg_values[sample2]
  }
  # Compute ICC(3,1): two-way mixed, consistency, single measures
  icc_out<-icc(cpg_pairs,model="twoway",type="consistency",unit="single")
  icc_results[i] <-icc_out$value
}
names(icc_results)<-row.names(bloodsalivaM)
saveRDS(icc_results,file=paste0(outdir,"icc_results_Mvalues.rds"))

acceptableCpGs<-names(icc_results)[icc_results > 0.5]
saveRDS(acceptableCpGs,file=paste0(outdir,"acceptableCpGs.RDS"))
##### 
# Calculate Bland Altman statistics
results_list<-vector("list", nrow(bloodofilt))
bloodofiltmeth<-log2((bloodofilt)/(1-bloodofilt))
salivaofiltmeth<-log2((salivaofilt)/(1-salivaofilt))

#Loop through each CpG site - use beta values

valid_rows <- rowSums(!is.na(bloodofilt) & !is.na(salivaofilt)) >= 2
bloodofilt2 <- bloodofilt[valid_rows,]
salivaofilt2<-salivaofilt[valid_rows,]

results_list <- vector("list", nrow(bloodofilt2))
for(i in 1:nrow(bloodofilt2)){
  blood_vals<-as.numeric(bloodofilt2[i,])
  saliva_vals<-as.numeric(salivaofilt2[i,])
  
}  
# Compute ICC
library(irr)
icc_results<-numeric(nrow(bloodsalivaM))

for(i in 1:nrow(bloodsalivaM)){
  cpg_values <- bloodsalivaM[i,]
  cpg_pairs<-as.data.frame(matrix(nrow=78, ncol=2))
  for (j in 1:78){
    sample1<-commonkey$blood[j]
    sample2<-commonkey$saliva[j]
    cpg_pairs[j,1]<-cpg_values[sample1]
    cpg_pairs[j,2]<-cpg_values[sample2]
  }
  # Compute ICC(3,1): two-way mixed, consistency, single measures
  icc_out<-icc(cpg_pairs,model="twoway",type="consistency",unit="single")
  icc_results[i] <-icc_out$value
}
names(icc_results)<-row.names(bloodsalivaM)
saveRDS(icc_results,file=paste0(outdir,"icc_results_Mvalues.rds"))

acceptableCpGs<-names(icc_results)[icc_results > 0.5]
saveRDS(acceptableCpGs,file=paste0(outdir,"acceptableCpGs.RDS"))
##### 
# Calculate Bland Altman statistics
results_list<-vector("list", nrow(bloodofilt))
bloodofiltmeth<-log2((bloodofilt)/(1-bloodofilt))
salivaofiltmeth<-log2((salivaofilt)/(1-salivaofilt))

#Loop through each CpG site - use beta values
valid_rows <- rowSums(!is.na(bloodofilt) & !is.na(salivaofilt)) >= 2
bloodofilt2 <- bloodofilt[valid_rows,]
salivaofilt2<-salivaofilt[valid_rows,]

for(i in 1:nrow(bloodofilt2)){
  blood_vals<-as.numeric(bloodofilt2[i,])
 saliva_vals<-as.numeric(salivaofilt2[i,])
}  
  
  
  # Skip CpG if all values are NA or lengths don't match
  # if(all(is.na(blood_vals)) || all(is.na(saliva_vals)) || length(blood_vals) != length(saliva_vals)){
  # Calculate differences and stats
  diff_vals<-salivaofilt2 - bloodofilt2
  non_na_diffs <- rowSums(!is.na(diff_vals))
  goodCpGs<-row.names(diff_vals)[non_na_diffs>=70]
  diff_vals<-diff_vals[row.names(diff_vals)%in% goodCpGs,]
  

  mean_diff <- apply(diff_vals,1,mean,na.rm=T)
  sd_diff <- apply(diff_vals,1,sd,na.rm=T)
  loa_lower <- mean_diff - 1.96 * sd_diff
  loa_upper <- mean_diff + 1.96 * sd_diff
  max_abs_diff <- abs(apply(diff_vals,1,mean,na.rm=T))
  
  # Optional: flag CpGs with big bias
  flagged <-  abs(mean_diff) > 0.1
  

  
  # Store results
  bland_altman_summary <- data.frame(
    CpG_ID = names(mean_diff),
    Mean_Diff = mean_diff,
    SD_Diff = sd_diff,
    LoA_Lower = loa_lower,
    LoA_Upper = loa_upper,
    Max_Abs_Diff = max_abs_diff,
    Flagged = flagged
  )
saveRDS(bland_altman_summary,paste0(outdir,"/bland_altman_summary.RDS"))
CairoPDF(paste0(outdir,"/bland_altman_meandiff.pdf",family="Arial"))
#ggplot(bland_altman_summary,aes(x=Mean_Diff))+geom_histogram() +theme_bw()+xlab("Mean Difference")
ggplot(bland_altman_summary,aes(x=Mean_Diff))+
  geom_vline(xintercept = c(-0.2, 0.2),
             #xintercept = c(interval[1], interval[2]), 
             linetype = "dashed",
             color = "darkgrey")+
  geom_histogram() +
  theme_classic()+
  xlab("Mean Difference in Beta Values") +
  ylab("Frequency") +
  ggtitle("Distribution of Difference in Beta Values", subtitle = "Between Saliva and Blood") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 13),
        axis.title = element_text(size = 12),
        legend.position = NULL) 
dev.off()
saveRDS(bloodofilt2,file=paste0(outdir,"bloodofilt.RDS"))
saveRDS(salivaofilt2,file=paste0(outdir,"salivaofilt.RDS"))

########################
fit <- lmFit(bloodsalivaM_order, design_tissue, block = block, correlation = corFit$consensus.correlation)
efit <- eBayes(fit)
top <- topTable(efit, number = Inf)
topf<-top[row.names(top)%in% bland_altman_summary$CpG_ID,]
topo<-topf[match(bland_altman_summary$CpG_ID, row.names(topf)),]
topo$Mean_Diff <- bland_altman_summary$Mean_Diff
saveRDS(topo, paste0(outdir,"/saliva_blood_limma_results.RDS"))
bland_altman_summary$newflag<-ifelse(topo$adj.P.Val < 0.001 & abs(topo$Mean_Diff) > 0.1,1,0)
saveRDS(dif, paste0(outdir,"/saliva_blood_limma_results.difMethOnly.RDS"))
saveRDS(rownames(dif), paste0(outdir,"/saliva_blood_comparision/saliva_blood_difMethCpGs.RDS"))

# Compute ICC
# library(psych)
library(irr)
icc_results<-numeric(nrow(bloodsalivaM))

for(i in 1:nrow(bloodsalivaM)){
  cpg_values <- bloodsalivaM[i,]
  cpg_pairs<-as.data.frame(matrix(nrow=78, ncol=2))
  for (j in 1:78){
    sample1<-commonkey$blood[j]
    sample2<-commonkey$saliva[j]
    cpg_pairs[j,1]<-cpg_values[sample1]
    cpg_pairs[j,2]<-cpg_values[sample2]
  }
  # Compute ICC(3,1): two-way mixed, consistency, single measures
  icc_out<-icc(cpg_pairs,model="twoway",type="consistency",unit="single")
  icc_results[i] <-icc_out$value
}
names(icc_results)<-row.names(bloodsalivaM)
saveRDS(icc_results,file="icc_results_Mvalues.rds")

acceptableCpGs<-names(icc_results)[icc_results > 0.5]
saveRDS(acceptableCpGs,file="acceptableCpGs.RDS")
##### 
# Calculate Bland Altman statistics
results_list<-vector("list", nrow(bloodofilt))
bloodofiltmeth<-log2((bloodofilt)/(1-bloodofilt))
salivaofiltmeth<-log2((salivaofilt)/(1-salivaofilt))

#Loop through each CpG site - use beta values

valid_rows <- rowSums(!is.na(bloodofilt) & !is.na(salivaofilt)) >= 2
bloodofilt2 <- bloodofilt[valid_rows,]
salivaofilt2<-salivaofilt[valid_rows,]

results_list <- vector("list", nrow(bloodofilt2))
for(i in 1:nrow(bloodofilt2)){
  blood_vals<-as.numeric(bloodofilt2[i,])
  saliva_vals<-as.numeric(salivaofilt2[i,])
  
  
  
  # Skip CpG if all values are NA or lengths don't match
  # if(all(is.na(blood_vals)) || all(is.na(saliva_vals)) || length(blood_vals) != length(saliva_vals)){
  # Calculate differnces and stats
  diff_vals<-saliva_vals - blood_vals
  non_na_diffs <- diff_vals[!is.na(diff_vals)]
  
  if( length(non_na_diffs) >= 2){
    mean_diff <- mean(non_na_diffs)
    sd_diff <- sd(non_na_diffs)
    loa_lower <- mean_diff - 1.96 * sd_diff
    loa_upper <- mean_diff + 1.96 * sd_diff
    max_abs_diff <- max(abs(non_na_diffs))
  } else {
    mean_diff <- NA
    sd_diff <-NA
    loa_lower <-NA
    loa_upper <- NA
    max_abs_diff <- NA
  }
  
  # Optional: flag CpGs with big bias
  flagged <- ! is.na(mean_diff) && abs(mean_diff) > 0.1
    
  # Store results
  results_list[[i]] <- data.frame(
    CpG_ID = rownames(bloodofilt2)[i],
    Mean_Diff = mean_diff,
    SD_Diff = sd_diff,
    LoA_Lower = loa_lower,
    LoA_Upper = loa_upper,
    Max_Abs_Diff = max_abs_diff,
    Flagged = flagged
  )
}

saveRDS(bloodofilt,file="bloodofilt.RDS")
saveRDS(salivaofilt,file="salivaofilt.RDS")


methyl_w7_saliva_common<-log2(salivaofilt2 / (1 - salivaofilt2))
methyl_w7_blood_common<-log2(bloodofilt2 / (1 - bloodofilt2))

# Calculate correlation between saliva and blood
corel<-cor(methyl_w7_blood_common,methyl_w7_saliva_common,use ="pairwise.complete.obs")

#sal_blood_cor$estimate #0.9672289 
#sal_blood_cor$p.value #0

# Check the difference in M-values between blood and saliva
delta_M <- (methyl_w7_blood_common - methyl_w7_saliva_common)
all_delta_M <- as.vector(as.matrix(delta_M))
sd_delta_M <- sd(all_delta_M, na.rm = T)
mean_delta_M <- mean(all_delta_M, na.rm = T)
df_delta_M <- data.frame(delta_M = all_delta_M)
ggplot(df_delta_M, aes(x = delta_M)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = c(mean_delta_M - 2 * sd_delta_M, mean_delta_M + 2 * sd_delta_M),
             color = "darkgreen", linetype = "dashed") +
  # theme_classic() +
  labs(x = "Distribution of Difference in M-Values",
       y = "Frequency") 

# PCA Analysis

library(stats)

methyl_w7_combined <- cbind(methyl_w7_blood_common, methyl_w7_saliva_common)
imputed <- imputePCA(t(methyl_w7_combined), npc = 5)
pca_imputed <- imputed$completeObs
pca_result <- prcomp(pca_imputed, scale. = T)
pc_scores <- as.data.frame(pca_result$x)
pc_scores$Sample_Type <- c(rep("Blood", ncol(methyl_w7_blood_common)), rep("Saliva", ncol(methyl_w7_saliva_common)))

CairoPDF(paste0(outdir,"/cross_tissue_PCA.pdf"),family="Arial",width=6,height=6)
ggplot(pc_scores, aes(x = PC1, y = -PC2, color = Sample_Type)) + 
  geom_point(size = 2) + 
  labs(title = "Cross-tissue PCA", color = "Tissue") + 
  theme_classic() +
  theme(axis.title.x = element_text(size = 14, face = "bold"), 
        axis.title.y = element_text(size = 14, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  scale_color_manual(values = c( "#FB9A99","#A6CEE3"))+
  xlab("PC1")+ylab("PC2")

dev.off()
newbad<-bland_altman_summary[abs(bland_altman_summary$Mean_Diff)>0.2,]
bland_altman_summary$Flagged<-ifelse(abs(bland_altman_summary$Mean_Diff)>0.2,"TRUE","FALSE")
saveRDS(newbad$CpG_ID,file=paste0(outdir,"/saliva_blood_difMethCpGs.RDS"))
write.table(newbad$CpG_ID,file=paste0(outdir,"/saliva_blood_difMethCpGs.txt"),sep="\t",row.names=F,quote=F)
write.table(bland_altman_summary,file=paste0(outdir,"/bland_altman_summary.txt"),sep="\t",row.names=F,quote=F)
saveRDS(bland_altman_summary,file=paste0(outdir,"/bland_altman_summary.RDS"))

        