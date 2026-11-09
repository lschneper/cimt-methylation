.libPaths("/drives/drive1/home/lweis/R/x86_64-redhat-linux-gnu-library/4.4")
library(ggplot2)
library(MASS)
library(colorBlindness)
library(Cairo)
# Correlation between race/ethnicity and ancestry
survey1 <- readRDS("/drives/drive1/ff_chaya/Lisa/nc_samplesheets/y9_450K_ses_toast_ancPC_032726.RDS")
survey2 <- readRDS("/drives/drive1/ff_chaya/Lisa/nc_samplesheets/y9_epic_ses_toast_ancPC_032726.RDS")
survey3 <- readRDS("/drives/drive1/ff_chaya/Lisa/nc_samplesheets/y15_450K_ses_toast_ancPC_032726.RDS")
survey4 <- readRDS("/drives/drive1/ff_chaya/Lisa/nc_samplesheets/y15_epic_ses_toast_ancPC_032726.RDS")
survey5<-readRDS("/drives/drive1/ff_chaya/Lisa/nc_samplesheets/y15_epic_ses_toast_ancPC_032726.RDS")
survey5<-survey5[!is.na(survey5$Mean_Mean_CIMT),]
survey5<-survey5[,c("CVID","raceEthY22","GPC1","GPC2","GPC3","GPC4","GPC5")]
all<-rbind(survey1,survey2,survey3[,1:29],survey4[,1:29])
all<-all[!is.na(all$Mean_Mean_CIMT),]
all<-all[,c(1,19,25:29)]
all<-rbind(all,survey5)
all<-all[!duplicated(all),]
#covar<-read.delim("/drives/drive1/ff_chaya/Lisa/MR_code/mrfiltering_topmed/covardf_t_snp.txt")
#names(covar)[1]<-"covariate"
#names(covar)<-gsub("X","",names(covar))
#row.names(covar)<-covar$covariate
#covar<-covar[,-1]
#covar<-covar[,names(covar) %in% all$CVID]
survey<-all
survey<-survey[!is.na(survey$GPC1),]
eur<-survey[which(survey$GPC1 < -0.019 & survey$GPC2> 0.04),]
afr_lm<-survey[which(survey$GPC1 > 0.04 & survey$GPC2< 0),]
his<-survey[which(survey$GPC1 < -0.018 & survey$GPC2 < 0.04 & survey$GPC2 > -0.05),]
#eur<-survey[which(survey$GPC1> 0.018 & survey$GPC2> -0.0075),]
#eur<-survey[which(survey$GPC1> 0.018 & survey$GPC2> -0.0075),]
#bio<-read.csv("/drives/drive1/HPA/HPA_group_datasets/PGS/FFCW_PGS_and_PC_for_PrincetonFF_FamilyID.csv",header=T,stringsAsFactors=F)
survey$test<-ifelse(survey$CVID %in% eur$CVID, "EUR",ifelse(survey$CVID %in% afr_lm$CVID,"AFR",ifelse(survey$CVID %in% his$CVID,"HIS","other")))
survey$race<-ifelse(survey$raceEthY22=="WHITE","NHW",ifelse(survey$raceEthY22=="BLACK","NHB",ifelse(survey$raceEthY22=="HISPANIC","HISPANIC","OTHER")))
survey$race<-as.factor(survey$race)
survey$race<-relevel(survey$race,ref = "NHW")

manova_model<-manova(cbind(GPC1,GPC2,GPC3)~race,data=survey)
summary(manova_model)


lda_model<-lda(race~GPC1+GPC2+GPC3,data=survey)
print(lda_model)
pred<-predict(lda_model)$class
conf_matrix<-table(Predicted=pred,Actual=survey$race)
accuracy<-sum(diag(conf_matrix))/sum(conf_matrix)
cat("LDA Classification Accuracy:",accuracy,"\n")
#install.packages("colorBlindness",contriburl="file:///mnt/modules/cran/src/contrib",repos=NULL)

cbpallette<-c("#E69F00","#56B4E9","#009E73","#CC79A7")
CairoPDF("/drives/drive1/HPA/bin/Methylation/Analysis/Troubleshooting/LS_eve_Final_code/y9y15_ancestry_race3.pdf",family="arial")
ggplot(survey,aes(x=GPC1,y=-GPC2, color=race))+geom_point()+ theme_minimal()+xlab("PC1")+ylab("PC2")+theme(legend.title=element_blank(),legend.text=element_text(family="arial",face="bold"),axis.title=element_text(family="arial",face="bold"),axis.text=element_text(family="arial")) + scale_color_manual(values=cbpallette)
dev.off()

CairoPDF("/drives/drive1/HPA/bin/Methylation/Analysis/Troubleshooting/LS_eve_Final_code/y9y15_ancestry_race2.pdf",family="arial")

cb<-ggplot(surveycc,aes(x=GPC1,y=-GPC2, color=race))+geom_point()+ theme_minimal()+xlab("PC1")+ylab("PC2")+theme(legend.title=element_blank(),legend.text=element_text(family="arial",face="bold"),axis.title=element_text(family="arial",face="bold"),axis.text=element_text(family="arial"))+scale_color_manual(values=c(cbpallette))
cvdPlot(cb)
dev.off()
