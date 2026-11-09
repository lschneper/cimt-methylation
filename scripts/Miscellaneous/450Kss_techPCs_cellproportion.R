# Recreate technical pcs and middleton pcs for 450K
# read in beta
beta <- readRDS("/drives/drive1/HPA/Methylation/RDS_files/k450/BetaMatrix.k450.n1724.qcfilter.020425.RDS")
.libPaths("~/R/x86_64-redhat-linux-gnu-library/4.5/")
library(ewastools)
cell_props <- estimateLC(beta, ref = 'saliva',constrained=F)
rownames(cell_props) <- colnames(beta)
cell_props$MethID <- colnames(beta)

bloodsurvey<-readRDS("/drives/drive1/HPA/Methylation/Analysis/CVH_survey.Year22.dummy-vars.032425.RDS")
survey <- readRDS("/drives/drive1/HPA/Methylation/Analysis/subset_CVH_survey.Year9-15.dummy-vars.032425.RDS")
daterec<-read_excel("/drives/drive1/HPA/samples_keys/ffcore_datereceived.xlsx")
corecv_cw<-read.table("/drives/drive1/HPA/samples_keys/Y22ID_CVID.txt",sep="\t",header=T,stringsAsFactors=F)
corecv_cw$ffam<-gsub("Y","",corecv_cw$Y22LabID)
corecv_cw<-corecv_cw[,c(1,4)]
ffsubm<-merge(corecv_cw,daterec,by.x="ffam",by.y="westatID")
ff<-readRDS("/drives/drive1/ff_chaya/R_versions/FF_allwaves_wY22_CVID_3.24/FF_allwaves_wY22_CVID_3.24_datafile.RDS")
ffsub<-ff[,c("CVID","cm1intmon","cm1intyr","cm1bsex")]
#ffsub$mintdate<-paste(ff$cm1intmon,ff$cmintyr,sep="/")
ffsub$cm1intmon<-as.character(ffsub$cm1intmon)
ffsub$intmon<-substr(ffsub$cm1intmon,3,6)

ffsub$cm1intyr<-as.character(ffsub$cm1intyr)
ffsub$mintdate<-paste(ffsub$intmon,"1",ffsub$cm1intyr,sep="/")
ffsub$startdate<-as.Date(ffsub$mintdate,format="%b/%d/%Y")
ffsub<-ffsub[,c("CVID","startdate","cm1bsex")]

ffsubm<-merge(corecv_cw,ffsub,by="CVID")
ffsubm_daterec<-merge(ffsubm,daterec,by.x="ffam",by.y="westatID")
ffsubm_daterec$recdate<-as.Date(ffsubm_daterec$`Date Received`,format="%m/%d/%Y")
ffsubm_daterec<-ffsubm_daterec[which(!is.na(ffsubm_daterec$`Lab Number`)),]
ffsubm_daterec$age<-as.numeric(as.character(ffsubm_daterec$recdate-ffsubm_daterec$startdate))/365
# Remove mothers
ffsubm_daterec<-ffsubm_daterec[!ffsubm_daterec$wave %in% "M",]
# fix y and Y
ffsubm_daterec$wave[ffsubm_daterec$wave=="y"]<-"Y"


methidlabid<-readRDS("/drives/drive1/HPA/samples_keys/methid_labid.RDS")
methidlabid<-methidlabid[methidlabid$platform %in% "450K",]
# Since methidlabid didn't match Eve's revert to Colter's fiel
methidlabid2<-read.csv("/drives/drive1/home/lweis/Methylation/finaldatasets_Feb2022/idnumremoved_allclocks_bothstudies_withcells.csv",header=T,stringsAsFactors=F)
methidlabid2<-methidlabid2[methidlabid2$study %in% '450k',]
ffusbm_daterec_methid<-merge(ffsubm_daterec,methidlabid,by.x="Lab Number",by.y="LabID",)
ffusbm_daterec_methidm2<-merge(ffusbm_daterec_methid,methidlabid2,by.x="MethID",by.y="MethID",all.x=T,all.y=T)
# Check to see if CVIDs match
# They don't, although the MethIDs do
checkM<-ffusbm_daterec_methidm2[ffusbm_daterec_methidm2$MethID %in% survey$MethID,]
#survey<-survey[survey$platform %in% "450K",]
checkMo<-checkM[match(survey$MethID,checkM$MethID),]
table(checkMo$MethID==survey$MethID)
table(checkMo$CVID==survey$CVID)
# 450K samples match
# check why so few - Mean-Mean cimt?
cvd<-readRDS("/drives/drive1/ff_chaya/R_versions/FFS_LE8_LS7_all_mergedFFallwvs_Y22sh_CIMT2/FFS_LE8_LS7_all_mergedFFallwvs_Y22sh_CIMT2_datafile.RDS")
cvd<-cvd[,c("CVID","Mean_Mean_CIMT")]

ffusbm_daterec_methidm3<-merge(ffusbm_daterec_methidm2,cvd,by="CVID")
ffusbm_daterec_methidm3<-ffusbm_daterec_methidm3[!ffusbm_daterec_methidm3$dupe %in% "duplicate",]


