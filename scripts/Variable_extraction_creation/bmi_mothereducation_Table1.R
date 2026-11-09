setwd("path to sample sheet files")
surv<-readRDS("path to survey data")
bmi<-names(surv)[grep("bmi",names(surv),ignore.case = T)]
surv<-surv[,c("CVID","ck7cflag",bmi)]

files<-list.files(pattern="062726")
filelist<-list()
for(file in 1:length(files)){
  tmp<-readRDS(files[file])
  filelist[[file]]<-merge(tmp,surv,by="CVID",all.x=T,all.y=F)
}
for(i in 1:length(filelist)){
  print(files[i])
  print(mean(filelist[[i]]$eduMom_yrs,na.rm=T))
  print(range(filelist[[i]]$eduMom_yrs,na.rm=T))
  print(mean(filelist[[i]]$anth_bmi_recl,na.rm=T))
  print(range(filelist[[i]]$anth_bmi_recl,na.rm=T))
}
y9<-rbind(filelist[[4]],filelist[[5]])
y15<-rbind(filelist[[1]],filelist[[2]])
y22<-filelist[[3]]

filelist2<-list(
  y9=rbind(filelist[[4]],filelist[[5]]),
  y15=rbind(filelist[[1]],filelist[[2]]),
  y22=filelist[[3]]
)
for (i in 1:length(filelist2)){
  print(names(filelist2)[i])
  print(mean(filelist2[[i]]$eduMom_yrs,na.rm=T))
  print(range(filelist2[[i]]$eduMom_yrs,na.rm=T))
  print(mean(filelist2[[i]]$anth_bmi_recl,na.rm=T))
  print(range(filelist2[[i]]$anth_bmi_recl,na.rm=T))
}
