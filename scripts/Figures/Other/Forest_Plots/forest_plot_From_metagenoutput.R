# Forest plots
library(dplyr)
library(ggplot2)
library(Cairo)
# extract CpG function

library(dplyr)
extract_cpg <- function(meta_obj,cpg_name){
  data.frame(
    CpG=cpg_name,
    study=meta_obj$studlab,
    beta=meta_obj$TE,
    se= meta_obj$seTE
  ) %>%
    mutate(
      lower = beta -1.96 * se,
      upper=beta +1.96*se
    )
}
df_all<-bind_rows(
  lapply(names(metalist),function(cpg){
    extract_cpg(metalist[[cpg]],cpg)
  })
)


# Subset CpGs to plot
cpgs_of_interest<-c("cg07569194","cg08140712","cg08359167","cg09878554","cg16685388")

df_plot<-df_all %>%
  filter(CpG %in% cpgs_of_interest)


  
extract_summary<-function(meta_obj,cpg_name){
    data.frame(
      CpG=cpg_name,
      study="Meta-analysis",
      beta=meta_obj$TE.fixed,
      lower=meta_obj$lower.fixed,
      upper=meta_obj$upper.random)}

df_summary<-bind_rows(
  lapply(names(metalist),function(cpg){
    extract_summary(metalist[[cpg]],cpg)
  })
)
df_summary<-df_summary[df_summary$CpG %in% df_plot$CpG,]
df_plot2<-bind_rows(df_plot,df_summary)
df_plot2$study[df_plot2$study=="epic"]<-"EPIC"
df_plot2$study[df_plot2$study=="450k"]<-"450K"
df_plot2$study<-factor(df_plot2$study,levels=c("Meta-analysis","450K","EPIC"))

#Forest plot - note may need to change xlim
ggplot(df_plot2,aes(x=beta,y=study,shape=study,color=study)) +
  geom_point(size=2) + scale_color_manual(values=c("blue","black","black"))+scale_size_manual(values=c(1.5,1,1))+
  geom_errorbarh(aes(xmin=lower,xmax=upper), height=0.2) +scale_shape_manual(values=c(23,19,19))+
  geom_vline(xintercept=0,linetype="dashed") +
    facet_wrap(~CpG) + 
  xlim(-0.23,0.23)+
  theme_bw()+
xlab("Effect size") + 
       ylab("Cohort")
