library(ggplot2)
library(Cairo)
library(patchwork)
library(showtext)
font_add(family="Arial",regular = "~/.fonts/arial.ttf",
         bold = "~/.fonts/ARIALBD.TTF")
showtext_auto()

setwd("output directory filepath")

y9promfisher<-read.delim("path to y9_DMRs_DMRcate_combp_promoter_fishertest.table.txt")
y15promfisher<-read.delim("path to y15_DMRs_DMRcate_combp_promoter_fishertest.table.txt")
y22promfisher<-read.delim("path to y22_DMRs_DMRcate_combp_promoter_fishertest.table.txt")


y9promfisher$age<-"Age 9"
y15promfisher$age<-"Age 15"
y22promfisher$age<-"Age 22"

promfisher<-rbind(y9promfisher, y15promfisher, y22promfisher)
promfisher$Row.names<-factor(promfisher$Row.names,
                               levels=c("Distal Intergenic","Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5' UTR","Exon","Intron","3' UTR"),
                               labels=c("Distal Intergenic","Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5' UTR","Exon","Intron","3' UTR"))
promfisher$Age<-factor(promfisher$age, levels=c("Age 9","Age 15", "Age 22"))

names(promfisher)[names(promfisher)=="Row.names"]<-"Feature"
# add 0.01 to zero values
promfisher$label_y_pos<-promfisher$OR + 0.1

##### CpG enrichment
y9fisher<-read.delim("path to y9_DMRs_DMRcate_combp_relationtoCpGIsland_fishertest.table.txt")
y15fisher<-read.delim("path to y15_DMRs_DMRcate_combp_relationtoCpGIsland_fishertest.table.txt")
y22fisher<-read.delim("path to y22_DMRs_DMRcate_combp_relationtoCpGIsland_fishertest.table.txt")
y9fisher$age<-"Age 9"
y15fisher$age<-"Age 15"
y22fisher$age<-"Age 22"

cpg<-rbind(y9fisher, y15fisher, y22fisher)


cpg$Feature<-factor(cpg$Row.names,
                      levels=c("Shelf","Shore","CpG Island","Open Sea"),
                      labels=c("Shelf","Shore","CpG Island","Open Sea"))
cpg$Age<-factor(cpg$age,levels=c("Age 9","Age 15","Age 22"))

# Make plots
promcolor<-c("Promoter (2-3kb)"="#A6CEE3","Promoter (1-2kb)"="#1F78B4","Promoter (<=1kb)"="#B2DF8A","5' UTR"="#33A02C","Exon"="#FB9A99","Intron"="#E31A1C","3' UTR"="#FDbF6F","Distal Intergenic"="#FF7F00")

p3<-ggplot(promfisher,aes(x=Age,y=OR,fill=Feature))+geom_hline(yintercept = 1,linetype=2,color="gray")+geom_col(position=position_dodge(width=0.9))+
  scale_fill_manual(values=promcolor)+theme_bw(base_family = "Arial")+ggtitle("DMRs in Relation to CpG Islands")+
  ggtitle("DMRs in Genomic Features")+
  ylab("Odds Ratio")+
  theme(axis.title.x=element_blank(), 
         #plot.margin = margin(5,2,5,2,unit ="pt"),
         axis.text.x = element_text(face="bold",size=10,family="Arial",hjust=0.5),
         axis.title=element_text(face="bold",family="Arial",size=10),
         plot.title=element_text(face="bold",family="Arial",size=9,hjust=0.5),
         panel.grid=element_blank(),
         legend.position="bottom",
         legend.direction="horizontal",
         legend.text = element_text(family="Arial",size=10),
         #legend.title=element_text(family="Arial",face="bold",size=7),
         legend.title=element_blank(),
         legend.key.height=unit(3,"mm"),
         legend.key.width=unit(3,"mm"),
         legend.box.spacing=unit(1,"mm"),
         legend.key.spacing.x=unit(1,"mm"),
         legend.key.spacing.y = unit(1,"mm")
  ) +
  guides(
    fill=guide_legend(nrow=4)
  )+
  geom_text(aes(y=label_y_pos,label=sig),color="#36454F",position=position_dodge(width=0.9), vjust=0, size=4,family="Arial")


cpgcolor<-c("Shelf"="#A6CEE3","Shore"="#B2DF8A","CpG Island"="#FF7F00","Open Sea"="#1F78B4")

p4<-ggplot(cpg,aes(x=Age,y=OR,fill=Feature))+geom_hline(yintercept = 1,linetype=2,color="gray")+
  geom_col(position=position_dodge(width=0.9))+
  scale_fill_manual(values=cpgcolor)+
  theme_bw(base_family = "Arial")+
  ggtitle("DMRs in Relation to CpG Islands")+
  ylim(0,5)+
  ylab("Odds Ratio")+
  theme(axis.title.x=element_blank(), 
        #plot.margin = margin(5,5,5,5,unit ="pt"),
        axis.text.x = element_text(face="bold",size=10,family="Arial",hjust=1),
        axis.title=element_text(face="bold",family="Arial",size=10),
        plot.title=element_text(face="bold",family="Arial",size=9,hjust=0.5),
        panel.grid=element_blank(),
        legend.position="bottom",
        legend.direction="horizontal",
        legend.text = element_text(family="Arial",size=10),
        #legend.title=element_text(family="Arial",face="bold",size=7),
        legend.title=element_blank(),
        legend.key.height=unit(3,"mm"),
        legend.key.width=unit(3,"mm"),
        legend.box.spacing=unit(1,"mm"),
        legend.key.spacing.x=unit(1,"mm"),
        legend.key.spacing.y = unit(1,"mm")
  ) +
  
  
  
  guides(
    fill=guide_legend(nrow=3)
  )+
  geom_text(aes(y=label_y_pos,label=sig),color="#36454F",position=position_dodge(width=0.9), vjust=0, size=4,family="Arial")+ ylab("Odds Ratio")


CairoPDF(file="Fig4d_e.pdf",width = 7,height=3.5)
tag_theme<-theme(plot.tag=element_text(
  face="bold",family="Arial",size=14))

p<-(p3+tag_theme)|
  (p4+tag_theme)
p+
  plot_annotation(tag_levels=list(c('d','e')
  ))


dev.off()
