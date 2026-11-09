
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(ggsignif)

setwd("path to working directory")
y22 <- readRDS("path to Age 22 blood sample sheet")
y22$Race<-ifelse(y22$raceEthY22=="BLACK","NHB","ALL OTHER")

library(Cairo)
CairoPDF("path to save CIMT_density.pdf",family="Arial")
ggplot(y22,aes(x=Mean_Mean_CIMT))+geom_density(aes(fill=Race),alpha=0.3)+scale_fill_manual(values=c("yellow","dodgerblue"))+theme_bw()+xlab("Mean Mean CIMT")+ylab("Density")+theme(legend.title = element_blank(),axis.title=element_text(face="bold",family="Arial"))
dev.off()
df<-y22
df$ethrace<-factor(df$raceEthY22,
levels = c("BLACK",
"HISPANIC",
"OTHER","WHITE"))

df$ethrace2<-factor(df$raceEthY22,levels = c("BLACK","WHITE",
                              "HISPANIC",
                              "ALL OTHER"))

library(FSA)
kruskal.test(df$Mean_Mean_CIMT ~ df$ethrace2)
dunn_out<-dunnTest(df$Mean_Mean_CIMT~df$ethrace2,method="bonferroni")
dunn_df<-dunn_out$res
dunn_df$group1<-gsub(" - [A-Z]*","",dunn_df$Comparison)
dunn_df$group2<-gsub("[A-Z]* - ","",dunn_df$Comparison)

dunn_df$p.adj.significant <- cut(
  dunn_df$P.adj,
  breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
  labels = c("***", "**", "*", ".", "ns")
  
)
library(Cairo)
CairoPDF("path to save Fig2a: cimt_by_ethrace.pdf",family="Arial",width=6,height=6)
ggplot(df, aes(x = ethrace2, y = Mean_Mean_CIMT, color = ethrace)) +
  geom_boxplot() +
  labs(x = "Self-reported Race", y = "Mean-mean CIMT") + 
  # scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  ggtitle("CIMT by Self-reported Race") +
  scale_color_manual(
    labels = c("WHITE" = "White",
               "BLACK" = "Black",
               "HISPANIC" = "Hispanic", 
               "OTHER" = "Other"),
    values = c("WHITE" = "#66C2A5", 
               "BLACK" = "#FC8D62",
               "HISPANIC" = "#8DA0CB",  
               "OTHER" = "#E78AC3")
  ) +
    scale_x_discrete(
      labels = c("WHITE" = "White
  n = 200",
                 "BLACK" = "Black
  n = 623",
                 "HISPANIC" = "Hispanic
  n = 354",
                 "OTHER" = "Other
  n = 55")
  ) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  stat_pvalue_manual(dunn_df[dunn_df$p.adj.significant != "ns", ], # not working. worked when using color and not fill
                     label = "p.adj.significant", 
                     xmin = "group1", xmax = "group2", 
                     y.position = c(0.75, 0.80))
dev.off()

# By sex
df$sex<-as.character(df$cm1bsex)
df$sex<-factor(df$sex,levels=c("1 Boy","2 Girl"))

wilcox.test(df$Mean_Mean_CIMT ~ df$sex)



CairoPDF("path to save Fig 2b: cimt_by_sex.pdf",family="Arial",height=6,width=6) 
ggplot(df, aes(x = sex, y = Mean_Mean_CIMT, fill =sex)) +
  geom_boxplot() +
  labs(x = "Sex", y = "Mean-mean CIMT") + 
  # scale_color_brewer(palette = "Set2") +
  theme_classic() +
  ggtitle("CIMT by Sex") +
   scale_fill_manual(
     labels = c("1 Boy" = "Male",
                "2 Girl" = "Female"),
     values = c("lightcyan3", "thistle")
   ) +
  scale_x_discrete(
    labels = c("1 Boy" = "Male
n = 549",
               "2 Girl" = "Female
n = 683")
  ) +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  geom_signif(comparisons = list(c("1 Boy", "2 Girl")), map_signif_level = T,
              test = "wilcox.test", color = "black") #+
 dev.off()

