library(ComplexHeatmap)
library(showtext)
font_add(family="Arial",regular = "~/.fonts/arial.ttf",
         bold = "~/.fonts/ARIALBD.TTF")
showtext_auto()
setwd("path to sample sheet directory")
all<-readRDS("all_samples_samplesheet.RDS")
samples<-list(`Y9 saliva 450K`=all[all$age=="y9" & all$platform=="450K","ID"],
              `Y9 saliva EPIC`=all[all$age=="y9" & all$platform=="EPIC","ID"],
              `Y15 saliva 450K`=all[all$age=="y15" & all$platform=="450K","ID"],
              `Y15 saliva EPIC`=all[all$age=="y15" & all$platform=="EPIC","ID"],
              `Y22 blood EPIC`=all[all$source=="blood" & all$age=="y22","ID"],
              `Y22_saliva EPIC`=all[all$source=="saliva" & all$age=="y22","ID"]
              )
samplesmat<-make_comb_mat(samples,mode="distinct")
set_size(samplesmat)
UpSet(samplesmat,
      set_order=order(set_size(samplesmat),decreasing=T)
)
cairo_pdf(file="path to upset_plot.pdf",width=10,height=8)
UpSet(samplesmat,
      comb_col=c("dodgerblue","purple","orange","forestgreen","blue")[comb_degree(samplesmat)],
      comb_order=order(comb_size(samplesmat),decreasing=T),
      row_names_gp=gpar(fontface="bold",fontfamily="Arial"),
      set_order=order(set_size(samplesmat),decreasing=T),
      top_annotation = 
        upset_top_annotation(samplesmat,add_numbers=T,height = unit(12,"cm"),
                             numbers_gp=gpar(fontsize=12,fontface="bold",fontfamily="Arial"),
                             axis_param=list(gp=gpar(fontfamily="Arial"))),
      right_annotation = upset_right_annotation(samplesmat,add_numbers=T,
                                                gp=gpar(fill="blue"),
                                                width = unit(8,"cm"),
                                                numbers_gp=gpar(fontsize=12,
                                                                fontface="bold",
                                                                fontfamily="Arial"),
                                                axis_param=list(gp=gpar(fontfamily="Arial"))))
dev.off()

