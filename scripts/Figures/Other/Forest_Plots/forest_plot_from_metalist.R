forestplot_metalist_fun<-function(metaoutputlist,group_var,pdftitle) {
  plot_list <- list()
  
  for(i in 1:length(metaoutputlist)){
    tmp<-metaoutputlist[[i]]
    df<-data.frame(cpg=rep(names(metaoutputlist)[i],2),var=tmp$studlab,beta=tmp$TE,se=tmp$seTE)
  df$lower<-df$beta - 1.96*df$se
  df$upper<- df$beta + 1.96*df$se
  
  
  cpg_name<-names(metaoutputlist)[i]
  
  Qp <- signif(tmp$Q, 2)
  I2 <- round(tmp$I2, 1)
  tau2 <- round(tmp$tau2,1)
  
  plot_list[[i]]<-
    ggplot(df,
                         aes(x = beta,
                             y = .data[[group_var]])) +
    
    geom_vline(xintercept = 0,
               linetype = "dashed",
               color = "gray50") +
    
    geom_point(size = 2.5) +
    
    geom_errorbarh(
      aes(xmin = lower,
          xmax = upper),
      height = 0.2
    ) +
    
    theme_bw(base_size = 10) +
    
    labs(
      title = cpg_name,
      x = "Effect estimate",
      y = NULL,
      subtitle = paste0(
        "Q p=", Qp,
        ", I²=", I2, 
        ", Tau²=",tau2,"%"
      )
    ) +
    
    theme(
      plot.title = element_text(size = 10,
                                face = "bold"),
      plot.subtitle = element_text(size = 8)
    )
  names(plot_list)[i]<-cpg_name
}



CairoPDF(pdftitle,
         width = 10,
         height = 8,
         family="Arial")

for(i in seq(1, length(plot_list), by = 6)){

  current_plots <- plot_list[i:min(i+5,
                                   length(plot_list))]

  combined <- wrap_plots(current_plots,
                         ncol = 2)

  print(combined)
}

dev.off()
}
