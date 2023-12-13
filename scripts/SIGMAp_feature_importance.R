var_imp_plot_sigmap <- function(dat){
  dat.plot <- dat
  dat.plot$variable <- factor(dat.plot$variable, levels = rev(dat.plot$variable))
  p<-ggplot(data = dat.plot, aes(x = variable, y = scaled_importance))+
    geom_bar(stat = "identity",width = 0.67, fill = "#447fb1", color = "black", size = 0.2) +
    ylab("Scaled importance")+ xlab("Predictors")+
    coord_flip()+
    theme_bw()+ 
    theme(axis.text.x=element_text(size=8, colour= "black"), 
          axis.text.y=element_text(size=8,face="plain",colour= "black"), 
          axis.title.y=element_text(size = 10,face="plain"), 
          axis.title.x=element_text(size = 10,face="plain"),
          axis.line = element_line(colour = "black", size=0.3), 
          axis.ticks = element_line(colour = "black", size = 0.2),
          legend.text=element_text(colour="black",  
                                   size=10),
          legend.title=element_blank(),
          legend.position = c(0.8,0.2),
          legend.key.size = unit(0.5,"cm"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) 
  return(p)
}
