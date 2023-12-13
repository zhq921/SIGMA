
# score distribution
score_distri_plot <- function(dat){
  dat.bar <- data.frame(
    score_bin = paste(seq(0,0.95, by = 0.05), seq(0.05, 1, by = 0.05), sep ="-"),
    Group = c(rep("Benign",20), rep("Pathogenic", 20)),
    count = c(cut(dat$SIGMA[dat$label=="Benign"], breaks = seq(0,1,length.out = 21)) %>% table,
              cut(dat$SIGMA[dat$label=="Pathogenic"], breaks = seq(0,1,length.out = 21)) %>% table),
    freq = c((cut(dat$SIGMA[dat$label=="Benign"], breaks = seq(0,1,length.out = 21)) %>% table)*100 / table(dat$label)[1],
             (cut(dat$SIGMA[dat$label=="Pathogenic"], breaks = seq(0,1,length.out = 21)) %>% table)*100 / table(dat$label)[2])
  )
  p.density <- ggplot(dat, aes(x = SIGMA, fill = label)) + 
    geom_density(size = 0.3) +
    scale_fill_manual(values = c("#377eb8dd", "#e41a1cdd")) +
    theme_bw()+
    xlab("SIGMA score") + 
    ylab("Density") +
    theme(
      axis.text.x = element_text(color="black",size = 8),
      axis.text.y = element_text(color="black",size = 8),
      axis.title.x = element_text(color="black",size = 10),
      axis.title.y = element_text(color="black",size = 10),
      legend.text = element_text(color="black",size =9),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.size = unit(0.2, "cm"),
      panel.border = element_blank(),panel.grid=element_blank(),
      axis.line = element_line(colour = "black", size = 0.3),
      axis.ticks = element_line(colour = "black", size = 0.3)
    )
  p.hist <- ggplot(dat.bar, aes(x = score_bin,y = freq,fill = Group))+
    geom_bar(stat ="identity",width = 0.7,position ="dodge", color = "black", size = 0.2)+
    scale_fill_manual(values=c("#377eb8", "#e41a1c"))+
    labs(x = "SIGMA score bin", y = "Percentage (%)")+
    theme_bw()+
    theme(
      axis.text.x = element_text(color="black",size = 7, angle = 45, vjust = 1,hjust = 1),
      axis.text.y = element_text(color="black",size = 8, margin = margin(r = 3, l =1)),
      axis.title.x = element_text(color="black",size = 10),
      axis.title.y = element_text(color="black",size = 10),
      legend.text = element_text(color="black",size = 8),
      legend.title =element_blank(),
      legend.position = c(0.3,0.75),
      legend.background = element_rect(fill = "white", color = "black",
                                       size = 0.3),
      legend.margin = margin(t= 0, b = 0.2,r = 0.2, l = 0.2, unit = "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm"),
      panel.border = element_blank(),panel.grid=element_blank(),
      axis.line = element_line(colour = "black", size = 0.3),
      axis.ticks = element_line(colour = "black", size = 0.3)
    )
  return(list(p.density = p.density, p.hist = p.hist))
}
# score distribution sigmap
score_distri_plot_sigmap <- function(dat){
  colnames(dat)[colnames(dat) == "Combined_score"] <- "SIGMAp"
  dat <- dat[!is.na(dat$SIGMAp),]
  dat.bar <- data.frame(
    score_bin = paste(seq(0,0.95, by = 0.05), seq(0.05, 1, by = 0.05), sep ="-"),
    Group = c(rep("Benign",20), rep("Pathogenic", 20)),
    count = c(cut(dat$SIGMAp[dat$label=="Benign"], breaks = seq(0,1,length.out = 21)) %>% table,
              cut(dat$SIGMAp[dat$label=="Pathogenic"], breaks = seq(0,1,length.out = 21)) %>% table),
    freq = c((cut(dat$SIGMAp[dat$label=="Benign"], breaks = seq(0,1,length.out = 21)) %>% table)*100 / table(dat$label)[1],
             (cut(dat$SIGMAp[dat$label=="Pathogenic"], breaks = seq(0,1,length.out = 21)) %>% table)*100 / table(dat$label)[2])
  )
  p.density <- ggplot(dat, aes(x = SIGMAp, fill = label)) + 
    geom_density(size = 0.3) +
    scale_fill_manual(values = c("#377eb8dd", "#e41a1cdd")) +
    theme_bw()+
    xlab("SIGMA+ score") + 
    ylab("Density") +
    theme(
      axis.text.x = element_text(color="black",size = 8),
      axis.text.y = element_text(color="black",size = 8),
      axis.title.x = element_text(color="black",size = 10),
      axis.title.y = element_text(color="black",size = 10),
      legend.text = element_text(color="black",size =9),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.size = unit(0.2, "cm"),
      panel.border = element_blank(),panel.grid=element_blank(),
      axis.line = element_line(colour = "black", size = 0.3),
      axis.ticks = element_line(colour = "black", size = 0.3)
    )
  p.hist <- ggplot(dat.bar, aes(x = score_bin,y = freq,fill = Group))+
    geom_bar(stat ="identity",width = 0.7,position ="dodge", color = "black", size = 0.2)+
    scale_fill_manual(values=c("#377eb8", "#e41a1c"))+
    labs(x = "SIGMA+ score bin", y = "Percentage (%)")+
    theme_bw()+
    theme(
      axis.text.x = element_text(color="black",size = 8, angle = 45, vjust = 1,hjust = 1),
      axis.text.y = element_text(color="black",size = 8),
      axis.title.x = element_text(color="black",size = 10),
      axis.title.y = element_text(color="black",size = 10),
      legend.text = element_text(color="black",size =9),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.size = unit(0.2, "cm"),
      panel.border = element_blank(),panel.grid=element_blank(),
      axis.line = element_line(colour = "black", size = 0.3),
      axis.ticks = element_line(colour = "black", size = 0.3)
    )
  return(list(p.density = p.density, p.hist = p.hist))
}
