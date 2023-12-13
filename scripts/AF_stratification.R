af_stra <- function(dat = dat.eva){
  dat <- data.frame(
    af = dat$vep.gnomAD_genomes_POPMAX_AF,
    score = dat$SIGMA,
    label = dat$label,
    stringsAsFactors = F
  )
  dat$af[is.na(dat$af)] <- 0
  dat <- data.frame(
    dat,
    af_group = cut(dat$af, breaks = c(-1,0, 0.0001, 0.001, 0.01, 0.1,1))
  )
  dat2 <- dat # all
  dat2$label <- "All"
  dat.plot <- rbind(dat, dat2) # benign + pathogenic + all
  dat.plot$label <- as.factor(dat.plot$label)
  stat_wilcox <- wilcox_test(group_by(dat.plot, af_group), score~label)
  stat_wilcox <- add_significance(stat_wilcox, 'p') 
   stat_wilcox.test <-  add_xy_position(stat_wilcox, x = 'af_group')
  stat_wilcox.test <- stat_wilcox.test[stat_wilcox.test$group1!="All",]
  # number
  nn <- table(dat$af_group)
  nn <- lapply(nn, function(x){
    x <- as.character(x)
    paste0(substr(x, start = 1, stop = nchar(x)-3),
           ",",
           substr(x, start = nchar(x)-2, stop = nchar(x)))
  }) %>% unlist()
  
  p <- ggboxplot(dat.plot, x = 'af_group', y = 'score', fill = 'label', 
                 color = 'black', width = 0.6, size = 0.2, legend = 'none',
                 outlier.shape = NA) +
    scale_fill_manual(values = pal_npg(alpha = 0.85)(4)[c(4,3,1)], name = "Label") +
    labs(x = 'Allele Frequency', y = 'SIGMA score') +
    ylim(c(0,1.02))+
    scale_x_discrete(labels = 
                       paste(c("0%", "0-0.01%", "0.01-0.1%",
                               "0.1-1%", "1-10%", ">10%"),
                             "\n(",
                             nn, ")", sep = "")
    ) +
    theme_bw()+
    theme(
      legend.position = "none",
      axis.text.x = element_text(color="black",size = 8,
                                 angle = 60, hjust = 1, vjust = 1),
      axis.text.y = element_text(color="black",size = 8),
      axis.title.y = element_text(color="black",size = 10),
      axis.title.x = element_text(color="black",size = 10),
      title = element_text(color = "black", size = 10),
      panel.border = element_blank(),panel.grid=element_blank(),
      axis.line = element_line(colour = "black", size = 0.2),
      axis.ticks = element_line(colour = "black", size = 0.2),
    )+ 
    stat_pvalue_manual(stat_wilcox.test, label = 'p.signif', tip.length = 0, y.position = 1.02,
                       bracket.size = 0.2,
                       label.size = 2.5)
  return(p)
}