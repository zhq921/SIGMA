# mutational constraint
constraint_stra <- function(dat, constraintFile){
  dat.cons <- fread(constraintFile, data.table = F, header = T)
  dat.cons <- dat.cons[,c("gene", "mis_z", "n_mis", "exp_mis")]
  dat <- dat[,c("alleleID", "label", "SIGMA", "vep.SYMBOL")]
  dat <- data.frame(
    dat,
    mis_z = dat.cons$mis_z[match(dat$vep.SYMBOL, dat.cons$gene)],
    stringsAsFactors = F
  ) 
  dat <- na.omit(dat) 
  dat <- data.frame(
    dat,
    group = dat$mis_z %>% cut(., breaks = c(-10, seq(-3, 6, by = 1), 20))
  )
  dat2 <- dat # all
  dat2$label <- "All"
  dat.plot <- rbind(dat, dat2) # benign + pathogenic + all
  dat.plot$label <- as.factor(dat.plot$label)
  stat_wilcox <- wilcox_test(group_by(dat.plot, group), SIGMA~label) 
  stat_wilcox <- add_significance(stat_wilcox, 'p')  
   stat_wilcox.test <-  add_xy_position(stat_wilcox, x = 'group')
  stat_wilcox.test <- stat_wilcox.test[stat_wilcox.test$group1!="All",]
   # number
  nn <- table(dat$group)
  nn <- lapply(nn, function(x){
    x <- as.character(x)
    paste0(substr(x, start = 1, stop = nchar(x)-3),
           ",",
           substr(x, start = nchar(x)-2, stop = nchar(x)))
  }) %>% unlist()
  
  p <- ggboxplot(dat.plot, x = 'group', y = 'SIGMA', fill = 'label', 
                 color = 'black', width = 0.7, size = 0.2, legend = 'none',
                 outlier.shape = NA) +
    scale_fill_manual(values = pal_npg(alpha = 0.85)(4)[c(4,3,1)], name = "Label") +
    labs(x = 'Mutational Constraint (Z-score)', y = 'SIGMA score') +
    ylim(c(0,1.02))+
    scale_x_discrete(labels = 
                       paste(dat$group %>% levels(),
                             "\n(",
                             nn, ")", sep = "")
    )+
    theme_bw()+
    theme(
      axis.text.x = element_text(color="black",size = 8,angle = 60, 
                                 hjust = 1, vjust = 1),
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
