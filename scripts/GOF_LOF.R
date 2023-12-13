# GOF LOF
gof_lof_stra <- function(dat, cutoff){
  dat <- dat[dat$label == "Pathogenic",]
  goflof.info <- read.csv("./data/other_data/goflof_ClinVar_v062021.csv", 
                          header = T, stringsAsFactors = F)
  goflof.info$AlleleID <- as.character(goflof.info$AlleleID) # 总共5085个突变，555 GOF and 4,530 LOF variants
  goflof.info <- goflof.info[goflof.info$AlleleID %in% dat$alleleID,] # 是否在clinvar中，193 GOF, 940 LOF，其余主要在有3000是非SNV，还有一些是没星
  all(goflof.info$AlleleID %in% dat$alleleID)
  dat <- dat[match(goflof.info$AlleleID, dat$alleleID),]
  dat.info <- data.frame(
    SIGMA = dat$SIGMA,
    label = dat$label,
    gofLof = goflof.info$LABEL,
    SIGMA_pred = ifelse(dat$SIGMA > cutoff, "Pathogenic", "Benign")
  )
  test <- table(dat.info$gofLof, dat.info$SIGMA_pred) %>% chisq.test()
  dat.plot <- table(dat.info$gofLof, dat.info$SIGMA_pred) %>% reshape2::melt()
  lof.num <- table(dat.info$gofLof)[2]
  gof.num <- table(dat.info$gofLof)[1]
  dat.plot <- data.frame(
    dat.plot,
    perc = (100*dat.plot$value / c(gof.num,lof.num,gof.num,lof.num)) %>% round(.,digits = 2)
  )
  dat.plot <- data.frame(
    dat.plot,
    label_y = paste0(dat.plot$perc, "%"),
    label_ypos = c(100, 100, dat.plot$perc[3:4])
  )
  colnames(dat.plot)[2] <- "SIGMA"
  col <- c("#377eb8", "#e41a1c") # pal_npg()(2) %>% rev()
  p <- ggplot(data = dat.plot, aes(x = Var1, y = perc, fill = SIGMA)) +
    geom_bar( stat="identity", width = 0.6, color = "black", size = 0.2)+
    geom_text(aes(y=label_ypos, label=label_y), vjust=1.2, 
              color="black", size=3)+
    geom_text(x = 1.5, y = 104, label = 
                paste("P =",signif(test$p.value, digits = 3)),
              size = 3, color = "white") +
    xlab(label = NULL) +
    ylab("Percentage (%)") +
    ylim(c(0,105))+
    scale_fill_manual(name = "Prediction", values = col)+
    theme_classic()+
    theme(#legend.position = "none",
      title = element_text(size = 4, colour = "black"),
      legend.title = element_text(size = 8, colour = "black"),
      axis.text.x = element_text(size = 8, colour = "black"),
      axis.text.y = element_text(size = 8, colour = "black"),
      axis.ticks.x = element_line(colour = "black", size = 0.32),
      axis.ticks.y = element_line(colour = "black", size = 0.32),
      axis.line = element_line(colour = "black", size = 0.32),
      axis.title.y = element_text(size = 10, colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(colour = "white", fill = "white"),
      strip.text = element_text(colour = "black", size = 6)
    )
  return(p)
}
