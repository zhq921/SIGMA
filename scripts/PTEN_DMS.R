vep.bar.pten2 <- function(dat){
  rownames(dat) <- gsub("_converted", "", rownames(dat))
  rownames(dat)[rownames(dat) == "Combined_score"] <- "SIGMA+"
  dat$Method[dat$Method == "Combined_score"] <- "SIGMA+"
  dat$Correlation <- round(dat$Correlation, digits = 2)
  dat$Method <- factor(rownames(dat), levels = rownames(dat))
  vep_clas <- read.table("./data/other_data/VEPs_classification.tsv", 
                         sep = "\t", header = T, quote = "",
                         stringsAsFactors = F)
  vep_clas$VEPs <- gsub("vep.","",vep_clas$VEPs)
  vep_clas$VEPs <- gsub("_converted","",vep_clas$VEPs)
  vep_clas$VEPs <- gsub("_rankscore","",vep_clas$VEPs)
  dat <- data.frame(dat, 
                    vep_clas[match(rownames(dat), vep_clas$VEPs),c("Label", "Group")])
  p1 <- ggplot(dat[dat$Group == "Individual predictor",], 
               aes(x = Method, y = Correlation)) +
    geom_bar(stat = "identity", position = "identity",
             width = 0.67, fill = "#447fb1", color = "black", size = 0.2) +
    geom_text(aes(label = Correlation), hjust = 0.5, vjust = -1, angle = 0,
              colour = "black", position = position_dodge(.9), size = 2)+
    xlab("Individual predictors")+ylab("Correlation")+
    coord_cartesian(ylim=c(0,0.6))+
    facet_wrap(~ Group, scales = "free", nrow = 2) +
    scale_fill_manual(values = pal_npg()(6)[c(1,2)])+
    theme_bw()+
    theme(axis.text.x=element_text(size=6.5, colour= "black",angle = 45, vjust = 1,hjust = 1),
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
  
  p2 <- ggplot(dat[dat$Group == "meta-predictor",], 
               aes(x = Method, y = Correlation)) +
    geom_bar(stat = "identity", position = "dodge",
             width = 0.67, fill = "#447fb1", color = "black", size = 0.2) +
    geom_text(aes(label = Correlation), hjust = 0.5, vjust = -0.6, angle = 0,
              colour = "black", position = position_dodge(.9), size = 2.5)+
    xlab("Meta-predictors")+ylab("Correlation")+
    coord_cartesian(ylim=c(0,0.5)) +
    facet_wrap(~ Group, scales = "free", nrow = 2) +
    scale_fill_manual(values = pal_npg()(6)[c(1,2)])+
    theme_bw()+
    theme(axis.text.x=element_text(size=6.5, colour= "black",angle = 30, vjust = 1,hjust = 1),
          axis.text.y=element_text(size=8,face="plain",colour= "black"),
          axis.title.y = element_text(size = 10, face="plain"),
          axis.title.x = element_text(size=10, face="plain", colour= "black", margin = margin(0,0,6,0)),
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
  return(list(p.ind = p1, p.meta = p2))
}

# DMS PTEN for another PTEN DMS data
cor4dms4pten2 <- function(dat, snv.only, vep.names, noclinvar = T){
  dat <- dat[dat$Gene %in% c("PTEN"),]
  if(noclinvar){
    dat <- dat[!dat$vep.HGVSp %in% dat.eva$vep.HGVSp,]
  }
  pten2 <- read.csv("./data/DMS/PTEN_matreyek_abundance.csv", stringsAsFactors = F,
                    skip = 4)
  pten2 <- pten2[,c("hgvs_pro", "score")]
  dat.hgvsp <- rownames(dat) %>% strsplit(., split = ":") %>%
    sapply(., function(x) x[2])
  dat$vep.DMS_score <- pten2$score[match(dat.hgvsp, pten2$hgvs_pro)]
  dat <- dat[!is.na(dat$vep.DMS_score),]
  dat$vep.DMS_score <-
    -dat$vep.DMS_score
  
  dat.sub <- dat
  # SIGMA and DMS correlation
  r <- cor(dat.sub$vep.DMS_score, dat.sub$SIGMA, method = "spearman") %>% round(., digits = 3)
  r.p <- cor.test(dat.sub$vep.DMS_score, dat.sub$SIGMA, method = "spearman")$p.value %>% signif(., digits = 3)
  if(r.p < 2.2e-16){
    r.p = 2.2e-16
    dat_text <- data.frame(
      label = paste0("r = ", r, ", P < ", r.p))
  }else{
    dat_text <- data.frame(
      label = paste0("r = ", r, ", P = ", r.p))
  }
  dat.sub$SIGMA <- rank(dat.sub$SIGMA)/nrow(dat.sub)
  dat.sub$vep.DMS_score <- rank(dat.sub$vep.DMS_score)/nrow(dat.sub)
  
  p <- ggplot(dat.sub,aes(x=SIGMA,y=vep.DMS_score))+
    geom_hex(binwidth = 1/30)+
    scale_fill_viridis_c(option = "rocket",direction = -1, name = "Count") +
    ggtitle("PTEN")+
    labs(y="DMS score",x="SIGMA score")+
    coord_cartesian(ylim=c(0,1.05))+
    theme_bw()+
    theme(
      axis.text.x = element_text(color="black",size = 8),
      axis.text.y = element_text(color="black",size = 8),
      plot.title = element_text(color = "black", size = 10, hjust = 0.5),
      axis.title.y = element_text(color="black",size = 10),
      axis.title.x = element_text(color="black",size = 10),
      panel.border = element_blank(),panel.grid=element_blank(),
      axis.line = element_line(colour = "black", size = 0.2),
      axis.ticks = element_line(colour = "black", size = 0.2),
      legend.key.size = unit(8, units ="pt"),
      legend.margin = margin(l=-0.25, unit = "cm"),
      legend.title = element_text(color="black",size = 8),
      legend.text = element_text(color="black",size = 8),
      strip.background = element_blank()
    )+
    geom_text(
      data    = dat_text,
      mapping = aes(x = 0, y = 1.07, label = label),
      size = 2.5,
      hjust   = 0
    )
  assign(paste0("p",1), value = p)
  
  # all VEPs correlation with DMS for each gene
  dat.sub <- dat.sub[,vep.names]
  colnames(dat.sub) <- gsub("vep.","",colnames(dat.sub))
  colnames(dat.sub) <- gsub("_rankscore","",colnames(dat.sub))
  res <- c()
  for(j in 1:ncol(dat.sub)){
    if(colnames(dat.sub)[j] == "DMS_score"){next}
    if(all(is.na(dat.sub[,j]))){
      res <- rbind(res, c(NA, NA))
      rownames(res)[nrow(res)] <- colnames(dat.sub)[j]
      next
    }
    corr <- cor.test(dat.sub[,j], dat.sub$DMS_score, method = "spearman")
    res <- rbind(res, c(corr$estimate, NA))
    rownames(res)[nrow(res)] <- colnames(dat.sub)[j]
  }
  colnames(res) <- c("Correlation", "Correlation_total")
  res <- res[order(res[,1], decreasing = T),]
  res <- data.frame(res, Method = rownames(res), Gene = "PTEN", stringsAsFactors = F)
  corr.all <- res
  return(list(p1 = p1, corr.all = corr.all))
}
