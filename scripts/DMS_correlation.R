# DMS
cor4dms <- function(dat, snv.only, vep.names, noclinvar = T){
  dat <- dat[!dat$Gene %in% c("SNCA"),]
  dat$Gene[dat$Gene == "TP53"] <- "P53"
  dat <- dat[!is.na(dat$vep.DMS_score),]
  if(snv.only){dat <- dat[grepl(">", dat$vep.HGVSc),]}
  if(noclinvar){
    dat <- dat[!dat$vep.HGVSp %in% dat.eva$vep.HGVSp,]
  }
  genes <- dat$Gene %>% unique() %>% as.character()
  for(i in 1:length(genes)){
    dat$vep.DMS_score[dat$Gene == genes[i]] <- 
      scale(dat$vep.DMS_score[dat$Gene == genes[i]])[,1]
    if(genes[i] %in% c("PTEN", "HRAS", "KCNH2", "VKORC1")){
      dat$vep.DMS_score[dat$Gene == genes[i]] <-
        -dat$vep.DMS_score[dat$Gene == genes[i]]}
  }
  
  corr.all <- c()
  for(i in 1:length(genes)){
    dat.sub <- dat[dat$Gene == genes[i],]
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
      ggtitle(genes[i])+
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
    assign(paste0("p",i), value = p)
    
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
    res <- data.frame(res, Method = rownames(res), Gene = genes[i], stringsAsFactors = F)
    corr.all <- rbind(corr.all, res)
  }
  # #  all VEPs correlation with DMS
  dat <- dat[,vep.names]
  colnames(dat) <- gsub("vep.","",colnames(dat))
  colnames(dat) <- gsub("_rankscore","",colnames(dat))
  res <- c()
  for(i in 1:ncol(dat)){
    if(colnames(dat)[i] == "DMS_score"){next}
    corr.mean <- corr.all$Correlation[corr.all$Method == colnames(dat)[i]] %>% mean(., na.rm = T)
    corr.total <- cor(dat$DMS_score, dat[,colnames(dat)[i]], use = "complete.obs",
                      method = "spearman") %>% round(., digits = 3)
    res <- rbind(res, c(corr.mean,corr.total))
    rownames(res)[nrow(res)] <- colnames(dat)[i]
  }
  colnames(res) <- c("Correlation", "Correlation_total")
  res <- res[order(res[,1], decreasing = T),]
  res <- data.frame(res, Method = rownames(res), Gene = "All", stringsAsFactors = F)
  corr.all <- rbind(res,corr.all)
  return(list(p1 = p1,p2 = p2,p3 = p3,p4 = p4,p5 = p5, 
              p6 = p6, corr.all = corr.all))
}
