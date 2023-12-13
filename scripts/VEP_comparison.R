
# vep comparison; return performance and correlation
vep_comparison <- function(dat, vep.names){
  dat.indv <- dat
  dat.meta <- dat[!is.na(dat$Combined_score),]
  vep_clas <- read.table("./data/other_data/VEPs_classification.tsv", 
                         sep = "\t", header = T, quote = "",
                         stringsAsFactors = F)
  
  res <- matrix(NA, nrow = length(vep.names), ncol = 4)
  colnames(res) <- c("AUC", "lower", "upper", "p_value")
  rownames(res) <- vep.names
  res <- data.frame(res)
  
  for(i in 1:length(vep.names)){
    if(vep.names[i] == "Combined_score"){
      dat <- dat.meta
      meta <- T
    }else{
      if(vep_clas$Group[vep_clas$VEPs==vep.names[i]] %in% "meta-predictor" ){
        dat <- dat.meta
        meta <- T
      }else{
        dat <- dat.indv
        meta <- F
      }
    }
    roc.vep <- pROC::roc(response = dat$label, predictor = dat[,vep.names[i]], auc = T, ci = T, plot = F)
    auc.vep <- round(roc.vep$auc,3)
    lower.vep <- round(roc.vep$ci[1],3)
    upper.vep <- round(roc.vep$ci[3],3)
    if(meta){
      p.vep <- roc.test(response = dat$label, 
                        predictor1 = dat$Combined_score, 
                        predictor2 = dat[,vep.names[i]]
      )$p.value
    }else{
      p.vep <- roc.test(response = dat$label, 
                        predictor1 = dat$SIGMA, 
                        predictor2 = dat[,vep.names[i]])$p.value
    }
    res[vep.names[i],] <- c(auc.vep, lower.vep, upper.vep, p.vep)
  }
  res <- res[order(res$AUC, decreasing = T),]
  rownames(res) <- gsub("vep.","",rownames(res))
  rownames(res) <- gsub("_rankscore","",rownames(res))
  rownames(res)[rownames(res) == "Combined_score"] <- "SIGMA+"
  
  dat.sub <- dat[,vep.names]
  colnames(dat.sub) <- gsub("vep.","",colnames(dat.sub))
  colnames(dat.sub) <- gsub("_rankscore","",colnames(dat.sub))
  colnames(dat.sub) <- gsub("_converted","",colnames(dat.sub))
  corr <- cor(dat.sub)
  for(i in 1:length(vep.names)){
    for(j in 1:length(vep.names)){
      if(j<i){next}
      tmp <- dat.sub[,c(i,j)]
      tmp <- na.omit(tmp)
      tmp.cor <- cor(tmp[,1], tmp[,2], method = "spearman")
      corr[i,j] <- tmp.cor
      corr[j,i] <- tmp.cor
    }
  }
  rownames(corr)[rownames(corr) == "Combined_score"] <- "SIGMA+"
  colnames(corr)[colnames(corr) == "Combined_score"] <- "SIGMA+"
  return(list(vep.perf = res, vep.corr = corr))
}

## performance AUC
vep.bar <- function(dat){
  rownames(dat) <- gsub("_converted", "", rownames(dat))
  dat$Method <- factor(rownames(dat), levels = rownames(dat))
  vep_clas <- read.table("./data/other_data/VEPs_classification.tsv", 
                         sep = "\t", header = T, quote = "",
                         stringsAsFactors = F)
  vep_clas$VEPs <- gsub("vep.","",vep_clas$VEPs)
  vep_clas$VEPs <- gsub("_converted","",vep_clas$VEPs)
  vep_clas$VEPs <- gsub("_rankscore","",vep_clas$VEPs)
  dat <- data.frame(dat, 
                    vep_clas[match(rownames(dat), vep_clas$VEPs),c("Label", "Group")])
  # add * **
  p.star <- dat[dat$Group == "Individual predictor","p_value"]
  p.star <- lapply(p.star, function(pValue){
    if(pValue < 0.0001){return("****")}
    if(pValue < 0.001){return("***")}
    if(pValue < 0.01){return("**")}
    if(pValue < 0.05){return("*")}
    if(pValue >= 0.05){return("")}
  }) %>% unlist
  
  p1 <- ggplot(dat[dat$Group == "Individual predictor",], aes(x = Method, y = AUC)) +
    geom_bar(stat = "identity", position = "identity",
             width = 0.67, fill = "#447fb1", color = "black", size = 0.2) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,size = 0.1,
                  position=position_dodge(.9))+ 
    geom_text(aes(label = AUC), hjust = 1.2, angle = 90,
              colour = "white", position = position_dodge(.9), size = 2)+
    geom_text(label = p.star, hjust = 0.5, angle = 0, vjust = -0.5,
              colour = "black", position = position_dodge(.9), size = 2)+
    xlab("Variant effect predictor")+ylab("AUC")+
    coord_cartesian(ylim=c(0.5,1))+
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
  
  # add * **
  p.star <- dat[dat$Group == "meta-predictor","p_value"]
  p.star <- lapply(p.star, function(pValue){
    if(pValue < 0.0001){return("****")}
    if(pValue < 0.001){return("***")}
    if(pValue < 0.01){return("**")}
    if(pValue < 0.05){return("*")}
    if(pValue >= 0.05){return("")}
  }) %>% unlist
  p2 <- ggplot(dat[dat$Group == "meta-predictor",], aes(x = Method, y = AUC)) +
    geom_bar(stat = "identity", position = "dodge",
             width = 0.67, fill = "#447fb1", color = "black", size = 0.2) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,size = 0.1,
                  position=position_dodge(.9))+ 
    geom_text(aes(label = AUC), hjust = 1.35, angle = 90,
              colour = "white", position = position_dodge(.9), size = 2.5)+
    geom_text(aes(label = p.star), hjust = 0.5, angle = 0, vjust = -1,
              colour = "black", position = position_dodge(.9), size = 2)+
    xlab("Variant effect predictor")+ylab("AUC")+
    coord_cartesian(ylim=c(0.5,1)) +
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
