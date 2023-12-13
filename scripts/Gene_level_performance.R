gene_stra <- function(dat, cutoff, cutoff.eve = 0.5, cutoff.revel = 0.5){
  dat <- dat[,c("label", "SIGMA", "vep.SYMBOL", "vep.EVE", "vep.REVEL_score")]
  colnames(dat)[colnames(dat) == "SIGMA"] <- "SIGMA"
  table(dat$vep.SYMBOL, dat$label)
   genes <- dat$vep.SYMBOL %>% unique
  dat <- dat[dat$vep.SYMBOL %in% genes,]
  gene.info <- data.frame(
    gene = genes,
    vrt.num = NA,
    b.num = NA,
    p.num = NA,
    TP = NA,
    TN = NA,
    SEN = NA,
    SPE = NA,
    ACC = NA,
    AUC = NA,
    vrt.num.EVE = NA,
    b.num.EVE = NA,
    p.num.EVE = NA,
    TP.EVE = NA,
    TN.EVE = NA,
    SEN.EVE = NA,
    SPE.EVE = NA,
    ACC.EVE = NA,
    AUC.EVE = NA,
    vrt.num.REVEL = NA,
    b.num.REVEL = NA,
    p.num.REVEL = NA,
    TP.REVEL = NA,
    TN.REVEL = NA,
    SEN.REVEL = NA,
    SPE.REVEL = NA,
    ACC.REVEL = NA,
    AUC.REVEL = NA,
    stringsAsFactors = F
  )
  for(i in 1:length(genes)){
    dat.sub <- dat[dat$vep.SYMBOL %in% genes[i],] 
    
    vrt.num <- nrow(dat.sub)
    b.num <- sum(dat.sub$label=="Benign")
    p.num <- sum(dat.sub$label=="Pathogenic")
    TP <- sum((dat.sub$SIGMA > cutoff) & (dat.sub$label=="Pathogenic"))
    TN <- sum((dat.sub$SIGMA <= cutoff) & (dat.sub$label=="Benign"))
    if(p.num == 0){
      SEN <- NA
    }else{
      SEN <- (TP/p.num) %>% round(., digits = 3)
    }
    if(b.num == 0){
      SPE <- NA
    }else{
      SPE <- (TN/b.num) %>% round(., digits = 3)
    }
    ACC <- (sum(ifelse(dat.sub$SIGMA > cutoff, 
                       "Pathogenic", "Benign") == dat.sub$label) / 
              nrow(dat.sub)) %>% round(., digits = 3) 
    AUC <- NA
    if((b.num > 0) & (p.num > 0)){ 
      AUC <- pROC::roc(dat.sub$label, dat.sub$SIGMA, direction = "<",
                       quiet = T)$auc %>% as.numeric() %>% round(.,digits = 3)
    }
    gene.info[i, 2:10] <- c(vrt.num, b.num, p.num, TP, TN, SEN, SPE, ACC, AUC)
    

    dat.sub <- dat[dat$vep.SYMBOL %in% genes[i], 1:4] %>% na.omit 
    if(nrow(dat.sub)!=0){
      vrt.num <- nrow(dat.sub)
      b.num <- sum(dat.sub$label=="Benign")
      p.num <- sum(dat.sub$label=="Pathogenic")
      TP <- sum((dat.sub$vep.EVE > cutoff.eve) & (dat.sub$label=="Pathogenic"))
      TN <- sum((dat.sub$vep.EVE <= cutoff.eve) & (dat.sub$label=="Benign"))
      if(p.num == 0){
        SEN <- NA
      }else{
        SEN <- (TP/p.num) %>% round(., digits = 3)
      }
      if(b.num == 0){
        SPE <- NA
      }else{
        SPE <- (TN/b.num) %>% round(., digits = 3)
      }
      ACC <- (sum(ifelse(dat.sub$vep.EVE > cutoff.eve, 
                         "Pathogenic", "Benign") == dat.sub$label) /
                nrow(dat.sub)) %>% round(., digits = 3) 
      AUC <- NA
      if((b.num > 0) & (p.num > 0)){ 
        AUC <- pROC::roc(dat.sub$label, dat.sub$vep.EVE, direction = "<",
                         quiet = T)$auc %>% as.numeric() %>% round(.,digits = 3)
      }
      gene.info[i, 11:19] <- c(vrt.num, b.num, p.num, TP, TN, SEN, SPE, ACC, AUC)
      
    }
    
  
    dat.sub <- dat[dat$vep.SYMBOL %in% genes[i], c(1:3,5)] %>% na.omit 
    if(nrow(dat.sub)!=0){ 
      vrt.num <- nrow(dat.sub)
      b.num <- sum(dat.sub$label=="Benign")
      p.num <- sum(dat.sub$label=="Pathogenic")
      TP <- sum((dat.sub$vep.REVEL_score > cutoff.revel) & (dat.sub$label=="Pathogenic"))
      TN <- sum((dat.sub$vep.REVEL_score <= cutoff.revel) & (dat.sub$label=="Benign"))
      if(p.num == 0){
        SEN <- NA
      }else{
        SEN <- (TP/p.num) %>% round(., digits = 3)
      }
      if(b.num == 0){
        SPE <- NA
      }else{
        SPE <- (TN/b.num) %>% round(., digits = 3)
      }
      ACC <- (sum(ifelse(dat.sub$vep.REVEL_score > cutoff.revel, 
                         "Pathogenic", "Benign") == dat.sub$label) /
                nrow(dat.sub)) %>% round(., digits = 3) 
      AUC <- NA
      if((b.num > 0) & (p.num > 0)){
        AUC <- pROC::roc(dat.sub$label, dat.sub$vep.REVEL_score, direction = "<",
                         quiet = T)$auc %>% as.numeric() %>% round(.,digits = 3)
      }
      gene.info[i, 20:28] <- c(vrt.num, b.num, p.num, TP, TN, SEN, SPE, ACC, AUC)
    }
  }
  gene.info <- gene.info[order(gene.info$vrt.num, decreasing = T),]
  p.hist <- ggplot(data = gene.info[gene.info$vrt.num>=50,], aes(x = ACC)) +
    geom_histogram(bins = 20, color = "black", fill = "#86aab8", size = 0.1) +
    xlab("Accuracy") +
    ylab("Frequency") +
    theme_pubclean()
  
  
  tmp <- gene.info[gene.info$vrt.num>=50,]
  r <- cor(tmp$vrt.num, tmp$ACC, method = "spearman") %>% round(., digits = 3)
  r.p <- cor.test(tmp$vrt.num, tmp$ACC, method = "spearman")$p.value %>% signif(., digits = 3)
  dat_text <- data.frame(
    label = paste0("r = ", r, ", P = ", r.p))
  p.point <- ggplot(data = tmp, aes(x = vrt.num, y = ACC)) +
    geom_point(colour = "#335b68",size = 0.8)+
    labs(y="Accuracy",x="# mutations")+
    geom_smooth(method = lm,se=T, colour = "red", size = 0.5)+
    theme_bw()+
    theme(
      axis.text.x = element_text(color="black",size = 8),
      axis.text.y = element_text(color="black",size = 8),
      axis.title.y = element_text(color="black",size = 10),
      axis.title.x = element_text(color="black",size = 10),
      title = element_text(color = "black", size = 10),
      panel.border = element_blank(),panel.grid=element_blank(),
      axis.line = element_line(colour = "black", size = 0.2),
      axis.ticks = element_line(colour = "black", size = 0.2),
      strip.background = element_blank()
    )+
    geom_text(
      data    = dat_text,
      mapping = aes(x = 400, y = 1, label = label),
      size = 2.5,
      hjust   = -0.1
    )
  out <- list(gene.info = gene.info, p.hist = p.hist, p.point = p.point)
  return(out)
}


genelevel_acc_stat <- function(a, b, c){
  a <- a[a$b.num>=10,]
  a <- a[a$p.num>=10,]
  b <- b[b$b.num>=10,]
  b <- b[b$p.num>=10,]
  c <- c[c$b.num>=10,]
  c <- c[c$p.num>=10,]
  a$ACC %>% cut(., breaks = seq(0,1,0.1)) %>% 
    table() %>% rev() %>% cumsum() %>%
    divide_by(., 200) %>% multiply_by(., 100) %>%
    rev() %>% print
  b$ACC %>% cut(., breaks = seq(0,1,0.1)) %>% 
    table() %>% rev() %>% cumsum() %>%
    divide_by(., 200) %>% multiply_by(., 100) %>%
    rev() %>% print
  c$ACC %>% cut(., breaks = seq(0,1,0.1)) %>% 
    table() %>% rev() %>% cumsum() %>%
    divide_by(., 200) %>% multiply_by(., 100) %>%
    rev() %>% print
}

