
min_max_scaling <- function(dat = dat){
  for(i in 1:ncol(dat)){
    x <- dat[,i]
    rg <- max(x) - min(x)
    if(rg==0){
      x <- x
    }else{
      x <- (x - min(x)) / (max(x) - min(x))
    }
    dat[,i] <- x
  }
  return(dat)
}


get_feature_importance <- function(dat, out.f){
  dat.train <- dat
  y <- "label"
  x <- setdiff(names(dat.train), y)
  dat.train <- as.h2o(dat.train);class(dat.train);h2o.describe(dat.train)
  h2o.levels(dat.train$label)
  nrow(dat.train)
  res.auc <- c()
  for(i in c("DRF")){
    if(i=="DeepLearning"){
      maxmd <- 1
    }else{
      maxmd <- 1
    }
    aml <- h2o.automl(y = y, x = x,
                      training_frame = dat.train,
                      max_models = maxmd, 
                      include_algos = i,
                      keep_cross_validation_predictions = T,
                      seed = 1, balance_classes = T)
    print(aml@leaderboard, n = 30)
    aml@leader %>% h2o.auc(., xval =T)
  }
  impt <- aml@leader@model$variable_importances %>% as.data.frame()
  h2o.performance(model = aml@leader, newdata = dat.train) %>% h2o.auc()
  pdf(file = out.f, height = 13, width = 6)
  va_plot <- h2o.varimp_plot(aml@leader, num_of_features = ncol(dat.train))
  dev.off()
  out <- list(rf_model = aml@leader, importance = impt)
  return(out)
}

# get OR for each feature
get_OR <- function(x, label, ft_name){
  if(all(unique(x) %in% c(0, 1))){
    or.dat <- table(x, label)
    or <- OR_calculator(or.dat)
  }else{
    logi <- glm(label ~ x, family = "binomial") %>% summary()
    or <- exp(logi$coefficients[2,1]) # OR logistic
    or.se <- logi$coefficients[2,2] # OR std err
    or.low <- exp(log(or) - 1.96*or.se)
    or.up <- exp(log(or) + 1.96*or.se)
    p <- logi$coefficients[2,4]
    
    if(p == 0){
      
      # Calculate the standard error of the log odds ratio
      SE_logOR <- (log(or.up) - log(or.low)) / (2 * 1.96)
      
      # Calculate the test statistic (Z)
      Z <- (log(or) - 0) / SE_logOR
      
      # Calculate the two-tailed p-value associated with the test statistic
      z <- mpfr(Z, precBits = 100)
      p <- 2 * pnorm(-abs(z))
      p <- format(p, digits = 3)
      
      # p <- logi$coefficients[2, 4]
      res <- c(signif(or, digits = 3), 
               signif(or.low,digits = 3), 
               signif(or.up,digits = 3),
               p)
    }else{
      p <- logi$coefficients[2, 4]
      res <- c(signif(or, digits = 3), 
               signif(or.low,digits = 3), 
               signif(or.up,digits = 3),
               signif(p,digits = 3))
    }
    or <- res
  }
  or <- matrix(or, nrow = 1)
  colnames(or) <- c("OR", "lower", "upper", "pValue")
  rownames(or) <- ft_name
  return(or)
}
# OR_calculator
OR_calculator <- function(dat){
  dat[1,1] <- as.numeric(dat[1,1]) 
  or <- (dat[2,2]*dat[1,1])/(dat[1,2]*dat[2,1])
  or.se <- sqrt(1/dat[1] + 1/dat[2] + 1/dat[3] +1/dat[4])
  or.low <- exp(log(or) - 1.96*or.se)
  or.up <- exp(log(or) + 1.96*or.se)
  chisq.stat <- (sum(dat) - 1) * (dat[4]*dat[1] - dat[2]*dat[3])^2 /
    ((dat[1]+dat[3]) * (dat[2]+dat[4]) * (dat[1]+dat[2]) * (dat[3]+dat[4]))
  # chisq.stat <- mpfr(chisq.stat, precBits = 100)
  p <- pchisq(chisq.stat, df = 1, lower.tail = F)
  if(p == 0){
    # Calculate the standard error of the log odds ratio
    SE_logOR <- (log(or.up) - log(or.low)) / (2 * 1.96)
    
    # Calculate the test statistic (Z)
    Z <- (log(or) - 0) / SE_logOR
    
    # Calculate the two-tailed p-value associated with the test statistic
    z <- mpfr(Z, precBits = 100)
    p <- 2 * pnorm(-abs(z))
    p <- format(p, digits = 3)
    
    res <- c(signif(or, digits = 3), 
             signif(or.low,digits = 3), 
             signif(or.up,digits = 3),
             p)
  }else{
    res <- c(signif(or, digits = 3), 
             signif(or.low,digits = 3), 
             signif(or.up,digits = 3),
             signif(p,digits = 3)
    )
  }
  
  return(res)
}

## forest plot
forest_plot <- function(or.ft){
  res <- as.data.frame(or.ft[,2:5], stringsAsFactors = F)
  colnames(res) <- c("HR", "HR_lower", "HR_upper", "P_value")
  res <- as.data.frame(res)[!res[,3]%in%Inf,]
  res <- res[order(as.numeric(res$HR)),]
  forest_table <- cbind(c("Feature", rownames(res)),
                        c("OR (95% CI)",paste0(res$HR," (",res$HR_lower,"-",res$HR_upper,")")),
                        c("P-value", res$P_value))
  
  csize <- data.frame(mean=c(NA, log10(as.numeric(res$HR))),
                      lower=c(NA, log10(as.numeric(res$HR_lower))),
                      upper=c(NA, log10(as.numeric(res$HR_upper))))
  p <- forestplot(labeltext = forest_table, 
                  csize,
                  graph.pos = 3,
                  graphwidth = unit(5, "cm"),
                  zero = 0,
                  lwd.zero = 1,
                  cex = 1,
                  lineheight = "auto",
                  boxsize = 0.1,
                  fn.ci_norm = fpDrawNormalCI,
                  lwd.ci = 1,
                  ci.vertices = F,
                  lwd.xaxis = 1,
                  xlab = "log10(OR)", 
                  ci.vertices.height = 0.2,
                  col = fpColors(box = "#e0171b", 
                                 line = "#e0171b",
                                 zero = "black"),
                  shapes_gp = fpShapesGp(zero = gpar(lty = "dashed")),
                  txt_gp = fpTxtGp(ticks = gpar(cex = 1),
                                   xlab = gpar(cex = 1))
  )
  return(p)
}


var_imp_plot <- function(Var_imp){
  dat.plot <- Var_imp[1:10,]
  dat.plot$variable <- gsub("Struc_", "", dat.plot$variable)
  dat.plot$variable <- gsub("FoldX_", "", dat.plot$variable)
  dat.plot$variable <- gsub("Total_delta", "Delta delta G", dat.plot$variable)
  dat.plot$variable <- gsub("_", " ", dat.plot$variable)
  dat.plot$variable <- factor(dat.plot$variable, levels = rev(dat.plot$variable))
  p<-ggplot(data = dat.plot, aes(x = variable, y = scaled_importance))+
    geom_bar(stat = "identity",width = 0.67, fill = "#447fb1", color = "black", size = 0.2) +
    ylab("Scale Importance") + xlab("Top 10 Features")+
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

