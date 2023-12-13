combination_pretest <- function(dat){
  dat <- data.frame(
    alleleID = dat$alleleID,
    DEOGEN2 = dat$vep.DEOGEN2_rankscore,
    SIGMA = dat$SIGMA,
    PROVEAN = dat$vep.PROVEAN_converted_rankscore,
    MutPred = dat$vep.MutPred_rankscore,
    label = dat$label,
    group = dat$train_test,
    EVE = dat$vep.EVE,
    stringsAsFactors = F
  )
  
  dat <- na.omit(dat) 
  table(dat$group, dat$label)
  
  predictors <- list(c("DEOGEN2"),
                     c("SIGMA"),
                     c("EVE"),
                     c("PROVEAN"),
                     c("MutPred"),
                     c("SIGMA", "DEOGEN2"),
                     c("SIGMA", "EVE"),
                     c("SIGMA", "PROVEAN"),
                     c("SIGMA", "MutPred"),
                     c("SIGMA", "DEOGEN2", "EVE", "PROVEAN", "MutPred")
  )
  stat.dat <- data.frame()
  for(i in 1:length(predictors)){
    predictor <- predictors[[i]]
    if(length(predictor) == 1){ 
      train.pred <- dat[dat$group == "Training", predictor]
      train.label <- dat[dat$group == "Training", "label"] %>% as.factor()
    }else{
      # model
      h2o.init()
      dat.train <- as.h2o(dat[dat$group == "Training",])
      dat.test1 <- as.h2o(dat[dat$group == "Test",])
      # Change the column type to a factor:
      dat.train["label"] <- as.factor(dat.train["label"])
      dat.test1["label"] <- as.factor(dat.test1["label"])
      # logistic
      response <- c("label")
      aml <- h2o.automl(y = response, x = predictor,
                        training_frame = dat.train,
                        max_models = 15, 
                        include_algos = "GLM",
                        keep_cross_validation_predictions = T,
                        balance_classes = T,
                        nfolds = 5,
                        seed = 1)
      print(aml@leaderboard, n = 30)
      md <- aml@leader
      train.label <- dat.train["label"] %>% as.data.frame() %>% extract(.,, 1)
      train.pred =  h2o.cross_validation_holdout_predictions(md) %>% 
        as.data.frame() %>% extract(.,,"Pathogenic")
    }
    ## get auc
    roc.tmp <- pROC::roc(response = train.label, predictor = train.pred, 
                         direction = "<", auc=T,ci =T,plot = F)
    auc = round(roc.tmp$auc,3)
    auc.ci = c(round(roc.tmp$ci[1],3),round(roc.tmp$ci[3],3))
    auc <- c(auc,auc.ci)
    ## get cutoff
    cutoff <- cutpointr(x = train.pred, 
                        class = train.label, 
                        metric = youden)$optimal_cutpoint
    ## get acc
    cm.dat <- confusionMatrix(data = ifelse(train.pred >= cutoff,
                                            "Pathogenic", "Benign") %>% 
                                factor(., levels = c("Benign", "Pathogenic")),
                              reference = train.label,
                              positive = "Pathogenic")
    acc <- cm.dat$overall[1] %>% round(., digits = 3)
    acc.ci <- c(
      cm.dat$overall[3] %>% round(., digits = 3),
      cm.dat$overall[4] %>% round(., digits = 3)
    )
    acc <- c(acc, acc.ci)
    ## get sen
    sen.tmp <- c(sum( (train.pred>=cutoff) & (train.label=="Pathogenic") ), # TP
                 sum( (train.pred>=cutoff) & (train.label=="Benign") ), # FP
                 sum( (train.pred<cutoff) & (train.label=="Pathogenic") ), # FN
                 sum( (train.pred<cutoff) & (train.label=="Benign") ) # TN
    ) %>% matrix(., nrow = 2, byrow = T) %>% epi.tests(.)
    sen <- sen.tmp$detail$se %>% as.numeric() %>% round(., digits = 3)
    spe <- sen.tmp$detail$sp %>% as.numeric() %>% round(., digits = 3)
    ppv <- sen.tmp$detail$pv.pos %>% as.numeric() %>% round(., digits = 3)
    npv <- sen.tmp$detail$pv.neg %>% as.numeric() %>% round(., digits = 3)
    tmp <- matrix(c(auc,acc,sen,spe,ppv,npv),
                  byrow = T, ncol = 3) %>% as.data.frame()
    colnames(tmp) <- c("estimate", "lower", "upper")
    tmp <- data.frame(
      tmp,
      measure = c("AUC", "ACC", "SEN", "SPE", "PPV", "NPV"),
      model = paste(predictor, collapse = "+")
    )
    stat.dat <- rbind(stat.dat, tmp)
  }
  
  stat.dat <- stat.dat[stat.dat$measure %in% c("AUC", "SEN", "SPE"),]
  stat.dat$measure <- factor(stat.dat$measure, levels = 
                               c("AUC", "SEN", "SPE"))
  stat.dat$model <- as.character(stat.dat$model)
  stat.dat$model[stat.dat$model == "SIGMA+DEOGEN2+EVE+PROVEAN+MutPred"] <- 
    "SIGMA+"
  stat.dat$model <- factor(stat.dat$model, levels = 
                             c("PROVEAN", "EVE", "DEOGEN2", 
                               "MutPred", "SIGMA", 
                               "SIGMA+DEOGEN2", 
                               "SIGMA+MutPred", "SIGMA+PROVEAN", 
                               "SIGMA+EVE", "SIGMA+"))
  p <- ggplot(data = stat.dat, aes(x = model, y = estimate, fill = measure)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.3,size = 0.1,
                  position=position_dodge(.9)) +
    geom_text(aes(label = estimate), hjust = 1.2,
              colour = "white", angle = 90,
              position = position_dodge(.9), size = 2)+
    coord_cartesian(ylim=c(0.5,1))+
    xlab("Evidence Combinations")+ylab("Estimate")+
    scale_fill_manual(values = pal_npg()(6)[c(1,3,4)])+
    theme_bw()+ 
    theme(axis.text.x=element_text(size=7, colour= "black",angle = 30, vjust = 1,hjust = 1),
          axis.text.y=element_text(size=8,face="plain",colour= "black"), 
          axis.title.y=element_text(size = 10,face="plain"), 
          axis.title.x=element_blank(), 
          axis.line = element_line(colour = "black", size=0.3), 
          axis.ticks = element_line(colour = "black", size = 0.2),
          legend.text=element_text(colour="black",  
                                   size=8),
          legend.title=element_blank(),
          legend.key.size = unit(0.5,"cm"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank())
  return(list(p = p, stat.dat = stat.dat))
}