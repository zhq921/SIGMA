# AUC plot
myroc <- function(res, pre, tit, text.size, title.margin.tb, title.x.margin.b){
  roc.tmp <- pROC::roc(response = res, predictor = pre, direction = "<", auc=T,ci =T,plot = F,
                       col = "#29476a", lwd = 2, legacy.axes = T, quiet = T, main = tit, 
                       mar =c(5,3.5,4,1.7), mgp = c(2,0.5,0), ci.type = "bars", cex.axis = 1) # axis label 8 pt
  auc = paste0("AUC = ",round(roc.tmp$auc,3))
  ci = paste0("95% CI = ",round(roc.tmp$ci[1],3),"-",round(roc.tmp$ci[3],3))
  p = "P-value < 2.2e-16"
  out <- ggroc(roc.tmp, legacy.axes = T, color = "#3073af") + 
    geom_abline(intercept = 0, slope = 1, color = "#a9a9a9", size = 0.25) +
    annotate("text", x = 0.3, y = 0.2, label = auc, hjust = 0, size = text.size)+
    annotate("text", x = 0.3, y = 0.12, label = ci, hjust = 0, size = text.size)+
    annotate("text", x = 0.3, y = 0.04, label = p,  hjust = 0, size = text.size)+
    theme_bw()+
    ggtitle(tit) + 
    ylab("Sensitivity") +
    theme(
      axis.text = element_text(colour = "black", size = 8),
      axis.title.x = element_text(colour = "black", size = 10, margin = margin(0,0,title.x.margin.b,0)),
      axis.title.y = element_text(colour = "black", size = 10),
      axis.ticks = element_line(colour = "black", size = 0.2),
      plot.title = element_text(colour = "black", size = 10, hjust = 0.5, margin = margin(title.margin.tb,0,title.margin.tb,0)),
      panel.grid = element_line(size = 0.15, linetype = "dashed")
    )
  return(out)
}