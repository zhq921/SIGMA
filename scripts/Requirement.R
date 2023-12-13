library(magrittr) # for pipe
library(reshape2)
library(corrplot)
library(h2o) # for machine learning
library(caret)
library(forestplot) # for forest plot
library(pROC) # for ROC curve
library(PRROC)
library(ggpubr)
library(verification)
library(pheatmap)
library(ggsci) # for color
library(viridis) # for color
library(data.table)
library(epiR)
library(cutpointr) # for cutoff selection
library(ggDCA) # for DCA curve
library(rstatix)

load("./data/clinvar_data_all.RData")
load("./data/DMS.RData") 
load("./data/dat_model.RData")
load("./data/model_performance.RData") 
load("./data/pathogenicityMatrix.RData")

aminoAcid_table <- matrix(
  c(
    "Alanine",       "丙氨酸",             "Ala", "A",
    "Arginine",      "精氨酸",             "Arg", "R",
    "Asparagine",    "天冬酰胺",           "Asn", "N",
    "Asparticacid",  "天冬氨酸",           "Asp", "D",
    "Cysteine",      "半胱氨酸",           "Cys", "C",
    "Glutamine",     "谷氨酰胺",           "Gln", "Q",
    "Glutamicacid",  "谷氨酸",             "Glu", "E",
    "Glycine",       "甘氨酸",             "Gly", "G",
    "Histidine",     "组氨酸",             "His", "H",
    "Isoleucine",    "异亮氨酸",           "Ile", "I",
    "Leucine",       "亮氨酸",             "Leu", "L",
    "Lysine",        "赖氨酸",             "Lys", "K",
    "Methionine",    "甲硫氨酸（蛋氨酸）", "Met", "M",
    "Phenylalanine", "苯丙氨酸",           "Phe", "F",
    "Proline",       "脯氨酸",             "Pro", "P",
    "Serine",        "丝氨酸",             "Ser", "S",
    "Threonine",     "苏氨酸",             "Thr", "T",
    "Tryptophan",    "色氨酸",             "Trp", "W",
    "Tyrosine",      "酪氨酸",             "Tyr", "Y",
    "Valine",        "缬氨酸",             "Val", "V"
  ),
  ncol = 4, byrow = T
) %>% as.data.frame(., stringsAsFactors = F)
colnames(aminoAcid_table) <- c("english_name","chinese_name","abbr3","abbr1")