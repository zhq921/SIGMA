DMS.veps.com_plot <- function(dat){
  dat$Method[dat$Method == "Combined_score"] <- "SIGMA+"
  dat$Method <- gsub("_converted", "", dat$Method)
  dat$Method <- factor(dat$Method, levels = unique(dat$Method))
  vep_clas <- read.table("./data/other_data/VEPs_classification.tsv", 
                         sep = "\t", header = T, quote = "",
                         stringsAsFactors = F)
  vep_clas$VEPs <- gsub("vep.","",vep_clas$VEPs)
  vep_clas$VEPs <- gsub("_converted","",vep_clas$VEPs)
  vep_clas$VEPs <- gsub("_rankscore","",vep_clas$VEPs)
  vep_clas <- data.frame(vep_clas,
                         highlight = 
                           ifelse(grepl("SIGMA", vep_clas$VEPs), "yes", "no"),
                         stringsAsFactors = F)
  dat <- data.frame(dat, 
                    vep_clas[match(dat$Method, vep_clas$VEPs),c("Label", "Group", "highlight")])
  dat$Gene <- factor(dat$Gene, levels = unique(dat$Gene))
  dat.ind <- dat[dat$Group == "Individual predictor",]
  p1 <- ggplot(data = dat[dat$Group == "Individual predictor",], 
               mapping = aes(x = Gene, y = Correlation, colour = Method, shape = highlight))+
    geom_jitter(width = 0.2)+
    # facet_wrap(~ Group, scales = "free", nrow = 2) +
    xlab("Protein") +
    scale_color_manual(values = 
                         c(pal_npg( )(10)[8],
                           pal_jama()(6)[c(2:3, 5:6)],
                           pal_npg( )(10)[c(2:7,9:10)],
                           pal_nejm()(8)[c(2,4:8)],
                           pal_aaas()(8)
                         ),
                       name = "VEPs"
    ) +
    scale_shape_manual(values = c(yes = 18, no = 18)) + 
    guides(colour=guide_legend(nrow=5,byrow = F, keyheight = 0.11,
                               label.theme = element_text(size = 6),
                               title.position = "top",
                               title.theme = element_text(size = 7),
    ),
    shape = "none")+
    theme_bw()+
    theme(legend.spacing = unit(0.2, units ="pt"),
          legend.key.size = unit(0.15, "cm"),
          legend.position = "bottom",
          legend.margin = margin(t=-0.4, unit = "cm"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text.x = element_text(color="black",size = 7),
          axis.text.y = element_text(color="black",size = 8),
          axis.title.y = element_text(color="black",size = 10),
          axis.title.x = element_text(color="black",size = 10),
          panel.border = element_blank(),panel.grid=element_blank(),
          axis.line = element_line(colour = "black", size = 0.3),
          axis.ticks = element_line(colour = "black", size = 0.2)
    )
  p2 <- ggplot(data = dat[dat$Group == "meta-predictor",], 
               mapping = aes(x = Gene, y = Correlation, colour = Method, shape = highlight))+
    geom_jitter(width = 0.2)+
    xlab("Protein") +
    scale_color_manual(values = 
                         c(pal_npg( )(10)[8],
                           pal_jama()(6)[c(2:3, 5:6)],
                           pal_npg( )(10)[c(2:7,9:10)],
                           pal_nejm()(8)[c(2,4:8)],
                           pal_aaas()(8)
                         ),
                       name = "VEPs"
    ) +
    scale_shape_manual(values = c(yes = 18, no = 18)) + 
    guides(colour=guide_legend(nrow = 3,byrow = F, keyheight = 0.1,
                               label.theme = element_text(size = 6),
                               title.position = "top",
                               title.theme = element_text(size = 7),
    ),
    shape = "none")+
    theme_bw()+
    theme(legend.spacing = unit(0.2, units ="pt"),
          legend.key.size = unit(0.15, "cm"),
          legend.position = "bottom",
          legend.margin = margin(t=-0.5, unit = "cm"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text.x = element_text(color="black",size = 7),
          axis.text.y = element_text(color="black",size = 8),
          axis.title.y = element_text(color="black",size = 10),
          axis.title.x = element_text(color="black",size = 10),
          panel.border = element_blank(),panel.grid=element_blank(),
          axis.line = element_line(colour = "black", size = 0.3),
          axis.ticks = element_line(colour = "black", size = 0.2)
    )
  p <- list(p1 = p1, p2 = p2)
  return(p)
}

DMS.veps.rankflow <- function(dat){
  dat$Method[dat$Method == "Combined_score"] <- "SIGMA+"
  dat$Method <- gsub("_converted", "", dat$Method)
  dat$Method <- factor(dat$Method, levels = unique(dat$Method))
  vep_clas <- read.table("./data/other_data/VEPs_classification.tsv", 
                         sep = "\t", header = T, quote = "",
                         stringsAsFactors = F)
  vep_clas$VEPs <- gsub("vep.","",vep_clas$VEPs)
  vep_clas$VEPs <- gsub("_converted","",vep_clas$VEPs)
  vep_clas$VEPs <- gsub("_rankscore","",vep_clas$VEPs)
  vep_clas <- data.frame(vep_clas,
                         highlight = 
                           ifelse(grepl("SIGMA", vep_clas$VEPs), "yes", "no"),
                         stringsAsFactors = F)
  dat <- data.frame(dat, 
                    vep_clas[match(dat$Method, vep_clas$VEPs),c("Label", "Group", "highlight")])
  dat$Gene <- factor(dat$Gene, levels = rev(unique(dat$Gene)))
  dat$Correlation[is.na(dat$Correlation)] <- 0
  library(dplyr)
  dat.ind <- dat[dat$Group == "Individual predictor",]
  dat.ind.rank <- dat.ind %>%
    group_by(Gene) %>% 
    arrange(Gene, desc(Correlation), Method) %>% 
    mutate(ranking = row_number(),
           gene = Gene) %>% 
    as.data.frame()
  # Create rankflow plot using ggplot2
  my_theme <- function() {
    
    # Colors
    color.background = "transparent"
    color.text = "#000000"
    
    # Begin construction of chart
    theme_bw(base_size=15) +
      
      # Format background colors
      theme(panel.background = element_rect(fill=color.background, color=color.background)) +
      theme(plot.background  = element_rect(fill=color.background, color=color.background)) +
      theme(panel.border     = element_rect(color=color.background)) +
      theme(strip.background = element_rect(fill=color.background, color=color.background)) +
      
      
      # Format the grid
      theme(panel.grid.major.y = element_blank()) +
      theme(panel.grid.minor.y = element_blank()) +
      theme(axis.ticks       = element_blank()) +
      
      # Format the legend
      theme(legend.position = "none") +
      
      # Format title and axis labels
      theme(plot.title       = element_text(color=color.text, size=20)) +
      theme(axis.title.x     = element_text(size=10, color="black")) +
      theme(axis.title.y     = element_text(size=10, color="black", vjust=1.25)) +
      theme(axis.text.x      = element_text(size=6.4, vjust=0.5, hjust=0.5, color = color.text)) +
      theme(axis.text.y      = element_text(size=7, vjust=0.5, hjust=0.5, color = color.text)) +
      theme(strip.text       = element_text(face = "bold")) +
      
      # Plot margins
      theme(plot.margin = unit(c(0.2, 0.1, 0.15, 0.2), "cm"))
  }
  p1 <- ggplot(data = dat.ind.rank, aes(x = gene, 
                                        y = ranking, 
                                        group = Method)) +
    geom_line(aes(color = Method), size = 1, alpha = 0.7) +
    geom_point(aes(color = Method), size = 3-0.4, alpha = 0.75) +
    geom_point(color = "#FFFFFF", size = 0.8-0.4) +
    scale_y_reverse(
      breaks = 1:18,
      sec.axis = sec_axis(~ ., 
                          breaks = 1:18,
                          labels = dat.ind.rank$Method[dat.ind.rank$Gene == "All"]))+ 
    scale_x_discrete(expand = c(0, 0.16)) +
    theme(legend.position = "none") +
    labs(x = "Protein",
         y = "Rank") +
    scale_color_manual(values = 
                         c(pal_npg( )(10)[8],
                           pal_jama()(6)[c(2:3, 5:6)],
                           pal_npg( )(10)[c(2:7,9:10)],
                           pal_nejm()(8)[c(2,4:8)],
                           pal_aaas()(8)
                         ),
                       name = "VEPs"
    )+
    my_theme()
  
  dat.meta <- dat[dat$Group == "meta-predictor",]
  dat.meta.rank <- dat.meta %>%
    group_by(Gene) %>% 
    arrange(Gene, desc(Correlation), Method) %>% 
    mutate(ranking = row_number(),
           gene = Gene) %>% 
    as.data.frame()
  
  p2 <- ggplot(data = dat.meta.rank, aes(x = gene, 
                                         y = ranking, 
                                         group = Method)) +
    geom_line(aes(color = Method, alpha = 0), size = 2) +
    geom_point(aes(color = Method, alpha = 0), size = 4) +
    geom_point(color = "#FFFFFF", size = 1) +
    scale_y_reverse(
      breaks = 1:12,
      sec.axis = sec_axis(~ ., 
                          breaks = 1:12,
                          labels = dat.meta.rank$Method[dat.meta.rank$Gene == "All"]))+ 
    scale_x_discrete(expand = c(0, 0.21)) +
    theme(legend.position = "none") +
    labs(x = "Protein",
         y = "Rank") +
    scale_color_manual(values = 
                         c(pal_npg( )(10)[8],
                           pal_jama()(6)[c(2:3, 5:6)],
                           pal_npg( )(10)[c(2:7,9:10)],
                           pal_nejm()(8)[c(2,4:8)],
                           pal_aaas()(8)
                         ),
                       name = "VEPs"
    )+
    my_theme()
  
  
  p <- list(p1 = p1, p2 = p2)
  return(p)
}