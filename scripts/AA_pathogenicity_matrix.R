# amino change bubble plot
amino_change_dotplot <- function(dat){
  dat <- dat[, c("label", "SIGMA", "vep.HGVSp")]
  dat <- data.frame(
    dat,
    ref = dat$vep.HGVSp %>% strsplit(., split = "p[.]") %>% sapply(., function(x) x[2]) %>%
      substr(., start = 1, stop = 3),
    alt = dat$vep.HGVSp %>% strsplit(., split = "p[.]") %>% sapply(., function(x) x[2]) %>%
      substr(., start = nchar(.)-2, stop = nchar(.)),
    stringsAsFactors = F
  )
  dat <- data.frame(
    dat,
    change = paste0(dat$ref, ">",dat$alt),
    stringsAsFactors = F
  )
  change <- unique(dat$change) 
  dat.score <- matrix(NA, nrow = 20, ncol = 20) 
  colnames(dat.score) <- unique(dat$ref) %>% sort()
  rownames(dat.score) <- unique(dat$ref) %>% sort()
  dat.vrtNum <- matrix(0, nrow = 20, ncol = 20) 
  colnames(dat.vrtNum) <- unique(dat$ref) %>% sort()
  rownames(dat.vrtNum) <- unique(dat$ref) %>% sort()
  for(i in 1:length(change)){
    ref <- change[i] %>% substr(., start = 1, stop = 3)
    alt <- change[i] %>% substr(., start = 5, stop = 7)
    change.dat <- dat[dat$change == change[i],]
    dat.score[ref, alt] <- change.dat$SIGMA %>% median()
    dat.vrtNum[ref, alt] <- nrow(change.dat)
  }
  

  dat.plot.score <- dat.score %>% reshape2::melt()
  dat.plot.vrtNum <- dat.vrtNum %>% reshape2::melt()
  dat.plot <- cbind(dat.plot.vrtNum, dat.plot.score$value)
  colnames(dat.plot) <- c("ref", "alt", "Number", "SIGMA")
  

  dat.plot <- na.omit(dat.plot)
  amino.ord <- c("Arg", "Lys", "His", # pos
                 "Asp", "Glu", # neg
                 "Asn", "Gln", 
                 "Trp", "Tyr", "Phe", 
                 "Gly", "Pro", "Ser", "Thr",
                 "Val", "Cys", "Ile", "Leu", "Met", "Ala" 
  )
  amino.ord <- c("Arg", "Lys",
                 "Asp", "Glu", 
                 "Asn", "Gln", "His", 
                 "Pro", "Ser", "Thr", "Gly","Tyr",
                 "Trp", "Ala", "Met", "Cys", "Phe", "Leu", 
                 "Val", "Ile" 
  )
  amino.ord <- c(amino.ord, setdiff(unique(dat.plot$ref), amino.ord))
  dat.plot$ref <- factor(dat.plot$ref, levels = amino.ord)
  dat.plot$alt <- factor(dat.plot$alt, levels = amino.ord)
  p <- ggplot(data = dat.plot, aes(x = ref, y = alt, size = Number, fill = SIGMA)) +
    geom_point(shape = 21, color = "black") +
    scale_size(range = c(2,7)) +
    xlab("Reference AA") +
    ylab("Alternative AA") +
    scale_fill_viridis_c(option = "magma", direction = -1) +
    theme_linedraw()
  return(p)
}
