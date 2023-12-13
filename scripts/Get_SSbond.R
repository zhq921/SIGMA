library(magrittr)
library(stringr)
library(getopt)

spec <- matrix(
  c("dis_cutoff",  "c", 2, "numeric",  "disulfide distance cutoff", 
    "pdbfile", "i", 1, "character", "the path for the input pdb file",
    "outfile", "o", 1, "character", "the path for the output file",
    "help", "h", 0, "logical",  "This is Help!"),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)

if( !is.null(opt$help) || is.null(opt$pdbfile) ||
    is.null(opt$outfile) ){
  cat(paste(getopt(spec=spec, usage = T), "\n"))
  quit()
}
if(is.null(opt$dis_cutoff)){
  dis_cutoff <- 3
}else{
  dis_cutoff <- opt$dis_cutoff
}

parse_pdb <- function(pdb){
  atom.num <- lapply(pdb, function(x) substr(x, start = 7, stop = 11)) %>% unlist() %>% as.numeric()
  atom <- lapply(pdb, function(x) substr(x, start = 13, stop = 16)) %>% unlist() %>% trimws(., which = "both")
  residue <- lapply(pdb, function(x) substr(x, start = 18, stop = 20)) %>% unlist() %>% trimws(., which = "both")
  residue.num <- lapply(pdb, function(x) substr(x, start = 23, stop = 26)) %>% unlist() %>% as.numeric()
  coord.x <- lapply(pdb, function(x) substr(x, start = 31, stop = 38)) %>% unlist() %>% as.numeric()
  coord.y <- lapply(pdb, function(x) substr(x, start = 39, stop = 46)) %>% unlist() %>% as.numeric()
  coord.z <- lapply(pdb, function(x) substr(x, start = 47, stop = 54)) %>% unlist() %>% as.numeric()
  plddt <- lapply(pdb, function(x) substr(x, start = 61, stop = 66)) %>% unlist() %>% as.numeric()
  out <- cbind(atom.num,atom,residue,residue.num,coord.x,coord.y,coord.z, plddt)
  return(out)
}

pdb <- readLines(opt$pdbfile, warn = F)
pdb <- grepl("^ATOM", x = pdb) %>% extract(pdb, .) 
pdb <- (substr(pdb, 18, 20)=="CYS") %>% extract(pdb, .) 
pdb <- parse_pdb(pdb) %>% as.data.frame(., stringsAsFactors = F)
pdb <- split(x = pdb, 
             f = factor(pdb[,4], 
                        levels = unique(pdb[,4]) %>% as.numeric() %>% sort() %>% as.character()
             )
) 

get_ss_dis <- function(s1, s2){
  (s1 - s2)^2 %>% sum() %>% return()
}

get_dihedral_angle <- function(c1, c2, s1, s2){
  x1 <- c1[1]; y1 <- c1[2]; z1 <- c1[3];
  x2 <- s1[1]; y2 <- s1[2]; z2 <- s1[3];
  x3 <- s2[1]; y3 <- s2[2]; z3 <- s2[3];
  A1 <- (y3 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1);
  B1 <- (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
  C1 <- (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
  
  x1 <- c2[1]; y1 <- c2[2]; z1 <- c2[3];
  x2 <- s1[1]; y2 <- s1[2]; z2 <- s1[3];
  x3 <- s2[1]; y3 <- s2[2]; z3 <- s2[3];
  A2 <- (y3 - y1)*(z3 - z1) - (z2 - z1)*(y3 - y1);
  B2 <- (x3 - x1)*(z2 - z1) - (x2 - x1)*(z3 - z1);
  C2 <- (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
  cos_theta <- (A1*A2+B1*B2+C1*C2) / (sqrt(A1^2+B1^2+C1^2) * sqrt(A2^2+B2^2+C2^2))
  theta <- (acos(cos_theta)/pi) * 180
  return(theta)
}

dis_cut_2 <- dis_cutoff^2
free_cys <- names(pdb)
disulfide <- c()
cys.num <- length(pdb)
if(cys.num %in% c(0,1)){
  write.table("No disulfide bond found!",
              file = opt$outfile, row.names = F, col.names = F, quote = F)
  q(save = "no")
}
for(i in 1:(cys.num-1)){
  if(is.na(free_cys[i])){next}
  for(j in (i+1):cys.num){
    if(is.na(free_cys[j])){next}

    s1 <- pdb[[i]][pdb[[i]]$atom=="SG", 5:7] %>% as.numeric()
    s2 <- pdb[[j]][pdb[[j]]$atom=="SG", 5:7] %>% as.numeric()
    ss_dis <- get_ss_dis(s1 = s1, s2 = s2)
    if(ss_dis >= dis_cut_2){next} 
    c1 <- pdb[[i]][pdb[[i]]$atom=="CB", 5:7] %>% as.numeric() 
    c2 <- pdb[[j]][pdb[[j]]$atom=="CB", 5:7] %>% as.numeric() 
    theta <- get_dihedral_angle(c1 = c1, c2 = c2, s1 = s1, s2 = s2)
    if( (theta < 60) | (theta > 120)){next} 
    disulfide <- rbind(disulfide, c(
      names(pdb)[i], names(pdb)[j], pdb[[i]][1, "plddt"], pdb[[j]][1, "plddt"]
    ))
    free_cys[i] <- NA; free_cys[j] <- NA
  }
}
if(!is.null(disulfide)){
  colnames(disulfide) <- c("Res1_pos", "Res2_pos", "Res1_pLDDT", "Res2_pLDDT")
  write.table(disulfide, file = opt$outfile, sep = "\t", row.names = F, col.names = T, quote = F)
}else{
  write.table("No disulfide bond found!",
              file = opt$outfile, row.names = F, col.names = F, quote = F)
}
