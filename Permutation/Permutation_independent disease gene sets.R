

setwd("") # set working directory
getwd()

diseaseGeneSetDir <- "" # set path to independent disease gene sets, saved as individual .csv files;
diseaseGeneSetNames <- list.files( diseaseGeneSetDir, pattern="^.*\\.csv$" )

path <- "" # set path to "zebrafish genes human homologs.csv" file, a list containing human homologs of all zebrafish genes 
genes <- read.csv(path, header=F)$V1 

for (i in diseaseGeneSetNames) {
  
  df <- data.frame( matrix(0, ncol = 2, nrow = 1) )
  colnames(df) <- c("permutation", "odds.ratio")
  
  filepath <- file.path(diseaseGeneSetDir, i)
  x <- read.csv(filepath, header=F)$V1
  x <- as.character(x)
  
  for (j in 1:1000) {
    y <- as.character(sample(genes, 5000))

    # odds ratio
    a <- length( intersect(x, y) )
    b <- length( unique(x) )
    c <- length( unique(y) )
    or <- ( a / b ) / ( ( c - a) / ( 20438 - b ) ) # 20438: total number of human coding genes, GRCh38.p13
    
    df <- rbind( df, list(j, or) )
    
  }
  
  df <- df[-1,]
  write.csv(df, paste("Permutation ", i, sep=""), row.names=FALSE)
  
}
