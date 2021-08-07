

setwd("") # set working directory
getwd()

path <- "" # set path to a GMT format GLAD4U gene set library .csv file;
data <- read.csv(path, header=T, stringsAsFactors = FALSE)

diseases <- c("Autistic Disorder", "Autism Spectrum Disorder", "Schizophrenia", "Bipolar Disorder", "Mental Retardation", "Depression")

path <- "" # set path to "zebrafish genes human homologs.csv" file, a list containing human homologs of all zebrafish genes 
genes <- read.csv(path, header=F)$V1 


for (i in diseases) {
  
  row <- which( unlist(unname(data[,1])) == i )
  x <- unlist(unname(data[row,2:ncol(data)]))
  x <- x[!is.na(x)] 
  x <- x[which(x != "")]
  
  df <- data.frame( matrix(0, ncol = 2, nrow = 1) )
  colnames(df) <- c("permutation", "odds.ratio")
  
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
  write.csv(df, paste("Permutation ", i, " GLAD4U.csv", sep=""), row.names=FALSE)
  
}
