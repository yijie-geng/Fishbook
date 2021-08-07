


setwd("") # set working directory
getwd()


myGeneSetDir <- "" # set path to query gene sets saved as individual .csv files, such as can4Dn.csv and can4Up.csv
myGeneSetNames <- list.files( myGeneSetDir, pattern="^.*\\.csv$" )

path <- "" # set path to a GMT format disease gene set library .csv file, such as the GLAD4U library
data <- read.csv(path, header=T, stringsAsFactors = FALSE)

df <- data.frame( matrix(0, ncol = 4, nrow = 1) )
colnames(df) <- c("disease.gene.set", "query.gene.set", "Fisher.pval", "odds.ratio")

for (j in myGeneSetNames) {
  
  filepath <- file.path(myGeneSetDir, j)
  y <- read.csv(filepath, header=F)$V1
  y <- as.character(y)
  
  for (i in 1:nrow(data)){
    x <- unlist(unname(data[i,2:ncol(data)]))
    x <- x[!is.na(x)] 
    x <- x[which(x != "")]
    disease <- unlist(unname(data[i,1]))
    
    a <- length( intersect(x, y) )
    b <- length( unique(x) ) - a
    c <- length( unique(y) ) - a
    d <- 20438 - a - b - c  # 20438: total number of human coding genes, GRCh38.p13
    ft <- fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2),alternative="greater")$p.value
    
    or <- ( a / length( unique(x) ) ) / ( (length( unique(y) ) - a) / ( 20438 - length( unique(x) ) ) ) 
    
    df <- rbind( df, list(disease, j, ft, or) )
    
  }
  
}

df <- df[-1,]
write.csv(df, "ORA.csv", row.names=FALSE)

