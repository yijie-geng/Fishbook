

setwd("") # set working directory
getwd()

diseaseGeneSetDir <- "" # set path to independent disease gene sets, saved as individual .csv files
queryGeneSetDir <- "" # set path to query gene sets, saved as individual .csv files, such as can4Dn.csv and can4Up.csv

diseaseGeneSetNames <- list.files( diseaseGeneSetDir, pattern="^.*\\.csv$" )
queryGeneSetNames <- list.files( queryGeneSetDir, pattern="^.*\\.csv$" )

df <- data.frame( matrix(0, ncol = 4, nrow = 1) )
colnames(df) <- c("disease.gene.set", "query.gene.set", "Fisher.pval", "odds.ratio")

for (i in queryGeneSetNames) {
  
  filepath <- file.path(queryGeneSetDir, i)
  y <- read.csv(filepath, header=F)$V1
  y <- as.character(y)
  y <- unique(y)
  
  print(i)
  
  for (j in diseaseGeneSetNames) {
    
    filepath <- file.path(diseaseGeneSetDir, j)
    x <- read.csv(filepath, header=F)$V1
    x <- as.character(x)
    x <- unique(x)
    
    # fisher test
    a <- length( intersect(x, y) )
    b <- length( x ) - a
    c <- length( y ) - a
    d <- 20438 - a - b - c # 20438: total number of human coding genes, GRCh38.p13
    if (d < 0) {
      ft <- "NA"
    } 
    else {
      ft <- fisher.test(matrix(c(a,b,c,d),nrow=2,ncol=2),alternative="greater")$p.value
    }
    
    # odds ratio
    or <- ( a / length( x ) ) / ( (length( y ) - a) / ( 20438 - length( x ) ) )
    
    df <- rbind( df, list(i, j, ft, or) )
    
  }
  
}

df <- df[-1,]
write.csv(df, "ORA.csv", row.names=FALSE)

