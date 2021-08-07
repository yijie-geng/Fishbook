


library(readr)
library(dplyr)
library(enrichR)


### use gene-sets-of-interest to query ChIP-X, ENCODE_HM, and Epigenomics_Roadmap_HM libraries in Enrichr 

setwd("") # set working directory
getwd()

dbs <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "ENCODE_Histone_Modifications_2015", "Epigenomics_Roadmap_HM_ChIP-seq")

PATH <- "" # set path to query gene sets in individual .csv files, such as can4Dn.csv and can4Up.csv; include a SFARI.csv file for subsequent ChIP-X UES-blast analysis
filenames <- list.files(PATH, pattern="*.csv")

for (i in filenames) {
  filepath <- file.path(PATH, i)
  x <- read.csv(filepath, header=F)$V1
  x <- as.character(x)
  
  enriched <- enrichr(x, dbs)
  
  for (j in dbs) {
    df <- enriched[[j]]
    ordered_df <- df[order(df$Adjusted.P.value), c(1,3,4)]
    write.csv(ordered_df, file = paste(as.character(j), " ", as.character(i), sep=""), row.names=FALSE)
  }
}


### copy and save results to separate folders named by the 3 libraries

fromDir <- getwd()
dataDir <- "" # set destination folder

for (j in dbs){
  
  toDir <- file.path(dataDir, j) 
  if (!(dir.exists(toDir))) dir.create(toDir)
  
  filenames <- list.files( fromDir, pattern=paste("^", j, ".*\\.csv$", sep="") )
  
  for (i in filenames) {
    from <- file.path(fromDir, i)
    file.copy(from, toDir)
  }
  
}

### ChIP-X UES-blast against SFARI.csv

dbs <- c("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

fileDir <- "" # set destination to Enrichr libraries folder: download ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X, 
              # ENCODE_Histone_Modifications_2015, and Epigenomics_Roadmap_HM_ChIP-seq library files from Enrichr 
              # website: https://maayanlab.cloud/Enrichr/#libraries, save in this folder

for (j in dbs){
  
  filepath <- file.path(fileDir, paste(j, ".txt", sep=""))
  fn <- read.delim(filepath, header=F, stringsAsFactors = FALSE, col.names = paste0("V",seq_len( max(count.fields(filepath)) )), fill = TRUE)
  fn <- as.character(fn[,1])
  fn <- c("Disease", fn)
  fn <- unique(fn)
  
  df <- data.frame(matrix(0, ncol = length(fn), nrow = 1))
  colnames(df) <- fn
  
  datapath <- file.path(dataDir, j)
  datanames <- list.files( datapath, pattern="^.*\\.csv$" )
  
  for (i in datanames) {
    
    datafilepath <- file.path(datapath, i)
    x <- read.csv(datafilepath)
    
    y <- as.list( rep(0, length(fn)) )
    y[[1]] <- i
    
    for (k in 2:length(fn)) {
      if ( any( x[[1]] == fn[[k]] ) == FALSE ) next
      y[[k]] <- x[ x[[1]] == fn[[k]], c(3) ]
    }
    
    df <- rbind(df,y)
    
  }
  
  df <- df[-1,]
  write.csv(df, paste(j, " padj.csv", sep=""), row.names=FALSE)
  
  
  # calculate UES, and blast
  
  sfariRow <- which( df$Disease == paste(j, " SFARI.csv", sep="") )

  data <- select(df, c(2:length(fn)))
  data <- replace(data, data >= 0.05, 0)
  data <- -log(data)
  data <- replace(data, data== Inf, 0)
  
  for(m in 1:nrow(data)){
    if (max(data[m,]) == 0) next
    data[m,] <- data[m,] / max(data[m,])
  }
  
  sfari_dist_list <- list()
  
  for (n in 1:nrow(data)) {
    sfari_dist_list[[n]] <- dist(rbind(data[sfariRow,], data[n,]))[1]
  }
  
  data$disease <- df[,1]
  data$distance.sfari <- unlist(sfari_dist_list)
  
  ordered_data <- data[ order(data$distance.sfari), ]
  write.csv(ordered_data, file = paste(j, " UES-blast.csv", sep=""), row.names=FALSE) # UES and UES-blast result
  
}


### H3K27me3 score calculation

dbs <- c("ENCODE_Histone_Modifications_2015", "Epigenomics_Roadmap_HM_ChIP-seq")

for (j in dbs){
  
  filepath <- file.path(fileDir, paste(j, ".txt", sep=""))
  fn <- read.delim(filepath, header=F, stringsAsFactors = FALSE, col.names = paste0("V",seq_len( max(count.fields(filepath)) )), fill = TRUE)
  fn <- as.character(fn[,1])
  fn <- c("Disease", fn)
  fn <- unique(fn)
  
  df <- data.frame(matrix(0, ncol = length(fn), nrow = 1))
  colnames(df) <- fn
  
  datapath <- file.path(dataDir, j)
  datanames <- list.files( datapath, pattern="^.*\\.csv$" )
  
  for (i in datanames) {
    
    datafilepath <- file.path(datapath, i)
    x <- read.csv(datafilepath)
    
    y <- as.list( rep(1, length(fn)) )
    y[[1]] <- i
    
    for (k in 2:length(fn)) {
      if ( any( x[[1]] == fn[[k]] ) == FALSE ) next
      y[[k]] <- x[ x[[1]] == fn[[k]], c(3) ]
    }
    
    df <- rbind(df,y)
    
  }
  
  df <- df[-1,]
  write.csv(df, paste(j, " padj.csv", sep=""), row.names=FALSE)
  
}


### ENCODE_HM

path <- "./ENCODE_Histone_Modifications_2015 padj.csv"

data <- read.csv(path, header = T, stringsAsFactors = FALSE)

cn <- colnames(data)
cp <- grep("27me3", cn)
c <- length(cp) # total number of H3K27me3 datasets

df <- data

score_list <- list()
a_list <- list()
b_list <- list()

for (i in 1:nrow(data)){
  a <- length(which(data[i,cp] < 0.05)) #total number of significant H3K27me3 hits
  b <- length(which(data[i,2:413] < 0.05)) #total number of significant hits out of 412 total in ENCODE_HM library
  a_list[[i]] <- a
  b_list[[i]] <- b
  if (b == 0){
    #score_list[[i]] <- "NA"
    score_list[[i]] <- 0
    next
  }
  score_list[[i]] <- (a/c) * (a/b)
}

df$H3k27me3.rank <- unlist(score_list)
df$H3k27me3.hits <- unlist(a_list)
df$all.hits <- unlist(b_list)

ordered_data <- df[ order(df$H3k27me3.rank, decreasing = T), ]

write.csv(ordered_data, "ENCODE_Histone_Modifications_2015 H3K27me3 score.csv", row.names=FALSE) # ENCODE_HM H3K27me3 score result


### Epigenomics_Roadmap_HM

path <- "./Epigenomics_Roadmap_HM_ChIP-seq padj.csv" 

data <- read.csv(path, header = T, stringsAsFactors = FALSE)

cn <- colnames(data)
cp <- grep("27me3", cn)
c <- length(cp) # total number of H3K27me3 datasets

df <- data

score_list <- list()
a_list <- list()
b_list <- list()

for (i in 1:nrow(data)){
  a <- length(which(data[i,cp] < 0.05)) #total number of significant H3K27me3 hits
  b <- length(which(data[i,2:384] < 0.05)) #total number of significant hits out of 383 total in Epigenomics_Roadmap_HM library
  a_list[[i]] <- a
  b_list[[i]] <- b
  if (b == 0){
    #score_list[[i]] <- "NA"
    score_list[[i]] <- 0
    next
  }
  score_list[[i]] <- (a/c) * (a/b)
}

df$H3k27me3.rank <- unlist(score_list)
df$H3k27me3.hits <- unlist(a_list)
df$all.hits <- unlist(b_list)

ordered_data <- df[ order(df$H3k27me3.rank, decreasing = T), ]

write.csv(ordered_data, "Epigenomics_Roadmap_HM_ChIP-seq H3K27me3 score.csv", row.names=FALSE) # Epigenomics_Roadmap_HM H3K27me3 score result

