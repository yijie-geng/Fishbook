

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler) 
library(readr)
library(GenomicRanges) 
library(BSgenome.Hsapiens.UCSC.hg38)
library(bumphunter)
library(org.Hs.eg.db)
genes <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)
library(stringr)
library(reshape2)


setwd("") # set working directory
getwd()

# read ChIPseq data .bed files, all based on hg38
TOP2A=readPeakFile( "" ) # GSE79593 TOP2A_DMSO sample
EZH2=readPeakFile( "" ) # GSM1003576 
H3K27me3=readPeakFile( "" ) # GSM733658 
SUZ12=readPeakFile( "" ) # GSM1003545 

overlap_gr<-EZH2
map <- nearest(overlap_gr, genes)
overlap_gr_ann <- matchGenes(overlap_gr[!is.na(map),], genes)
x <- data.frame(overlap_gr[!is.na(map),])
overlap_genes <- cbind(x[, c(1:3)], overlap_gr_ann[,-c(2,13,15)])
write.csv(overlap_genes, file="EZH2.csv", row.names=FALSE)

overlap_gr<-H3K27me3
map <- nearest(overlap_gr, genes)
overlap_gr_ann <- matchGenes(overlap_gr[!is.na(map),], genes)
x <- data.frame(overlap_gr[!is.na(map),])
overlap_genes <- cbind(x[, c(1:3)], overlap_gr_ann[,-c(2,13,15)])
write.csv(overlap_genes, file="H3K27me3.csv", row.names=FALSE)

overlap_gr<-SUZ12
map <- nearest(overlap_gr, genes)
overlap_gr_ann <- matchGenes(overlap_gr[!is.na(map),], genes)
x <- data.frame(overlap_gr[!is.na(map),])
overlap_genes <- cbind(x[, c(1:3)], overlap_gr_ann[,-c(2,13,15)])
write.csv(overlap_genes, file="SUZ12.csv", row.names=FALSE)


base=EZH2
overlap <-  findOverlaps(base, TOP2A)
overlap_gr <- base[queryHits(overlap),]
map <- nearest(overlap_gr, genes)
overlap_gr_ann <- matchGenes(overlap_gr[!is.na(map),], genes)
x <- data.frame(overlap_gr[!is.na(map),])
overlap_genes <- cbind(x[, c(1:3)], overlap_gr_ann[,-c(2,13,15)])
write.csv(overlap_genes, file="TOP2A_EZH2.csv", row.names=FALSE)

base=H3K27me3
overlap <-  findOverlaps(base, TOP2A)
overlap_gr <- base[queryHits(overlap),]
map <- nearest(overlap_gr, genes)
overlap_gr_ann <- matchGenes(overlap_gr[!is.na(map),], genes)
x <- data.frame(overlap_gr[!is.na(map),])
overlap_genes <- cbind(x[, c(1:3)], overlap_gr_ann[,-c(2,13,15)])
write.csv(overlap_genes, file="TOP2A_H3K27me3.csv", row.names=FALSE)

base=SUZ12
overlap <-  findOverlaps(base, TOP2A)
overlap_gr <- base[queryHits(overlap),]
map <- nearest(overlap_gr, genes)
overlap_gr_ann <- matchGenes(overlap_gr[!is.na(map),], genes)
x <- data.frame(overlap_gr[!is.na(map),])
overlap_genes <- cbind(x[, c(1:3)], overlap_gr_ann[,-c(2,13,15)])
write.csv(overlap_genes, file="TOP2A_SUZ12.csv", row.names=FALSE)




