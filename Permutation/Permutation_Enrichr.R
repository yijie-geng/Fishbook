

setwd("") # set working directory
getwd()


library(enrichR)

dbs <- c("DisGeNET", "KEGG_2019_Human", "Reactome_2016", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018")

geneSets <- c("Autistic Disorder", "Autism Spectrum Disorders", "Schizophrenia", "Bipolar Disorder", "Intellectual Disability", "Major Depressive Disorder", 
              "Axon guidance", "Glutamatergic synapse", "Synaptic vesicle cycle", "axonogenesis (GO:0007409)", "axon guidance (GO:0007411)", 
              "nervous system development (GO:0007399)", "chemical synaptic transmission (GO:0007268)", "modulation of chemical synaptic transmission (GO:0050804)", 
              "glutamate receptor signaling pathway (GO:0007215)", "axon (GO:0030424)", "ionotropic glutamate receptor complex (GO:0008328)", 
              "Axon guidance Homo sapiens R-HSA-422475", "Transmission across Chemical Synapses Homo sapiens R-HSA-112315", 
              "Neurotransmitter Receptor Binding And Downstream Transmission In The  Postsynaptic Cell Homo sapiens R-HSA-112314", "Neuronal System Homo sapiens R-HSA-112316")

df <- data.frame(matrix(0, ncol = length(geneSets)+1, nrow = 1))
colnames(df) <- c("Permutations", geneSets)

path <- "" # set path to "zebrafish genes human homologs.csv" file, a list containing human homologs of all zebrafish genes 
genes <- read.csv(path, header=F)$V1 

for (i in 1:1000){
  x <- as.character(sample(genes, 5000)) # random sample 5000 genes

  enriched <- enrichr(x, dbs)
  temp_df <- Reduce(rbind, list(enriched[[1]]))
  
  y <- as.list( rep(0, length(geneSets)) )
  y[[1]] <- i
  
  for (k in 1:length(geneSets)) {
    if ( any( temp_df[[1]] == geneSets[[k]] ) == FALSE ) next
    y[[k+1]] <- temp_df[ temp_df[[1]] == geneSets[[k]], c(7) ]    # column 7 is odds ratio
  }
  
  df <- rbind(df,y)
}

df <- df[-1,]

write.csv(df, "Enrichr permutation.csv", row.names=FALSE)

