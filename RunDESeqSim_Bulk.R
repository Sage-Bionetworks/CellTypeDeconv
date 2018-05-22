setwd('E:/SageDocs/CellTypeDeconv/')
source('DESeq_Helper.R')

#Load data 
Dat <- read.csv('BulkData/NeuronalBulk.csv')
GeneNames <- Dat$gene.symbol
Dat <- Dat[,2:dim(Dat)[2]]
Lab <- colnames(Dat)
Lab <- factor(gsub('.{1}$', '', Lab))
Lab <- as.numeric(Lab)

#Convert to counts 
Dat <- data.frame(Convert2Counts(Dat, Lab))

#Estimate parameters 
l <- Estimate.All.CellTypes(Dat, Lab)

#Simulate data for pure cell types 
NoSamp <- 10
PureData <- Generate.Pure(l$M, l$Disp, NoSamp = NoSamp)

#Simulate data for random mixtures 
MixData <- Generate.Mix(l$M, l$Disp, NoSamp = NoSamp)
