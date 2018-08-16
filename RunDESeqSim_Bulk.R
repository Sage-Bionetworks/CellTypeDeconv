#Path to folder containing code 
setwd('E:/SageDocs/CellTypeDeconv/')

#load required packages 
require(DESeq)
require(monocle)
require(MASS)

#load codes 
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
GeneIn <- l$In_genes

#Simulate data for pure cell types 
NoSamp <- 10
PureData <- Generate.PureDF(l$M, l$Disp)
PureData$Gene <- GeneNames[GeneIn]

#Simulate data for random mixtures 
MixData <- Generate.Mix(l$M, l$Disp, NoSamp = NoSamp)
MixData$Dat$Gene <- GeneNames[GeneIn]
