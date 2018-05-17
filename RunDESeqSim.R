setwd('E:/SageDocs/CellTypeDeconv/')
library(R.matlab)
source('DESeq_Helper.R')

#loading data
temp <- readMat('E:/UWPhDWork/SCS/DropBarcoding/CellClassificationData/10xPooled/Pooled400_data.mat')
Dat <- data.frame(as.matrix(temp$Dat))
Lab <- temp$Lab

#Estimate parameters 
l <- Estimate.All.CellTypes(Dat, Lab)

#Simulate data for pure cell types 
NoSamp <- 10
PureData <- Generate.Pure(l$M, l$Disp, NoSamp = NoSamp)

#Simulate data for random mixtures 
MixData <- Generate.Mix(l$M, l$Disp, NoSamp = NoSamp)