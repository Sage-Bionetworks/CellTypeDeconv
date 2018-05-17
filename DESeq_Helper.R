Estimate.Parameters <- function(Dat, Lab){
  
  library(DESeq)
  cds = newCountDataSet( Dat, as.factor(as.character(Lab)) )
  cds@phenoData@data$sizeFactor <- rep(1,length(Lab))
  cds = estimateDispersions( cds )
  
  
  l <- list()
  
  temp <- fitInfo(cds)
  l$M <- rowMeans(Dat)
  l$Disp <- temp$fittedDispEsts
  
  return(l)
  
}


Estimate.All.CellTypes <- function(Dat, Lab, Thresh = 1){
  
  In <- which(rowSums(Dat) >= Thresh)
  Dat <- Dat[In,]
  
  NoTypes <- length(unique(Lab))
  
  M <- matrix(rep(0,length(In)*NoTypes),nrow=length(In),ncol=NoTypes)

  l <- list()
  
  Main <- Estimate.Parameters(Dat, Lab)
  
  for (i in 1:NoTypes){
    
    In <- which(Lab==i)
    tmp <- Estimate.Parameters(Dat[,In], Lab[In])
    l[[i]] <- tmp
    M[,i] <- tmp$M
    
  }
  
  l2 <- list()
  l2$l <- l 
  l2$Main <- Main
  l2$M <- M 
  l2$Disp <- Main$Disp
  
  return(l2)
}

Gen.Data.Given.Params <- function(M, Disp, NoSamp){
  library(MASS)
  
  Dat <- replicate(NoSamp,rnegbin(M, theta = 1/Disp))
  
  return(Dat)
}



Generate.Params.Given.Mix <- function(M, w, Disp){
  
  l <- list()
  
  w_norm <- w/sum(w)
  M_new <- M %*% w 
  l$M <- M_new
  l$Disp <- Disp
  
  return(l)
  
}


Generate.Pure <- function(M, Disp, NoSamp){
  
  s <- dim(M)
  
  l <- list()
  
  for (i in 1:s[2]){
    l[[i]] <- Gen.Data.Given.Params(M[,i], Disp, NoSamp)
    
  }
  
  return(l)
}


Generate.Mix <- function(M, Disp, NoSamp, W = NULL){
  
  s <- dim(M)
  
  if (is.null(W)){
    W <- replicate(NoSamp, runif(s[2]))
  }
  
  s2 <- dim(W)
  
  #normalize coloumns of W 
  
  
  l <- list()
  
  Dat <- matrix(rep(0,s[1]*s2[2]),nrow=s[1],ncol=s2[2])
  
  for (i in 1:s2[2]){
    W[,i] <- W[,i]/sum(W[,i])
    tmp <- Generate.Params.Given.Mix(M, W[,i], Disp)
    Dat[,i] <- Gen.Data.Given.Params(tmp$M, tmp$Disp, 1)
  }
  
  l$Dat <- Dat 
  l$W <- W
  
  return(l)
  
}

