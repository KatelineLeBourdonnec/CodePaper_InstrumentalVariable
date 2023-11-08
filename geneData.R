
### Function "geneData" for generate data simulated ###

#' @param nXU number of unmeasured confounder
#' @param nXiv number of instrumental variable
#' @param n sample size
#' @param k number of visit
#' @param sdXu standard error XU
#' @param moyXu mean XU 
#' @param alpha Coefficient vector for the first stage - linear or logistic regression for cross-sectionnal data
#' @param Beta  Coefficient vector for the second stage - linear mixed regression for longitudinal data
#' @param sdEpY error term for the linear mixed model
#' @param moyU0 mean of baseline random effect
#' @param moyU1 mean slope for random effect
#' @param matrixSdU0_U1 matrix of correlation for the random effect
#' @param sdTemps sd error of time
#' @param pctMissing missing value
#' @param binary boolean variable for exposure variable : methodology adapted for binary (binary=TRUE) or continuous (binary=FALSE) exposure



geneData <- function(nXu,nXiv,n,k,p,sdXu,moyXu, alpha, Beta, sdEpY,moyU0, moyU1,matrixSdU0_U1, sdTemps, pctMissing, binary=F){
  # max row in longitudinal base :
  N <- n*k 
  
  # Generate instrumental variable, unmeasured confounder, epsilon error 
  Xiv <- replicate(nXiv,rnorm(n), simplify=T)
  Xu <- replicate(nXu, rnorm(n,moyXu, sdXu),simplify=T)
  EpsX <- rnorm(n) # error for the first step
  
  # Data generation if continuous exposure to create a dataset
  if (binary==F){
    XforXe <- cbind(rep(1,n),Xiv,Xu)
    Xe <- XforXe%*%alpha + EpsX
    dataXe <- data.frame("ID"=1:n,Xe,Xiv,Xu,EpsX)
    
    t <- NULL
    for(i in 1:k){
      tps <- rnorm(n,i,sd=0.05)
      t <- append(t,tps)
    }
    
    base_obs <- cbind("ID" = rep(1:n,k),"int"=rep(1,n*k), "t"=t)
    base_obs <-base_obs[order(base_obs[,c("ID")]),]
    base_obs <- merge(base_obs,dataXe[,-length(dataXe)], by='ID')
    base_obs[,(length(base_obs)+1):(length(base_obs)+nXiv+nXu+1)] <- base_obs[,4:length(base_obs)]*base_obs[,"t"]
    
    U <- mvrnorm(n, c(moyU0,moyU1), Sigma =matrixSdU0_U1)
    
    EpsY <- rnorm(nrow(base_obs), mean=0, sd=sdEpY)
    
    base_obs <- merge(base_obs,cbind("ID"=rep(1:n),U),by="ID")
    base_obs$y <- as.matrix(base_obs[,c(2,3,4,(5+nXiv):(4+nXiv+nXu+1),(5+2*nXiv+nXu+1):(5+2*nXiv+2*nXu))])%*%Beta + base_obs[,length(base_obs)-1] + base_obs[,length(base_obs)]*base_obs$t + EpsY 
    
    # Data generation if binary exposure to create a dataset
  }else{
    
    XforXe <- cbind(rep(1,n),Xiv,Xu)
    p <- 1/(1+exp(-(XforXe%*%alpha)))
    Xe <- rbinom(n,1,p)
    dataXe <- data.frame("ID"=1:n,Xe,Xiv,Xu)
    
    t <- NULL
    for(i in 1:k){
      tps <- rnorm(n,i,sd=0.05)
      t <- append(t,tps)
    }
    
    base_obs <- cbind("ID" = rep(1:n,k),"int"=rep(1,n*k), "t"=t)
    base_obs <-base_obs[order(base_obs[,c("ID")]),]
    base_obs <- merge(base_obs,dataXe, by='ID')
    base_obs[,(length(base_obs)+1):(length(base_obs)+nXiv+nXu+1)] <- base_obs[,4:length(base_obs)]*base_obs[,"t"]
    
    U <- mvrnorm(n, c(moyU0,moyU1), Sigma =matrixSdU0_U1)
    
    EpsY <- rnorm(nrow(base_obs), mean=0, sd=sdEpY)
    
    base_obs <- merge(base_obs,cbind("ID"=rep(1:n),U),by="ID")
    base_obs$y <- as.matrix(base_obs[,c(2,3,4,(5+nXiv),(6+1),9)])%*%Beta + base_obs[,length(base_obs)-1] + base_obs[,length(base_obs)]*base_obs$t + EpsY 
    
  }
  
  return(base_obs[,-c(2,(4+nXiv+nXu+1):(length(base_obs)-1))])
  
}
