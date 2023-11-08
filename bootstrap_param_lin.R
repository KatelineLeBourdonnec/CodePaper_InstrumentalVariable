#### Function for bootstrap for linear first stage: #### 

#The method used is a two-stage methodology, which is why a parametric bootstrap 
#method is necessary to obtain a valid confidence interval.


#' @param K number of bootstrap
#' @param data_suj cross-sectionnal dataset (one row by subject)
#' @param data_num longitudinal dataset (one row by visit for each subject)

bootstrap_param <- function(K,data_suj,data_num){
  
  # STEP 1 : Get Xhat by linear or logistic regression (depending on the nature of the exposure)
  # here it's for linear regression
  lmfit <- lm(Xe ~ Xiv,data=data_suj)
  coefAlphaIV <- cov(data_suj$Xiv,data_suj$Xe)/var(data_suj$Xiv)
  coefAlpha0 <- mean(data_suj$Xe)-coefAlphaIV*mean(data_suj$Xiv)
  
  Sy_X <- sqrt(sum((data_suj$Xe-predict(lmfit))**2)/(length(data_suj$ID)-2))
  
  Siv <- Sy_X/sqrt(sum((data_suj$Xiv-mean(data_suj$Xiv))**2))
  S0 <- Sy_X*sqrt(sum(data_suj$Xiv**2)/(n*sum((data_suj$Xiv-mean(data_suj$Xiv))**2)))
  
  # STEP 2 : Generate data according to Xhat value 
  alphaB0 <- rnorm(n=K,coefAlpha0,S0)
  alphaBIV <- rnorm(n=K,coefAlphaIV,Siv)
  alphaB  <- cbind(alphaB0,alphaBIV)
  
  # STEP 3 & STEP 4  : Correction of the standard error by bootstrapping
  coefB <- NULL
  coefintB <- NULL
  seB <- NULL
  seintB <- NULL
  for(i in 1:K){
    
    data_num$XeHat <- cbind(rep(1,length(data_num$ID)),data_num$Xiv)%*%alphaB[i,] # STEP 3
    lmXeHat <- hlme(Y ~ Temps*XeHat,random=~ 1+Temps, subject="ID",data=data_num) # STEP 4
    
    coefB <-append(coefB,lmXeHat$best[3])
    coefintB <- append(coefintB, lmXeHat$best[4])
    
    # CORRECTION STANDARD ERROR
    
    l<- sum(lmXeHat$idea0)
    B <- matrix(0,nrow=l,ncol=l)
    B[upper.tri(B,diag=TRUE)] <- lmXeHat$best[(4+1):(length(lmXeHat$best)-1)] 
    B[lower.tri(B,diag=FALSE)] <- t(B)[lower.tri(B,diag=FALSE)]
    sigma2 <- (lmXeHat$best[length(lmXeHat$best)])**2
    
    
    Vbeta <- matrix(data=0,nrow=(length(lmXeHat$best)-4),ncol=(length(lmXeHat$best)-4))
    Vrond <- matrix(data=0,nrow=(length(lmXeHat$best)-4),ncol=(length(lmXeHat$best)-4))
    
    for(j in unique(data_num$ID)){ 
      
      #Variance of Yi
      Z <- cbind("int" =rep(1,k), "t"=data_num$Temps[which(data_num$ID==j)]) # matrice des variables associC)es effets alC)atoires
      Vi <- as.matrix(Z)%*%B%*%t(as.matrix(Z))+diag(sigma2,nrow =k)
      Vi1 <- solve(Vi)
      
      #Matr var exp
      
      X <- cbind(rep(1,k),data_num[which(data_num$ID==j),c(2,(length(data_num)))],data_num[which(data_num$ID==j),c(2)]*data_num[which(data_num$ID==j),(length(data_num))])
      colnames(X) <- c("int",colnames(X[2:3]),"tXhat")
      Var <- as.matrix(t(X))%*%Vi1%*%as.matrix(X)
      
      # covariance : 
      
      Xtrue <- cbind("int"=rep(1,k),data_num[which(data_num$ID==j),c(2:3)],data_num[which(data_num$ID==j),2]*data_num[which(data_num$ID==j),3])
      colnames(Xtrue) <- c("int",colnames(Xtrue[2:3]),"tXe")
      M <- as.matrix(Xtrue)%*%lmXeHat$best[1:4]
      Ma <- data_num$Y[which(data_num$ID==j)] - M
      # VY <- cov(Ma,t(Ma))
      
      covY <- Ma %*% t(Ma)
      VarCorr <- (as.matrix(t(X))%*%Vi1%*%covY%*%Vi1%*%as.matrix(X))
      
      Vbeta <- Vbeta+ Var
      Vrond <- Vrond + VarCorr
      
      
    }
    VarBeta <- solve(Vbeta)
    varcovNew0 <- sqrt(diag(VarBeta)) 
    summary(lmXeHat)
    
    
    var_corrig <- VarBeta %*% Vrond %*% VarBeta
    se_corr <- sqrt(diag(var_corrig))
    
    
    seB <- append(seB, se_corr[3])
    seintB <- append(seintB,se_corr[4])
    
    
  }
  estimationB <- list(coefB, coefintB,seB,seintB)
  return(estimationB)
}
