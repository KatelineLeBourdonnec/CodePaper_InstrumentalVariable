################# REPLICATION OF A SIMULATION SCENARIO ##################

# This code runs a scenario of simulations for a continuous exposure
#
# Techniques for binary estimation are given for: 
# The residual estimation technique is given in file 
# "CODE_BJ_Function_Simu_Log-Res.R". 
# The linear first-stage is given in file
# "CODE_BJ_Function_Simu_Linear-Sub.R"
# The logistic first-stage is given in file
# "CODE_BJ_Function_Simu_Log-Sub.R"


#This code was originally intended to be submitted to a computing server 
# (e.g., using slurm) for faster computations. Server arguments are put in
# comments here, and alternative values for rep and job parameters are given
# as an example. 


# The file is split in sequential parts indicated below with a "#####". 

# This script needs the sourcing of two functions:
# "geneData" used to simulate a dataset
# "bootstrap_param" used to correct the standard error

#The results of this script are output in folder "cont_log_sub". 


##### Load packages #####
library(MASS)
library(lcmm)
library(fmsb)
library(doParallel)
library(marqLevAlg)
library(mvtnorm)
library(here)

##### For reproducibility : #####
set.seed(560)


##### Define directory for sources #####
directory <- here() # directory of the project - can be changed by the user
setwd(directory)
##### Functions loading  #####
source("bootstrap_param_lin.R")
source("geneData.R")

##### Specification of the simulation scenario #####
# Values of the parameters and arguments used in our simulation.  
# The different scenarios are given in Supplementary Table S1. 
# Here is an example of setting for n=6000 subjects with a continuous 
# exposure variable and a strength of IV of alphaIV = 0.5

n <- 6000 # TO BE CHANGED BY THE USER FOR OTHER SAMPLE SIZES
nXiv <- 1 # number of IV variables
nXu <- 1 # number of unmeasured confounders U
moyXu<- 0 # mean of U
sdXu <- 1 # sd of U
k <- 6 # number of repeated visits
binary=F  # type of exposure
alpha <- c(3,0.5,1) # Vector of parameters for the exposure: 
# intercept, alphaIV, alphaU (unmeasured confounder). 
# The second value is to be changed for other scenarios of the manuscript
moyU0 <- 0 # mean of the random intercept
moyU1 <- 0 # mean of the random slope
matrixSdU0_U1 <- matrix(c(1,0,0,1),2,2) # variance covariance of the random effects
sdEpY <- 1 # SD of the measurement error
Beta <- c(5,1,1,1,1,1) # coefficients of the linear regression (2nd stage)
K <- 200 # number of bootstrap replicates for SE computation

# Make 500 draws to simulate (on computing server)
res <-  foreach(i=1:500, .combine='c', .packages=c('MASS','lcmm','doParallel','foreach','marqLevAlg','fmsb','ggplot2','mvtnorm')) %dopar%
  { 
##### Generate data ######
    # create the dataset
    Data_num <- geneData(nXu,nXiv,n,k,p,sdXu,moyXu, alpha, Beta, sdEpY,moyU0, moyU1,matrixSdU0_U1, sdTemps, pctMissing, binary)
    colnames(Data_num) <- c("ID","Temps","Xe",colnames(Data_num[,4:(3+nXiv+nXu)]),"Y")
    # the corresponding data with 1 row per subject
    data_sujet<- Data_num[seq(1,n*k,k),]
    data_sujet<- data_sujet[,c(-2,-length(Data_num))]
    
##### Estimation of the models ##### 
    
    # Naive model without Confounders
    # formula
    formulaY_naif0 <-  Y~Temps*(Xe)
    formulaY_naif0 <- as.formula(formulaY_naif0)
    # estimation
    hlme_N0 <- hlme(formulaY_naif0,random=~ 1+Temps, subject="ID",data=Data_num)
    summaryN0 <- summary(hlme_N0)
    CoefNaif0 <- hlme_N0$best[[3]]
    CoefNaifInter0 <- hlme_N0$best[[4]]
    sdNaif0 <- summaryN0[,2][3]
    sdNaifInter0 <- summaryN0[,2][4]
    
    # true model (i.e., as generated with confounders)
    # Formula
    formulaY_T <- paste0(c("Y~Temps*(Xe+",paste(colnames(Data_num)[5],collapse ="+"),")"),sep="")
    formulaY_T <- as.formula(formulaY_T)
    # Estimation
    hlme_T <- hlme(formulaY_T,random=~ 1+Temps, subject="ID",data=Data_num)
    summaryT <- summary(hlme_T)
    CoefTrue <- hlme_T$best[[3]]
    CoefTrueInter <-hlme_T$best[[3+nXu+1]]
    sdTrue <- summaryT[,2][3]
    sdTrueInter <- summaryT[,2][3+nXu+1]
    
    
    # 1st stage model
    lm1 <- lm(Xe~Xiv,data=data_sujet)
    # prediction computation
    Data_num$XeHat <- as.matrix(cbind(rep(1,length(Data_num$ID)),Data_num[,4]))%*%lm1$coefficients
    
    #  R2 statistic computation
    r2_1 <- cor(data_sujet$Xe,predict(lm1))**2
    F_1 <- (r2_1/(nXiv))/((1-r2_1)/lm1$df.residual)
    
    # 2nd stage model 
    # formula
    formulaY_2stage <- Y~Temps*(XeHat)
    formulaY_2stage <- as.formula(formulaY_2stage)
    # estimation
    hlme_XeHat <- hlme(formulaY_2stage,random=~ 1+Temps, subject="ID",data=Data_num)
    summary <- summary(hlme_XeHat)
    Coef2stageM <- hlme_XeHat$best[[3]]
    Coef2stageInter <- hlme_XeHat$best[[3+1]]
    sd2stageM <- summary[,2][[3]]
    sd2stageInter <- summary[,2][[3+1]]
    
##### Correction standard error for the 2 stage model #####
    
    l<- sum(hlme_XeHat$idea0)
    B <- matrix(0,nrow=l,ncol=l)
    B[upper.tri(B,diag=TRUE)] <- hlme_XeHat$best[(4+1):(length(hlme_XeHat$best)-1)] 
    B[lower.tri(B,diag=FALSE)] <- t(B)[lower.tri(B,diag=FALSE)]
    sigma2 <- (hlme_XeHat$best[length(hlme_XeHat$best)])**2
    
    
    Vbeta <- matrix(data=0,nrow=(length(hlme_XeHat$best)-4),ncol=(length(hlme_XeHat$best)-4))
    Vrond <- matrix(data=0,nrow=(length(hlme_XeHat$best)-4),ncol=(length(hlme_XeHat$best)-4))
    
    for(i in unique(Data_num$ID)){ 
      
      #Variance of Yi
      Z <- cbind("int" =rep(1,k), "t"=Data_num$Temps[which(Data_num$ID==i)]) # matrice des variables associC)es effets alC)atoires
      Vi <- as.matrix(Z)%*%B%*%t(as.matrix(Z))+diag(sigma2,nrow =k)
      Vi1 <- solve(Vi)
      
      #Mat var exp
      
      X <- cbind(rep(1,k),Data_num[which(Data_num$ID==i),c(2,(length(Data_num)))],Data_num[which(Data_num$ID==i),c(2)]*Data_num[which(Data_num$ID==i),(length(Data_num))])
      colnames(X) <- c("int",colnames(X[2:3]),"tXhat")
      Var <- as.matrix(t(X))%*%Vi1%*%as.matrix(X)
      
      # covariance : 
      
      Xtrue <- cbind("int"=rep(1,k),Data_num[which(Data_num$ID==i),c(2:3)],Data_num[which(Data_num$ID==i),2]*Data_num[which(Data_num$ID==i),3])
      colnames(Xtrue) <- c("int",colnames(Xtrue[2:3]),"tXe")
      M <- as.matrix(Xtrue)%*%hlme_XeHat$best[1:4]
      Ma <- Data_num$Y[which(Data_num$ID==i)] - M

      covY <- Ma %*% t(Ma)
      
      VarCorr <- (as.matrix(t(X))%*%Vi1%*%covY%*%Vi1%*%as.matrix(X))
      
      Vbeta <- Vbeta+ Var
      Vrond <- Vrond + VarCorr
      
      
    }
    VarBeta <- solve(Vbeta)
    varcovNew0 <- sqrt(diag(VarBeta)) 
    summary(hlme_XeHat)
    
    
    var_corrig <- VarBeta %*% Vrond %*% VarBeta
    se_corr <- sqrt(diag(var_corrig))
    
##### Bootstrap Procedure #####
    resu <- bootstrap_param(K,data_sujet,Data_num)
    
    Bchapmoy <- mean(resu[[1]])
    Bchapintmoy <- mean(resu[[2]])
    # variance
    VarIntra <- mean(resu[[3]])
    VarInter <- sum((resu[[1]]-Bchapmoy)**2)/K
    VarTot <- VarIntra + VarInter 
    
    VarIntraInt <- mean(resu[[4]])
    VarInterInt <- sum((resu[[2]]-Bchapintmoy)**2)/K
    VarTotInt <- VarIntraInt + VarInterInt 
    
    ###confidence interval
    
    borneinf <- Bchapmoy - 1.96*(sqrt(VarTot))
    bornesup <- Bchapmoy + 1.96*(sqrt(VarTot))
    
    borneinfInt <- Bchapintmoy - 1.96*(sqrt(VarTotInt))
    bornesupInt <- Bchapintmoy + 1.96*(sqrt(VarTotInt))
    
    
    # Coverage rate of the 95% confidence interval
    TC <- ifelse(borneinf <= Beta[2] & bornesup > Beta[2],1,0)
    TCint <- ifelse(borneinf <= Beta[4] & bornesup > Beta[4],1,0)
    
    d <- list(CoefTrue,CoefTrueInter,sdTrue, sdTrueInter,CoefNaifInter0,CoefNaif0,sdNaif0,sdNaifInter0,Coef2stageM,Coef2stageInter,sd2stageM,sd2stageInter,r2_1,F_1,var_corrig,se_corr,TC,TCint,VarIntra,VarInter,VarTot,VarIntraInt,VarInterInt,VarTotInt,resu[[1]],resu[[2]])
    
}
##### Output #####

save(res,file="simu_boots_6000_05.Rdata")





