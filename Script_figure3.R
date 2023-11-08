########################## SCRIPT APPLICATION ##################################
################################################################################
################################################################################
###### This script replicates the analysis of the application section. 
###### The script requires the loading of two functions: 
###### predIC0 and Wald2. 
######
###### The real data were not available for sharing. The script thus uses 
###### simulated data (data_long.txt) that mimic as best as possible 
###### the actual application data. 
######  For simulating the data, we used the parameters from the estimated IV model  
###### plus the coefficients of the naive model for the adjustment variables. 
###### We added an arbitrary random unmeasured confounder as this information  
###### was not available in the real data. Consequently,
###### the results obtained for the application in the manuscript and 
###### in this replication script will not be exactly the same. 
######
###### In this script, different models are run to compare the results
###### of different methods.  
###### Trajectories (corresponding to Figure 3) are plotted in file
###### "Figure3.pdf"
################################################################################


##### Loading packages #####
library(ggplot2)
library(lcmm)
library(splines)
library(here)

##### Define directory #####
directory <- here()
setwd(directory)

##### Loading data & functions #####
# Data
read.table("data_long.txt",sep="")
# functions
source("predIC0.R")
source("Wald2.R")

##### Estimation ######
# Naive regression : 
regNaif1Xu <- hlme(Y ~ ns(Temps, knots=c(4.19,11.71), Boundary.knots = c(0,17.12))*Xe+Primo, random = ~1+ns(Temps, knots=c(4.19,11.71), Boundary.knots = c(0,17.12)), subject = "ID", data=data_long, verbose=F)

# Naive regression with adjustment :
regNaif1ajusXu <- hlme(Y ~ ns(Temps, knots=c(4.19,11.71), Boundary.knots = c(0,17.12))*(Xe+AGE75+SEX+niv2+niv3+niv4+niv5)+Primo, random = ~1+ns(Temps, knots=c(4.19,11.71), Boundary.knots = c(0,17.12)), subject = "ID", data=data_long, verbose=F, maxiter = 200)

# IV method with spline : 
reg_pred_spline1 <- hlme(Y ~ ns(Temps, knots=c(4.19,11.71), Boundary.knots = c(0,17))*predGLM+Primo,random=~ 1+ns(Temps,  knots=c(4.19,11.71), Boundary.knots = c(0,17)),subject = "ID", data= data_long, verbose=F)

##### Correction for confidence interval #####

nbEA <- reg_pred_spline1$N[3]+1
coef2SLS <- reg_pred_spline1$best[1:(length(reg_pred_spline1$best)-nbEA)]
binit <- regNaif1Xu$best
binit[1:(length(reg_pred_spline1$best)-nbEA)] <- coef2SLS

regMM_N <- hlme(Y ~ ns(Temps, knots=c(4.19,11.71), Boundary.knots = c(0,17.12))*Xe+Primo,random=~ 1+ns(Temps, knots=c(4.19,11.71), Boundary.knots = c(0,17.12)), subject = "ID", data= data_long, B = binit, posfix = 1:reg_pred_spline1$N[2], verbose=F) 

l<- sum(regNaif1Xu$idea0)
B <- matrix(0,nrow=l,ncol=l)
B[upper.tri(B,diag=TRUE)] <- regNaif1Xu$best[(length(regNaif1Xu$best)-nbEA+1):(length(regNaif1Xu$best)-1)] 
B[lower.tri(B,diag=FALSE)] <- t(B)[lower.tri(B,diag=FALSE)]
sigma2 <- (regNaif1Xu$best[length(regNaif1Xu$best)])**2


Vbeta <- matrix(data=0,nrow=regNaif1Xu$N[2],ncol=regNaif1Xu$N[2])

random <-~ 1+ns(Temps, knots=c(4.19,11.71), Boundary.knots = c(0,17.12))
for(i in unique(data_long$ID)){ 
  #Variance Yi
  Z <- model.matrix(random, data_long)
  Z <- cbind(Z, data_long$ID)
  Z <- matrix(Z[which(Z[,dim(Z)[2]]==i),-dim(Z)[2]],ncol=ncol(B))
  # Matrix random effect
  Vi <-Z%*%B%*%t(Z)+diag(sigma2,nrow =length(which(data_long$ID==i)))
  Vi1 <- solve(Vi)
  
  #Matrix var expli 
  formula <- Y ~ ns(Temps, knots=c(4.19,11.71), Boundary.knots = c(0,17.12))*predGLM +Primo
  options(na.action='na.pass') 
  X <- model.matrix(formula,data_long)
  X <- cbind(X, data_long$ID)
  X <- matrix(X[which(X[,dim(X)[2]]==i),-dim(X)[2]], ncol=reg_pred_spline1$N[2])
  
  Var <- t(X)%*%Vi1%*%X
  Vbeta <- Vbeta+ Var
}
VarBeta <- solve(Vbeta)
varcovNew <- sqrt(diag(VarBeta)) 

##### Data preparation for graphics #######

t <- seq(0,10.9,by=0.11)
x <- reg_pred_spline1

group = 0 # No diabetes
group = 1 # Diabetes
Primo <- 0
# CI for diabetics
diabeteCI <- t(sapply(t,function(t){predIC0(x,t,group=1,seCorrec=VarBeta)}))
colnames(diabeteCI) <- c("time","group","coef","se","binf","bsup")

# Ci for non diabetics
nondiabeteCI <- t(sapply(t,function(t){predIC0(x,t,group=0,seCorrec=VarBeta)}))
colnames(nondiabeteCI) <- c("time","group","coef","se","binf","bsup")

# Prediction for non diabetics : 
datnew <- data.frame("Temps"=seq(0,11,length=100))
datnew$predGLM <- 0
datnew$Xe <- 0
datnew$Primo <- 0
datnew$niv2 <- 0
datnew$niv3 <- 1
datnew$niv4 <- 0
datnew$niv5 <- 0
datnew$SEX <- 1
datnew$AGE75 <- mean(data_long$AGE75)

nondiabetenaif <- predictY(regNaif1Xu, newdata=datnew, var.time="Temps", draws=T)
nondiabetenaifConf <- predictY(regNaif1ajusXu, newdata=datnew, var.time="Temps", draws=T)
nondiabete_spline <- predictY(reg_pred_spline1, newdata=datnew, var.time="Temps", draws=T)

# Prediction for diabetics : 
datnew$predGLM <- 1
datnew$Xe <- 1

diabetenaif <- predictY(regNaif1Xu, newdata=datnew, var.time="Temps", draws=T)
diabetenaifConf <- predictY(regNaif1ajusXu, newdata=datnew, var.time="Temps", draws=T)
diabete_spline <- predictY(reg_pred_spline1, newdata=datnew, var.time="Temps", draws=T)

# Correction for CI
diabeteCorr <- diabete_spline
diabeteCorr$pred[,2] <- diabeteCI[,5]
diabeteCorr$pred[,3] <- diabeteCI[,6]
nondiabeteCorr <- nondiabete_spline
nondiabeteCorr$pred[,2] <- nondiabeteCI[,5]
nondiabeteCorr$pred[,3] <- nondiabeteCI[,6]

##### Final graph #####

setwd("Figure")

pdf(file ="Figure3.pdf",width=14)
par(mfrow=c(1,3))
# method naive
plot(diabetenaif, lwd=c(4,1), type="l", col="darkslategray4", ylim=c(30,50), xlab="Time in year",ylab="IST score",bty="l", legend=NULL, shades = TRUE,  main="Naive method")
a <- plot(nondiabetenaif, add=TRUE, col="black", lwd=c(4,1), shades=TRUE)
# method naive with ajustement
plot(diabetenaifConf, lwd=c(4,1), type="l", col="darkslategray4", ylim=c(30,50), xlab="Time in year",ylab="IST score",bty="l", legend=NULL, shades = TRUE,  main="Adjusted naive method")
b <- plot(nondiabetenaifConf, add=TRUE, col="black", lwd=c(4,1), shades=TRUE)
# method IV with correct CI
plot(diabeteCorr, lwd=c(4,1), type="l", col="darkslategray4", ylim=c(30,50), xlab="Time in year",ylab="IST score",bty="l", legend=NULL, shades = TRUE,  main="IV method")
c <- plot(nondiabeteCorr, add=TRUE, col="black", lwd=c(4,1), shades=TRUE)
legend(x="topleft", bty="n", ncol=1, lty=c(1,1), col=c("darkslategray4","black"), legend=c("Diabetic","No diabetic"), lwd=c(2,2)) 
dev.off()




