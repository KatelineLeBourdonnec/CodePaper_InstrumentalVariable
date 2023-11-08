########################### SCRIPT TABLE 2 ####################################
#####
##### This script is used to summarize the simulation results to reproduce 
##### Table 2, in order to obtain the biases and coverage rates for each scenario
##### in the binary scenarios. 
##### A similar script is available in the file "script_table1.R" to reproduce table 1. 
#####
##### This script requires the "table.R" function to be loaded.
#####
################################################################################



##### Load packages #####
library(here)

##### Define directory for sources #####
directory <- here() # directory of the project - can be changed by the user
setwd(directory)

##### Functions loading  #####
source("table.R")

##### concatenate simulation results:#####
#For logistic/substitution
setwd("./log_sub")
pattern <- c("simu_boots_binaire_2000_")
files <- list.files(pattern=pattern)
pattern <- "simu_boots_binaire_6000_"
files1 <- list.files(pattern=pattern)
pattern <- "simu_boots_binaire_20000_"
files2 <- list.files(pattern=pattern)
# For linear substitution
setwd("../linear_sub")
pattern <- "simu_boots_2000_LINEAR"
files3 <- list.files(pattern=pattern)
pattern <- "simu_boots_6000_LINEAR"
files4 <- list.files(pattern=pattern)
pattern <- "simu_boots_20000_LINEAR"
files5 <- list.files(pattern=pattern)
# For logistic/Residual
setwd("../log_res")
pattern <- "simu_boots_2sri_2000_"
files6 <- list.files(pattern=pattern)
pattern <- "simu_boots_2sri_6000_"
files7 <- list.files(pattern=pattern)
pattern <- "simu_boots_2sri_20000_"
files8 <- list.files(pattern=pattern)
files <- c(files,files1,files2,files3,files4,files5,files6,files7,files8)

##### Combine results : #####
# baseline coefficients :
coefN <- NULL
sdN <- NULL
coefIV <- NULL
se_corr_b <- NULL
# overtime coefficients :
coefN_inter <- NULL
sdN_inter <- NULL
coefIV_inter <- NULL
se_corr_int <- NULL
R0 <- NULL
F0 <- NULL
RB <- NULL
setwd("../log_sub")
r <- 0
for(f in files[1:4500]){
  r <- r+1
  load(f)
  
  coefN <- append(coefN,as.numeric(res[seq(6,length(res),25)]))
  sdN <- append(sdN,as.numeric(res[seq(7,length(res),25)]))
  
  coefIV <- append(coefIV,as.numeric(res[seq(9,length(res),25)]))
  se_corr <- res[seq(15,length(res),25)]
  se_corr_b[r] <- se_corr[[1]][3]
  
  coefN_inter <- append(coefN_inter,as.numeric(res[seq(5,length(res),25)]))
  sdN_inter <- append(sdN_inter,as.numeric(res[seq(8,length(res),25)]))
  
  coefIV_inter <- append(coefIV_inter,as.numeric(res[seq(10,length(res),25)]))
  se_corr_int[r] <- se_corr[[1]][4]
  
  R0 <- append(R0,as.numeric(res[seq(13,length(res),25)]))
}
setwd("../linear_sub")
for(f in files[4501:9000]){
  r <- r+1
  load(f)
  
  coefIV <- append(coefIV,as.numeric(res[seq(9,length(res),26)]))
  se_corr <- res[seq(16,length(res),26)]
  se_corr_b[r] <- se_corr[[1]][3]
  
  
  coefIV_inter <- append(coefIV_inter,as.numeric(res[seq(10,length(res),26)]))
  se_corr_int[r] <- se_corr[[1]][4]
  
  R0 <- append(R0,as.numeric(res[seq(13,length(res),26)]))
  F0 <- append(F0,as.numeric(res[seq(14,length(res),26)]))
  
}
setwd("../log_res")
for(f in files[9001:13500]){
  r <- r+1
  load(f)
  coefIV <- append(coefIV,as.numeric(res[seq(5,length(res),21)]))
  se_corr <- res[seq(11,length(res),21)]
  se_corr_b[r] <- se_corr[[1]][2]
  
  coefIV_inter <- append(coefIV_inter,as.numeric(res[seq(6,length(res),21)]))
  se_corr_int[r] <- se_corr[[1]][5]
  R0 <- append(R0,as.numeric(res[seq(9,length(res),21)]))
}

##### Create table ####
# for naive model & log_substitution : 
diff <- seq(1,4500,by=500)
tab <- rbind(table(diff,coefN,sdN),table(diff,coefN_inter,sdN_inter),table(diff,coefIV,se_corr_b),table(diff,coefIV_inter,se_corr_int))
# for linear_substition :  
diff <- seq(4501,9000,by=500)
tab <- rbind(tab,table(diff,coefIV,se_corr_b),table(diff,coefIV_inter,se_corr_int))
# for logistic_residual :
diff <- seq(9001,13500,by=500)
tab <- rbind(tab,table(diff,coefIV,se_corr_b),table(diff,coefIV_inter,se_corr_int))

##### final table : #####
setwd(directory)
setwd("tables")
tab <- rbind(c(rep(2,3),rep(3,3),rep(4,3)),tab)
row.names(tab) <- c("AlphaIV=","RB_naive","CR_naive","RB_naive_slope","CR_naive_slope",
                    "RB_log_sub","CR_log_sub","RB_log_sub_slope","CR_log_sub_slope",
                    "RB_lin_sub","CR_lin_sub","RB_lin_sub_slope","CR_lin_sub_slope",
                    "RB_log_res","CR_log_res","RB_log_res_slope","CR_log_res_slope")

colnames(tab) <- c(rep(c("2000","6000","20K"),3))
write.table(tab,file="Table2.txt")
