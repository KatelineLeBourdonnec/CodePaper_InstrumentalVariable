########################### SCRIPT TABLE 1 ####################################
#####
##### This script is used to summarize the simulation results and reproduce 
##### Table 1, in order to obtain the biases and coverage rates for each scenario
##### in the continuous scenario. 
##### A similar script is available in the file "script_table2.R" to reproduce table 2. 
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
setwd("./cont_log_sub")
pattern <- "simu_boots_2000.R"
files_c_1 <- list.files(pattern=pattern)
pattern <- "simu_boots_6000_05.R"
files_c_2 <- list.files(pattern=pattern)
pattern <- "simu_boots_05_20000_"
files_c_3 <- list.files(pattern = pattern)
pattern <- "simu_boots_2000_4"
files_c_4 <- list.files(pattern=pattern)
pattern <- "simu_boots_6000.R"
files_c_5 <- list.files(pattern=pattern)
pattern <- "simu_boots_20000_"
files_c_6 <- list.files(pattern = pattern)
files_c <- c(files_c_1,files_c_2,files_c_3,files_c_4,files_c_5,files_c_6)

##### Combine results : #####
# baseline coefficients :
coefN <- NULL
sdN <- NULL
coefIV <- NULL
se_corr_b <- NULL
# coefficients for association with time  :
coefN_inter <- NULL
sdN_inter <- NULL
coefIV_inter <- NULL
se_corr_int <- NULL
R0 <- NULL
F0 <- NULL
r <- 0
for(f in files_c){
  r <- r+1
  load(f)
  
  coefN <- append(coefN,as.numeric(res[seq(6,length(res),26)]))
  sdN <- append(sdN,as.numeric(res[seq(7,length(res),26)]))
  
  coefN_inter <- append(coefN_inter,as.numeric(res[seq(5,length(res),26)]))
  sdN_inter <- append(sdN_inter,as.numeric(res[seq(8,length(res),26)]))
  
  coefIV <- append(coefIV,as.numeric(res[seq(9,length(res),26)]))
  se_corr <- res[seq(16,length(res),26)]
  se_corr_b <- append(se_corr_b,lapply(se_corr, function(obj) as.numeric(obj[3])))
  
  
  coefIV_inter <- append(coefIV_inter,as.numeric(res[seq(10,length(res),26)]))
  se_corr_int <- append(se_corr_int,lapply(se_corr, function(obj) as.numeric(obj[4])))
  
  R0 <- append(R0,as.numeric(res[seq(13,length(res),26)]))
  F0 <- append(F0,as.numeric(res[seq(14,length(res),26)]))
}

##### Create table ####
diff <- c(1,501,601,1101,1601,2101)
tab <- rbind(table(diff,coefN,sdN),table(diff,coefN_inter,sdN_inter))
tabIV <- rbind(table(diff,coefIV,unlist(se_corr_b)),table(diff,coefIV_inter,unlist(se_corr_int)))

##### final table : #####
setwd(directory)
setwd("Tables")
tab <- rbind(c(rep(0.5,3),rep(1,3)),tab,tabIV)
row.names(tab) <- c("AlphaIV=","RB_naive","CR_naive","RB_naive_slope","CR_naive_slope",
                    "RB_cont_log_sub","CR_cont_log_sub","RB_cont_log_sub_slope","CR_cont_log_sub_slope")

colnames(tab) <- c(rep(c("2000","6000","20K"),2))
write.table(tab,file="Table1.txt")
