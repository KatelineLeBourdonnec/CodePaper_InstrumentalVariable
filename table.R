################################### TABLE ######################################
#
# Function used to create table in order to obtain :
#         Relative bias 
#         coverage rate 
# 
#
#' @param diff 
table <- function(diff,coef,sd){
  j <- 1
  RB <- NULL
  TxC <- NULL
  
  for(i in diff){
    RB1 <- mean((coef[i:c(i+(diff[j+1]-diff[j]-1))]-1)*100)
    RB <- cbind(RB,RB1)
    borninf <- coef[i:(i+(diff[j+1]-diff[j]-1))]-1.96*sd[i:(i+(diff[j+1]-diff[j]-1))]
    bornsup <- coef[i:(i+(diff[j+1]-diff[j]-1))]+1.96*sd[i:(i+(diff[j+1]-diff[j]-1))]
    TxC <- cbind(TxC,length(which(borninf<= 1 & bornsup>=1))/c((diff[j+1]-diff[j])))
  }
  j <- j+1
  tab1 <- rbind(RB,TxC)
  return(tab1)}
