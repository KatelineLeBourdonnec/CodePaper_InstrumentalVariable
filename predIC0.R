##### Function to obtain confidence interval with correct standard error
#' @param x : object regression model
#' @param t : time
#' @param group : indicator of exposure
#' @param seCorrect : correct standard error
#' 
##### This function returns confidence interval for coefficients of "x" i.e. the regression model
predIC0 <- function(x,t,group, seCorrec){
  t1 <-ns(t, knots=c(4.19,11.71), Boundary.knots = c(0,17.12))[,1]
  t2 <- ns(t, knots=c(4.19,11.71), Boundary.knots = c(0,17.12))[,2]
  t3 <-ns(t, knots=c(4.19,11.71), Boundary.knots = c(0,17.12))[,3]
  
  tmp <- Wald2(x,pos=c(1,2,3,4,5,7,8,9),seCorrec=VarBeta, contrast=c(1,t1,t2,t3,group,t1*group,t2*group,t3*group))[c(1,2)]
  binf <- tmp[1]-1.96*tmp[2]
  bsup <- tmp[1]+1.96*tmp[2]
  return(c(t,group,tmp, binf, bsup))
}