#### Function for Wald statistic test

#' @param Mod : regression model
#' @param pos : position of coefficient
#' @param name : NULL (expect if you want to rename columns of the results)
#' @param seCorrec : standard error
Wald2 <- function (Mod, pos = NULL, contrasts = NULL, name = NULL,seCorrec) 
{
  value <- NULL
  if (!(class(Mod) %in% c("hlme", "lcmm", "multlcmm", 
                          "Jointlcmm"))) 
    stop("applies to \"hlme\" or \"lcmm\" or \"multlcmm\" or \"Jointlcmm\" objects only")
  if (inherits(Mod, "hlme") | inherits(Mod, "lcmm")) {
    nea <- sum(Mod$idea)
    nef <- Mod$N[2]
    nvc <- Mod$N[3]
    nprob <- Mod$N[1]
    idiag <- ifelse(Mod$idiag == 1, TRUE, FALSE)
  }
  if (inherits(Mod, "multlcmm")) {

    nea <- sum(Mod$idea0)
    nef <- Mod$N[3]
    nvc <- Mod$N[4]
    nprob <- 0
    idiag <- ifelse(Mod$idiag == 1, TRUE, FALSE)
  }
  if (inherits(Mod, "Jointlcmm")) {
    nea <- sum(Mod$idea)
    nef <- Mod$N[4]
    nvc <- Mod$N[5]
    nprob <- sum(Mod$N[1:3])
    idiag <- ifelse(Mod$idiag == 1, TRUE, FALSE)
  }
  if (nvc > 0) {
    debut <- nprob + nef + 1
    fin <- nprob + nef + nvc
    cholesky <- Mod$cholesky
    if (inherits(Mod, "multlcmm")) 
      cholesky[1] <- NA
    if (isTRUE(idiag)) 
      cholesky[setdiff(1:(nea * (nea + 1)/2), 1:nea * (1:nea + 
                                                         1)/2)] <- NA
    Mod$best[debut:fin] <- na.omit(cholesky)
  }
  l <- length(Mod$best)
  V <- matrix(0, nrow = l, ncol = l)
  V[upper.tri(V, diag = TRUE)] <- Mod$V
  V[lower.tri(V, diag = FALSE)] <- t(V)[lower.tri(V, diag = FALSE)]
  if (is.null(pos)) {
    stop("pos must be specified")
  }
  else {
    if (!is.vector(pos)) 
      stop("Error : pos must be a numeric vector")
    Mat <- matrix(0, nrow = length(pos), ncol = length(pos))
    #Mat <- V[pos, pos]
    Mat <- seCorrec[pos,pos]
    Vect <- Mod$best[pos]
    
    if (is.null(contrasts)) {
      if (!is.null(value)) {
        if (!is.vector(value)) 
          stop("Error : value must be a numeric vector")
        if (length(value) != length(pos)) 
          stop("value must have the same length as the vector pos")
        Vect <- Mod$best[pos] - value
      }
      Wald <- t(Vect) %*% solve(Mat) %*% Vect
      ddl <- length(pos)
      p_value <- 1 - pchisq(Wald, df = ddl)
      Results <- matrix(NA, nrow = 1, ncol = 2)
      colnames(Results) <- c("Wald Test", "p_value")
      if (is.null(name)) {
        if (!is.null(value)) {
          rownames(Results) <- paste(names(Mod$best[pos]), 
                                     " = ", value, collapse = " and ", 
                                     sep = "")
        }
        else {
          rownames(Results) <- paste(paste(names(Mod$best[pos]), 
                                           collapse = " = "), "= 0")
        }
      }
      else {
        rownames(Results) <- name
      }
      Results[, 1] <- round(Wald, 5)
      Results[, 2] <- round(p_value, 5)
    }
    else {
      if (length(contrasts) != length(pos)) {
        stop("contrasts must have the same length as the vector pos")
      }
      if (sum(abs(contrasts)) == 0) {
        stop("The absolute value of the sum of contratsts components must be different from 0")
      }
      Scalaire <- sum(Vect * contrasts)
      if (!is.null(value)) {
        if (!is.vector(value)) 
          stop("value must be a numeric vector")
        if (length(value) != 1) 
          stop("value must be a vector with a unique argument")
        Scalaire <- sum(Vect * contrasts) - value
      }
      Var <- t(contrasts) %*% Mat %*% contrasts
      Wald <- Scalaire/sqrt(Var)
      p_value <- 2 * (1 - pnorm(abs(Wald)))
      Results <- matrix(NA, nrow = 1, ncol = 4)
      colnames(Results) <- c("coef", "Se", 
                             "Wald Test", "p_value")
      if (is.null(name)) {
        if (is.null(value)) 
          value <- 0
        rownames(Results) <- paste(paste(names(Mod$best[pos]), 
                                         "*", contrasts, collapse = " + "), 
                                   "= ", value)
      }
      else {
        rownames(Results) <- name
      }
      Results[, 1] <- round(sum(Vect * contrasts), 5)
      Results[, 2] <- round(sqrt(Var), 5)
      Results[, 3] <- round(Wald, 5)
      Results[, 4] <- round(p_value, 5)
    }
    return(Results)
  }
}

