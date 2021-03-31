##### Parameter tuning and Cross validation ######
source("Utilities.R")

CV_Delta <- function(X, deltaGrids, diagonal, se_est, merge) {
  # Cross validation for choosing delta. For each delta from the given grids, 
  # first split the data into two datasets. Obtain I, AI and C. Then calculate 
  # AICAI^T and choose delta which minimizes the criterion 
  #                   A_I C A_I' - Sigma(dataset 2)
  #
  # Args: 
  #   X: n by p matrix.
  #   deltaGrids: vector of numerical constants.
  #   diagonal: TRUE or FALSE for the diagonal structure of C.
  # 
  # Returns: 
  #   the selected optimal delta
  n <- nrow(X); p <- ncol(X)
  sampInd <- sample(n, floor(n / 2))
  X1 <- X[sampInd, ]
  X2 <- X[-sampInd, ]
  Sigma1 <- crossprod(X1) / nrow(X1);  
  diag(Sigma1) <- 0
  Sigma2 <- crossprod(X2) / nrow(X2)
  
  result_Ms <- FindRowMax(abs(Sigma1))
  Ms <- result_Ms$M
  arg_Ms <- result_Ms$arg_M
  
  loss <- c()
  for (i in 1:length(deltaGrids)) {
    resultFitted <- CalFittedSigma(Sigma1, deltaGrids[i], Ms, arg_Ms, se_est, 
                                   diagonal, merge)
    fittedValue <- resultFitted$fitted   
    estPureVec <- resultFitted$pureVec
    if (is.null(dim(fittedValue)) && fittedValue == -1)
      loss[i] <- Inf
    else {
      denom <- length(estPureVec) * (length(estPureVec) - 1)
      # loss[i] <- 2 * offSum(Sigma2[estPureVec, estPureVec], fittedValue, 1) / denom
      loss[i] <- 2 * offSum(Sigma2[estPureVec, estPureVec], fittedValue, se_est[estPureVec]) / denom
    }
  } 
  # cat(loss)
  return(deltaGrids[which.min(loss)])
}


CalFittedSigma <- function(Sigma, delta, Ms, arg_Ms, se_est, diagonal, merge) {
  # Calculate the fitted value of A_ICA_I^T for given Sigma and delta.
  # 
  # Args: 
  #   Sigma: p by p covariance matrix. 
  #   delta: given parameter. 
  #   Ms: the calculated maximal values of Sigma per row.
  #   diagonal: TRUE or FALSE.
  #
  # Returns: 
  #   a list including:
  #     pureVec: vector of the indices of estimated pure variables.
  #     fitted: fitted value of A_ICA_I^T 
  resultPureNode <- FindPureNode(abs(Sigma), delta, Ms, arg_Ms, se_est, merge)
  
  estPureIndices <- resultPureNode$pureInd
  # lapply(estPureIndices, function(x) cat(x, "\n"))
  
  if (singleton(estPureIndices))
    return(list(pureVec = NULL, fitted = -1))
  
  estSignPureIndices <- FindSignPureNode(estPureIndices, Sigma)
  AI <- RecoverAI(estSignPureIndices, length(se_est))
  C <- EstC(Sigma, AI, diagonal)
  
  if (length(estPureIndices) == 1)
    fitted <- -1
  else {
    subAI <- AI[resultPureNode$pureVec, ]
    fitted <- subAI %*% C %*% t(subAI)
  }
  return(list(pureVec = resultPureNode$pureVec, fitted = fitted))
}


CV_lbd <- function(X, lbdGrids, AI, pureVec, diagonal) {
  # Use Cross-validation to select lambda for estimating Omega. Split the data
  # into two parts. Estimating C on two datasets. Then, for each lambda, calcu-
  # late Omega on the first dataset and calculate the loss on the second dataset.
  # Find the lambda which gives the smallest loss of 
  #                 <C,Omega> - log(det(Omega))
  #
  # Args: 
  #   X: n by p matrix.
  #   lbdGrids: vector of numerical constants.
  #   AI: p by K matrix.
  #   pureVec: vector of indices of pure variables.
  #   diagonal: TRUE or FALSE for the diagonal structure of C.
  # 
  # Returns: 
  #   the selected optimal lambda
  sampInd <- sample(nrow(X), floor(nrow(X) / 2))
  X1 <- X[sampInd, ]
  X2 <- X[-sampInd, ]
  Sigma1 <- crossprod(X1) / nrow(X1)
  Sigma2 <- crossprod(X2) / nrow(X2)
  C1 <- EstC(Sigma1, AI, diagonal)
  C2 <- EstC(Sigma2, AI, diagonal)
  loss <- c()
  for (i in 1:length(lbdGrids)) {
    Omega <- estOmega(lbdGrids[i], C1)
    det_Omega <- det(Omega)
    loss[i] <- ifelse(det_Omega <= 0, Inf, sum(Omega * C2) - log(det_Omega))
  }
  return(lbdGrids[which.min(loss)])
}
