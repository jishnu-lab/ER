#################################################################################
#####                                                                       #####
##### This is the code for Step 2: regression to estimate the non-pure rows #####
#####                                                                       #####
#################################################################################
library(linprog)

EstC <- function(Sigma, AI, diagonal) {
  # Estimate C. If diagonal=True, estimate only diagonal elements.
  #
  # Args: 
  #   Sigma: p by p covariance matrix.
  #   AI: p by K matrix.
  #   diagonal: TRUE or FALSE.
  #
  # Returns: 
  #   K by K estimated C_hat 
  K <- ncol(AI)
  C <- diag(0, K, K)
  for (i in 1:K) {
    groupi <- which(AI[ ,i] != 0)
    sigmai <- as.matrix(abs(Sigma[groupi,groupi]))
    tmpEntry <- sum(sigmai) - sum(diag(sigmai))
    C[i,i] <- tmpEntry / (length(groupi) * (length(groupi) - 1))
    if (!diagonal && i < K) {
      for (j in (i+1):K) {
        groupj <- which(AI[ ,j]!=0)
        # adjust the sign for each row
        sigmaij <- AI[groupi,i] * as.matrix(Sigma[groupi, groupj]) 
        sigmaij <- t(AI[groupj, j] * t(sigmaij))
        C[i,j] <- C[j,i] <- sum(sigmaij) / (length(groupi) * length(groupj))
      }
    }
  }
  return(C)
}

EstY <- function(Sigma, AI, pureVec) {
  # Estimate the Sigma_{TJ}_hat
  # 
  # Args: 
  #   Sigma: p by p matrix.
  #   AI: p by K matrix.
  #   pureVec: the vector of pure node indices.
  #
  # Returns: 
  #   Sigma_{TJ}_hat: K by |J| matrix.
  SigmaS <- AdjustSign(Sigma, AI)
  SigmaJ <- matrix(SigmaS[ ,-pureVec], nrow = nrow(Sigma))
  SigmaTJ <- matrix(0, ncol(AI), nrow(AI) - length(pureVec))
  for (i in 1:ncol(AI)) {
    groupi <- which(AI[ ,i] != 0)
    SigmaiJ <- as.matrix(SigmaJ[groupi, ])
    SigmaTJ[i, ] <- apply(SigmaiJ, 2, mean) # Average columns along the rows.
  }
  return(SigmaTJ)
}

AdjustSign <- function(Sigma, AI) {
  # Sign operation on matrix {@code Sigma} according to the sign of {@code AI}
  #
  # Args: 
  #   Sigma: p by p matrix; AI: p by K matrix
  #
  # Returns: 
  #   p by p matrix
  SigmaS <- matrix(0, nrow(AI), nrow(AI))
  for (i in 1:nrow(AI)) {
    index <- which(AI[i, ] != 0)
    if (length(index) != 0)
      SigmaS[i, ] = sign(AI[i,index]) * Sigma[i, ]
  }
  return(SigmaS)
}

EstAJInv <- function(Omega, Y, lbd) {
  # This function estimates the |J| by K matrix A_J by using soft thresholding.
  #
  # Args: 
  #   Omega: K by K estimated C^{-1}.
  #   Y: K by |J| reponse matrix.
  #   lbd: the chosen constant of the RHS constraint for soft-thresholding.
  #
  # Returns: 
  #   |J| by K matrix A_J.
  AJ <- matrix(0, ncol(Y), nrow(Y))
  for (i in 1:ncol(Y)) {
    Atilde <- Omega %*% as.matrix(Y[ ,i])
    AJ[i, ] <- LP(Atilde, lbd)
    if (sum(abs(AJ[i, ])) > 1)
      AJ[i,] <- AJ[i,] / sum(abs(AJ[i, ]))
  }
  return(AJ)
}


LP <- function(y, lbd) {
  # Use LP to solve problem:
  #   beta^+ - beta^- \leq lbd + y
  #   - beta^+ + beta^- \leq lbd - y
  #   beta^+ \ge 0; beta^- \ge 0
  #
  # Args: 
  #   y: K by 1 vector. 
  #   lbd: positive contant.
  #
  # Returns:
  #   K by 1 vector
  K <- length(y)
  cvec <- rep(1, 2 * K)
  bvec <- c(lbd + y, lbd - y, rep(0, 2 * K))
  C <- matrix(0, K, 2 * K) 
  for (i in 1:K) {
    indices <- c(i,i + K)
    C[i,indices] = c(1,-1)
  }
  Amat <- rbind(C, -C, diag(-1, nrow = 2 * K))
  LPsol <- solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  beta <- LPsol[1:K] - LPsol[(K + 1):(2 * K)]
  return(beta)
}

EstAJDant <- function(C_hat, Y_hat, lbd, se_est_J) {
  AJ <- matrix(0, ncol(Y_hat), nrow(Y_hat))
  for (i in 1:ncol(Y_hat)) {
    AJ[i, ] <- Dantzig(C_hat, Y_hat[,i], lbd * se_est_J[i])
    if (sum(abs(AJ[i, ])) > 1)
      AJ[i,] <- AJ[i,] / sum(abs(AJ[i, ]))
    # cat("Finishing estimating the", i, "th row...\n")
  }
  return(AJ)
}

Dantzig <- function(C_hat, y, lbd) {
  K <- length(y)
  cvec <- rep(1, 2 * K)
  bvec <- c(lbd + y, lbd - y, rep(0, 2 * K))
  new_C_hat <- matrix(0, K, 2 * K)
  for (i in 1:K) {
    new_C_hat[i, ] <- c(C_hat[i,], -C_hat[i,])
  }
  Amat <- rbind(new_C_hat, -new_C_hat, diag(-1, nrow = 2 * K))
  LPsol <- solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  beta <- LPsol[1:K] - LPsol[(K + 1):(2 * K)]
  return(beta)
}





