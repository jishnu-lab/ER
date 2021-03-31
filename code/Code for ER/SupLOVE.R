########################################################################
########                Algorithm ingredients                ###########
########################################################################
source("Tuning.R")
source("EstNonpure.R")
source("EstPure.R")
source("Utilities.R")
source("EstOmega.R")
source("ER-inference.R")
source("Est_beta_dz.R")
library(MASS)


ER <- function(Y, X, delta, beta_est = "NULL", CI = F, pred = T, lbd = 0.1, 
               rep_CV = 50, diagonal = F, merge = F, equal_var = F,
               alpha_level = 0.05, support = NULL, correction = "Bonferroni",
               verbose = F) {
  n <- nrow(X);  p <- ncol(X)
  if (equal_var) 
    se_est <- rep(1, p)
  else 
    se_est <- apply(X, 2, sd)
  
  deltaGrids <- delta * sqrt(log(max(p, n)) / n)
  
  # if (length(deltaGrids) > 1) {
  #   selected_deltas <- replicate(rep_CV, CV_Delta(X, deltaGrids, diagonal, se_est, merge))
  #   cat(table(selected_deltas))
  #   optDelta <- median(selected_deltas)
  # } else
  #   optDelta <- deltaGrids
  
  optDelta <- ifelse(length(deltaGrids) > 1,
                     median(replicate(rep_CV, CV_Delta(X, deltaGrids, diagonal, se_est, merge))),
                     deltaGrids)
  if (verbose)
    cat("Selecting the delta =", delta[min(which(deltaGrids >= optDelta))], "\n")
  
  Sigma <- crossprod(X) / n
  resultAI <- EstAI(Sigma, optDelta, se_est, merge)
  
  pure_numb <- sapply(resultAI$pureSignInd, FUN = function(x) {length(c(x$pos, x$neg))})
  if (sum(pure_numb == 1) > 0) {
    cat("Changing ``merge'' to ``union'' and reselect delta ... \n")
    optDelta <- ifelse(length(deltaGrids) > 1,
                       median(replicate(rep_CV, CV_Delta(X, deltaGrids, diagonal, se_est, merge = F))),
                       deltaGrids)
    resultAI <- EstAI(Sigma, optDelta, se_est, merge = F)
  }
  
  A_hat <- resultAI$AI
  I_hat <- resultAI$pureVec
  I_hat_ind <- resultAI$pureSignInd
  
  if (is.null(I_hat)) {
    cat("Algorithm fails due to non-existence of pure variable.\n")
    stop()
  }
  
  C_hat <- EstC(Sigma, A_hat, diagonal)
  Gamma_hat <- rep(0, p)
  Gamma_hat[I_hat] <- diag(Sigma[I_hat, I_hat]) - diag(A_hat[I_hat,] %*% C_hat %*% t(A_hat[I_hat,]))
  Gamma_hat[Gamma_hat < 0] <- 1e-2
  
  if (pred) {
    pred_result <- Prediction(Y, X, A_hat, Gamma_hat, I_hat)
    
    # the matrix to predict Z
    Theta_hat <- pred_result$Theta
    Q <- try(Theta_hat %*% solve(crossprod(X %*% Theta_hat) / n, crossprod(Theta_hat)), silent = T)
    if (class(Q)[1] == "try-error")
      Q <- Theta_hat %*% ginv(crossprod(X %*% Theta_hat) / n) %*% crossprod(Theta_hat)
  } else
    pred_result <- Q <- NULL
  

  
  if (beta_est == "NULL" || beta_est == "Dantzig") {
    beta_hat <- beta_CIs <- beta_var <- NULL
    if (beta_est == "Dantzig") {
      beta_hat <- Est_beta_dz(Y, X, A_hat, C_hat, I_hat, optDelta, 0.5, 0.5)
    }
  } else {
    if (CI) {
      if (length(resultAI$pureVec) != nrow(Sigma)) {
        Y_hat <- EstY(Sigma, A_hat, resultAI$pureVec)
        if (lbd > 0) {
          AI_hat <- abs(A_hat[I_hat, ])
          sigma_bar_sup <- max(solve(crossprod(AI_hat), t(AI_hat)) %*% se_est[I_hat])
          AJ <- EstAJDant(C_hat, Y_hat, lbd * optDelta * sigma_bar_sup, sigma_bar_sup + se_est[-I_hat])
        } else 
          AJ <- t(solve(C_hat, Y_hat))
        A_hat[-resultAI$pureVec, ] <- AJ  
      }
      
      Gamma_hat[-I_hat] <- diag(Sigma[-I_hat, -I_hat]) - diag(A_hat[-I_hat,] %*% C_hat %*% t(A_hat[-I_hat,]))
      Gamma_hat[Gamma_hat < 0] <- 1e2
    }
    res_beta <- Est_beta(Y, X, A_hat, C_hat, Gamma_hat, I_hat, I_hat_ind, CI = CI,
                         alpha_level = alpha_level, correction = correction)
    beta_hat <- res_beta$beta
    beta_CIs <- res_beta$CIs
    beta_var <- res_beta$beta_var
  }
  return(list(K = ncol(A_hat), A = A_hat, C = C_hat, I = I_hat, I_ind = I_hat_ind, Gamma = Gamma_hat,
              beta = beta_hat, beta_CIs = beta_CIs, pred = pred_result,
              optDelta = delta[min(which(deltaGrids >= optDelta))],
              Q = Q, beta_var = beta_var))
}





Prediction <- function(Y, X, A_hat, Gamma_hat, I_hat, option = "Theta") {
  K <- ncol(A_hat)
  p <- ncol(X); n <- nrow(X)
  
  R <- matrix(0, nrow = K, ncol = p)
  Sigma <- crossprod(X) / n
  
  BI <- solve(t(A_hat[I_hat,]) %*% A_hat[I_hat,], t(A_hat[I_hat,]))
  R[, I_hat] = BI %*% (Sigma[I_hat, I_hat] - diag(Gamma_hat[I_hat]))
  R[, -I_hat] = BI %*% Sigma[I_hat, -I_hat]
  
  
  Q <- X %*% t(R)
  theta_hat <- try(t(R) %*% solve(crossprod(Q), crossprod(Q, Y)), silent = T)
  if (class(theta_hat)[1] == "try-error")
    theta_hat <- t(R) %*% ginv(crossprod(Q)) %*% crossprod(Q, Y)
  
  pred_val <- X %*% theta_hat
  return(list(theta = theta_hat, fit = pred_val, Theta = t(R)))
}







################################################################################
###           Functions for re-fitting by using the support of beta          ###
################################################################################


Pred_Z_BLP <- function(X, A_hat, C_hat, est_Gamma, S_beta) {
  est_Gamma_inv <- diag(est_Gamma ** (-1))
  G_hat <- crossprod(A_hat, est_Gamma_inv) %*% A_hat + solve(C_hat)
  Z_hat <- X %*% est_Gamma_inv %*% A_hat %*% ginv(G_hat)
  Z_hat[,S_beta,drop = F]
}

Pred_Z_BLP_avg <- function(X, A_hat, C_hat, S_beta) {
  G_hat <- crossprod(A_hat) + solve(C_hat)
  Z_hat <- X %*% A_hat %*% ginv(G_hat)
  Z_hat[,S_beta,drop = F]
}

# Pred_avg <- function(Y, X, A_hat, S, S_beta) {
#   B_hat <- t(solve(crossprod(A_hat[S,]), t(A_hat[S, ])))[, S_beta]
#   Z_hat <- X[ ,S] %*% B_hat
#   eta_hat <- try(solve(crossprod(Z_hat), crossprod(Z_hat, Y)), silent = T)
#   if (class(eta_hat) == "try-error")
#     eta_hat <- ginv(crossprod(Z_hat)) %*% crossprod(Z_hat, Y)
#   theta_hat <- B_hat %*% eta_hat
#   return(list(pred = X[ ,S] %*% theta_hat, theta = theta_hat))
# }








