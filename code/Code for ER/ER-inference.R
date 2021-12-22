########################################################################
########                Algorithm ingredients                ###########
########################################################################


Est_beta <- function(Y, X, A_hat, C_hat, Gamma_hat, I_hat, I_ind_hat, CI = T, 
                     alpha_level = 0.05, correction = NULL, support = NULL) {
  n <- nrow(X); p <- ncol(X)
  R <- matrix(0, nrow = p, ncol = ncol(A_hat))
  Sigma <- crossprod(X) / n
  
  BI <- t(solve(crossprod(A_hat[I_hat,]), t(A_hat[I_hat,])))
  R[I_hat, ] = (Sigma[I_hat, I_hat] - diag(Gamma_hat[I_hat])) %*% BI
  R[-I_hat, ] =  Sigma[-I_hat, I_hat] %*% BI
  
  if (is.null(support))
    beta_est <- solve(crossprod(R), t(R) %*% crossprod(X, Y) / n)
  else {
    beta_est <- rep(0, ncol(A_hat))
    beta_est[support] <- solve(crossprod(R[, support]), t(R[,support]) %*% crossprod(X, Y) / n) 
  }
  
  sigma2_hat <- Est_sigma2(Y, t(BI) %*% crossprod(X[,I_hat], Y) / n, beta_est, C_hat)
  Omega_hat <- solve(C_hat)
  
  
  if (CI) {  
    ### (1 - alpha / 2) CI of beta (for each individual component)
    
    beta_var <- Comp_ASV(sigma2_hat, BI, R, Gamma_hat, beta_est, Omega_hat, I_hat, I_ind_hat)
    beta_var[beta_var < 0] = 0
    
    if (!is.null(correction)) {
      if (correction == "Bonferroni") 
        alpha_level <- alpha_level / ncol(A_hat)
    }
    
    CIs_beta <- cbind(lower = beta_est - qnorm(1 - alpha_level / 2) * sqrt(beta_var) / sqrt(n),
                      upper = beta_est + qnorm(1 - alpha_level / 2) * sqrt(beta_var) / sqrt(n))
  } else
    CIs_beta <- beta_var <- NULL
  
  return(list(beta = beta_est, CIs = CIs_beta, beta_var = beta_var))
}




Est_sigma2 <- function(Y, h, beta_est, C_hat) {
  n <- length(Y)
  sigma2_hat <- crossprod(Y) / n - 2 * t(beta_est) %*% h + t(beta_est) %*% C_hat %*% beta_est
  ifelse(sigma2_hat < 0, 0, sigma2_hat)
}




Comp_ASV <- function(sigma2_hat, BI, R, Gamma_hat, beta_hat, Omega_hat, I_hat, I_ind_hat) {
  D_tau_bar <- t(BI) %*% diag(Gamma_hat[I_hat]) %*% BI
  V1 <- as.numeric(sigma2_hat + t(beta_hat) %*% D_tau_bar %*% beta_hat)
  Q <- t(solve(crossprod(R), t(R)))
  V2 <- Omega_hat + t(Q) %*% diag(Gamma_hat) %*% Q

  K_hat <- length(I_ind_hat)


  D <- c()
  ms <- c()
  for (a in 1:K_hat) {
    group_a <- unlist(I_ind_hat[[a]])
    m_a <- length(group_a)
    ms[a] <- m_a

    D1 <- (2 * m_a - 1) * (D_tau_bar[a,a] ^ 2) * (m_a ^ 2) - sum(Gamma_hat[group_a] ** 2)
    D[a] <- D1 * (beta_hat[a] ^ 2) / m_a / ((m_a - 1) ^ 2)
  }

  V3 <- t(Q[I_hat,]) %*% diag(rep(D, ms)) %*% Q[I_hat,]
  return(diag(V1 * V2 + V3))
}










