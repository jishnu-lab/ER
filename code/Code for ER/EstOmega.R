################################################################################
######                                                                   #######
######   These are the code which convert the original problem into LP   #######
######                                                                   #######
################################################################################
library(linprog)

solve_row <- function(col_ind, C, lbd) {
  K <- nrow(C)
  cvec <- c(1, rep(0, 2*K))
  Amat <- -cvec
  Amat <- rbind(Amat, c(-1, rep(1, 2*K)))
  tmp_constr <- C %x% t(c(1,-1))
  Amat <- rbind(Amat, cbind(-1 * lbd, rbind(tmp_constr, -tmp_constr)))
  tmp_vec <- rep(0, K)
  tmp_vec[col_ind] <- 1
  bvec <- c(0, 0, tmp_vec, -tmp_vec)
  
  lpResult <- solveLP(cvec, bvec, Amat, lpSolve = T)$solution
  while (length(lpResult) == 0) {
    cat("The penalty lambda =", lbd, "is too small and increased by 0.01...\n")
    lbd <- lbd + 0.01
    Amat[-(1:2), 1] <- lbd
    lpResult <- solveLP(cvec, bvec, Amat, lpSolve = T)$solution[-1]
  }
  ind <- seq(2, 2*K, 2)
  return(lpResult[ind] - lpResult[ind + 1])
}

estOmega <- function(lbd, C) {
  # For a given lbd and C, solve the C^{-1}. 
  # Require: C should be symmetric and square. 
  K <- nrow(C)
  omega <- matrix(0, K, K)
  for (i in 1:K) {
    omega[,i] <- solve_row(i, C, lbd)
  }
  return(omega)
}

