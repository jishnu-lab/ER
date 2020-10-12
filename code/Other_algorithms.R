library(glmnet)



LASSO <- function(X, Y, new_X) {
  cvfit = cv.glmnet(X, Y, alpha = 1, type.measure = "mse", nfolds = 10, intercept = F)
  predict(cvfit, newx = new_X, s = "lambda.min")
}


Est_K <- function(X, K_max = 30) {
  n <- nrow(X); p <- ncol(X)
  U <- svd(X, nv = 0)$u
  penalty <- (n + p) / n / p * log(p * n / (p + n))
  PCs <- t(U) %*% X
  loss <- vector("numeric", K_max)
  for (k in 1:K_max) {
    loss[k] <- log(norm(X - as.matrix(U[,1:k]) %*% PCs[1:k,], "F") ** 2 / n / p) + k * penalty
  }
  which.min(loss)
}


Est_K_ratio <- function(X, K_max = 30) {
  eigens <- (svd(X, nu = 0, nv = 0)$d[1:K_max]) ** 2
  L <- length(eigens)
  ratios <- eigens[1:(L-1)] / eigens[2:L]
  which.max(ratios)
}


PCR <- function(X, Y, K = NULL, option = "additive") {
  n <- nrow(X); p <- ncol(X)
  K_max <- min(n, p) %/% 2
  if (is.null(K)) {
    if (option == "additive")
      K_hat <- Est_K(X, K_max)
    else 
      K_hat <- Est_K_ratio(X, K_max)
  } else
    K_hat <- K
  
  K_hat <- ifelse(K_hat == 0, 1, K_hat)  # if the selected K is 0, set it to 1 instead
  
  U <- svd(X, nu = 0, nv = K_hat)$v
  Z_hat <- X %*% U / p
  beta_hat <- solve(crossprod(Z_hat), crossprod(Z_hat, Y))
  return(list(theta = U %*% beta_hat / p, fit = Z_hat %*% beta_hat, K = K_hat,
              Z = Z_hat, A = U / p))
}




PLS <- function(trainX, trainY, validX) {

  whole_data <- data.frame(Y = c(trainY, rep(NA, nrow(validX))), X = rbind(trainX, validX), row.names = NULL)
  training_index <- 1:nrow(trainX)
  fit_pls <- plsr(Y ~ ., data = whole_data[training_index, ], validation = "CV", segments = 5)
  n_comp <- selectNcomp(fit_pls, method = "randomization")
  n_comp <- ifelse(n_comp == 0, 1, n_comp)
  res_pred <- predict(fit_pls, comps = n_comp, newdata = whole_data[-training_index, ],
                      type = "response")
  return(res_pred)

}








