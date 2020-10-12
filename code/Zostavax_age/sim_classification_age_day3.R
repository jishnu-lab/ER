################################################################################
##########                                                            ##########
##########                      ER  application                      ###########
##########           prediction / classification  performance        ###########
##########                                                            ##########
################################################################################
rm(list=ls())
setwd("~/Documents/Mike/Projects/Application of Essential Regression/code/Code for ER")
source("SupLOVE.R")
setwd("~/Documents/Mike/Projects/Application of Essential Regression/code")
source("K-CV.R")
library(ROCR)
library(readr)

data <- read_delim("../dataset/age_imputed_5NN.txt", " ",
                   escape_double = FALSE, trim_ws = TRUE)
Y <- data$age
X <- data[,-(1:2)]



CV_binary <- function(Y, X, methods, delta_grid = NULL, ER_sparse = F) {
  n <- length(Y)
  pred_val <- matrix(NA, n, length(methods))
  colnames(pred_val) <- methods
  
  Y_label <- rep(0, n)
  Y_label[Y >= 50] = 1
  
  Y_perm <- Y[sample(1:n)]
  
  for (i in 1:n) {
    trainX <- X[-i,]
    trainY <- Y[-i]
    validX <- X[i,]
    
    for (j in 1:length(methods)) {
      method_j <- methods[j]
      if (method_j == "Perm-ER" | method_j == "Perm-Lasso") 
        trainY <- Y_perm[-i]
      
      pred_val[i, j] <- calPred(validX, trainX, trainY, method_j, 
                                delta_grid = delta_grid, ER_sparse = ER_sparse)
    }
  }
  return(pred_val)
}

Compute_metric <- function(pred_val, Y_label) {
  methods <- colnames(pred_val)
  error_list <- auc_vec <- c()
  for (l in 1:length(methods)) {
    pred_roc <- prediction(pred_val[,l], Y_label)
    perf <- performance(pred_roc, "tpr", "fpr")
    error_list[[l]] <- rbind("tpr" = perf@y.values[[1]], "fpr" = perf@x.values[[1]])
    perf_auc <- performance(pred_roc, "auc")
    auc_vec[l] <- as.numeric(perf_auc@y.values)
  }
  names(error_list) <- methods
  return(list(error = error_list, auc = auc_vec))
}

Plot_ROC <- function(mean_errors, methods) {
  par(mar = c(4,4,1,1))
  n_method <- length(mean_errors)
  sel_names <- c()
  for (i in 1:n_method) {
    if (! names(mean_errors)[i] %in% methods)
      next
    sel_names <- c(sel_names, names(mean_errors)[i])
    count <- length(sel_names)
    Y_i <- mean_errors[[i]]['tpr',]
    X_i <- mean_errors[[i]]['fpr',]
    if (count == 1)
      plot(X_i, Y_i, type = "l", col = count + 1, lty = count, main = "", xlab = "False positive rate",
           ylab = "True positive rate", ylim = c(0, 1), xlim = c(0, 1), lwd = 2)
    else
      points(X_i, Y_i, type = "l", col = count + 1, lty = count, lwd = 2)
  }
  abline(a = 0, b = 1, col = 1, lty = 2, lwd = 1)
  legend("topleft", legend = sel_names, col = 2:(count+1), lty = 1:count, lwd = rep(2, count), bty = "n")
}


################################################################################
######                    Classification performance                      ######
################################################################################
X <- scale(X, T, T)

methods = c("ER", "PFR", "PLS", "Lasso", "Perm-Lasso")


set.seed(20200316)

# pred_val <- CV_binary(Y, X, methods, seq(0.25, 0.5, 0.02), T)

# write.table(pred_val, file = "../output/age_dataset/pred_val_binary.txt", row.names = F, col.names = T)


pred_val <- read.table("../output/age_dataset/pred_val_binary.txt", header = T)
colnames(pred_val) <- methods

true_label <- rep(0, length(Y))
true_label[Y >= 50] = 1

result <- Compute_metric(pred_val, true_label)

# auc  
round(result$auc, 2)    
# 0.66 0.34 0.50 0.75 0.53


Plot_ROC(result$error, methods)
Plot_ROC(result$error, setdiff(methods, "Perm-Lasso"))














