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
source("Other_algorithms.R")
source("Helper.R")

library(readr)
data <- read_delim("../dataset/age_imputed_5NN.txt", " ",
                   escape_double = FALSE, trim_ws = TRUE)

# data <- read_delim("../dataset/age_D7_imputed_5NN.txt", " ", 
#                    escape_double = FALSE, trim_ws = TRUE)

Y <- data$age
X <- data[,-(1:2)]

feature_names <- names(X)

################################################################################
#######                       ER on the whole dataset                    #######
################################################################################

X <- scale(X, T, T)
Y <- Y - mean(Y)
# 
# set.seed(20200316)
# delta_grid <- seq(0.25, 0.75, 0.01)
#
# res <- ER(Y, X, delta_grid, pred = F, beta_est = "NULL", rep_CV = 200, verbose = T)
#
# res$K
## selected delta = 0.27

res <- ER(Y, X, 0.27, pred = F, beta_est = "LS", rep_CV = 50, verbose = T,
          merge = F, CI = T)
res$K
round(res$beta_CIs, 3)
p_vals <- pnorm(abs(res$beta / sqrt(res$beta_var / length(Y))), lower.tail = F)

which(p_vals <= 0.05)





Z_ind <- which(p_vals <= 0.05)
A_hat <- res$A
feature_ind <- which(rowSums(abs(A_hat[,Z_ind, drop = F])) != 0)

set.seed(20200415)

cvfit = cv.glmnet(X[,feature_ind,drop = F], Y, standardize = F, intercept = F)
sign_coef <- predict(cvfit, s = "lambda.min", type = "coefficient")[-1]
coef_ind <- feature_ind[which(sign_coef != 0)]


Lasso_cvfit = cv.glmnet(X, Y, standardize = F, intercept = F)
Lasso_sign_coef <- predict(Lasso_cvfit, s = "lambda.min", type = "coefficient")[-1]
Lasso_coef_ind <- which(Lasso_sign_coef != 0)

cat("ER-Lasso selects", length(coef_ind), "features while Lasso selects",
    length(Lasso_coef_ind), "features with", length(intersect(Lasso_coef_ind, coef_ind)),
    "features in common.")


file_name <- "../output/age_dataset/day_3_var_selection.txt"
write(feature_names[coef_ind], file = file_name)
write("\n", file = file_name, append = T)
write(feature_names[Lasso_coef_ind], file = file_name, append = T)




################################################################################
######                       Prediction performance                       ######
################################################################################
X <- scale(X, T, T)

set.seed(20200309)

methods <- c("ER", "ER-Lasso-Dantzig", "PFR", "PLS", "Lasso", "Perm-Lasso")
# methods <- c("PFR", "PLS", "Lasso", "Perm-Lasso")

result <- replicate(50, KfoldCV(5, Y = Y, X = X, methods = methods,
                                pt_option = "perm_X_Y", metric = "mse",
                                delta_grid = seq(0.2, 0.5, 0.02), ER_sparse = F))

CV_result <- t(sapply(result['MSE',], function(v) v))
fitted_array <- sapply(result['fitted',], function(v) v, simplify = "array")

fitted_mat <- apply(fitted_array, 1:2, mean)
CV_result <- CV_result / mean(Y ** 2)

# write.table(CV_result, file = "../output/age_dataset/D3_CV_result_04_11_mse.txt", row.names = F, col.names = T)
# write.table(fitted_mat, file = "../output/age_dataset/D3_CV_result_04_11_fit.txt", row.names = F, col.names = T)

CV_result <- read.table("../output/age_dataset/D3_CV_result_04_11_mse.txt", header = T)
fitted_mat <- read.table("../output/age_dataset/D3_CV_result_04_11_fit.txt", header = T)

# wilcox.test(CV_result$ER, CV_result$Perm.Lasso, alternative = "less")  ##  p-value < 2.2e-16



Violin_plot(methods, CV_result, MSE_flag = T, ylab_name = "1 - standardized mean square error",
            label_names = methods[c(1,2,5,3,4,6)], ylims = c(0.6, max(1 - CV_result)))

Fitted_plot(Y, fitted_mat, methods, order = c(1,2,5,3,4,6))






