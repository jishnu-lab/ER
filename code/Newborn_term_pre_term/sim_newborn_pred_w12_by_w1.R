
rm(list=ls())
setwd("~/Documents/Mike/Projects/Application of Essential Regression/code/Code for ER")
source("SupLOVE.R")
setwd("~/Documents/Mike/Projects/Application of Essential Regression/code")
source("K-CV.R")
source("Other_algorithms.R")
source("Helper.R")
library(readr)
library(ROCR)


child_data <- read.csv("../dataset/Newborndata/data_child.csv")
# View(child_data)

### Week 12    dims: 46 282
subdata_w12 <- subset(child_data, child_data$Time.point == "Week 12", c("Subject", "Neutrophils"))
# head(subdata_w12)

### Week 1     dims: 56  282
subdata_w1 <- subset(child_data, child_data$Time.point == "Week 1")

combined_data <- merge(subdata_w12, subdata_w1, by = "Subject")
na_row <- which(rowSums(is.na(combined_data)) > 0)

combined_data <- combined_data[-na_row, ]
# View(combined_data)

X <- combined_data[,-c(1:8, 283)]   # the 282th column is the NA flag 
# dims: 43 274
Y <- combined_data$Neutrophils.x

feature_names <- names(combined_data)[-c(1:8, 283)]


################################################################################
#######                       ER on the whole dataset                    #######
################################################################################

X <- scale(X, T, T)
Y <- Y - mean(Y)

# set.seed(20200316)
# delta_grid <- seq(0.05, 0.6, 0.01)
# 
# res <- ER(Y, X, delta_grid, pred = F, beta_est = "NULL", rep_CV = 100, verbose = T, merge = F)
# 
# res$K
# selected delta = 0.23
# 
lbd = 0.5
delta = 0.23

res <- ER(Y, X, delta = delta, pred = T, beta_est = "LS", rep_CV = 50, verbose = F,
          merge = F, CI = T, correction = "Bonferroni", lbd = lbd)
res$K
round(res$beta_CIs, 3)
p_vals <- 2 * pnorm(abs(res$beta) / sqrt(res$beta_var / length(Y)), lower.tail = F)
# Z_ind <- which(p_vals <= 0.1)


# clusters <- recoverGroup(res$A)
# list_clusters_pure <- lapply(Z_ind, FUN = function(x, data_list, feature_name) {
#   feature_name[unlist(data_list[[x]])]
# }, data_list = res$I_ind, feature_name = colnames(X))
#
#
# list_clusters <- lapply(Z_ind, FUN = function(x, data_list, feature_name) {
#   feature_name[unlist(data_list[[x]])]
# }, data_list = clusters, feature_name = colnames(X))
#
# sapply(list_clusters_pure, FUN = function(x) {length(x)})
# sapply(list_clusters, FUN = function(x) {length(x)})
#
# output cluster results
# save_LOVE_result(list_clusters_pure, delta, lbd, only.pure = T, "newborn_dataset/pred_week12")
# save_LOVE_result(list_clusters_pure, delta, lbd, only.pure = F, "newborn")



Z_ind <- which(p_vals <= 0.05)   #  group 4, 5, 9
A_hat <- res$A
feature_ind <- which(rowSums(abs(A_hat[,Z_ind, drop = F])) != 0)

set.seed(20200316)

cvfit = cv.glmnet(X[,feature_ind,drop = F], Y, standardize = F, intercept = F)
sign_coef <- predict(cvfit, s = "lambda.min", type = "coefficient")[-1]
coef_ind <- feature_ind[which(sign_coef != 0)]


Lasso_cvfit = cv.glmnet(X, Y, standardize = F, intercept = F)
Lasso_sign_coef <- predict(Lasso_cvfit, s = "lambda.min", type = "coefficient")[-1]
Lasso_coef_ind <- which(Lasso_sign_coef != 0)

cat("ER-Lasso selects", length(coef_ind), "features while Lasso selects",
    length(Lasso_coef_ind), "features with", length(intersect(Lasso_coef_ind, coef_ind)),
    "features in common.")


file_name <- "../output/newborn_dataset/var_selection.txt"
write(feature_names[coef_ind], file = file_name)
write("\n", file = file_name, append = T)
write(feature_names[Lasso_coef_ind], file = file_name, append = T)



################################################################################
######                      Prediction   performance                      ######
################################################################################
X <- scale(X, T, T)

methods <- c("ER", "PFR", "PLS", "Lasso", "Perm-Lasso", "ER-Lasso-LS")

set.seed(20200309)

delta_grid <- seq(0.05, 0.6, 0.02)

result <- replicate(50, KfoldCV(10, Y = Y, X = X, methods = methods, 
                                pt_option = "perm_Y", metric = "adv", 
                                delta_grid = delta_grid, ER_sparse = F))


CV_result <- t(sapply(result['MSE',], function(v) v))
fitted_array <- sapply(result['fitted',], function(v) v, simplify = "array")

fitted_mat <- apply(fitted_array, 1:2, mean)

summary(CV_result, digits = 4)

# write.table(CV_result, file = "../output/newborn_dataset/pred_w12_by_w1_mse.txt", row.names = F, col.names = T)
# write.table(fitted_mat, file = "../output/newborn_dataset/pred_w12_by_w1_fit.txt", row.names = F, col.names = T)


CV_result <- read.table("../output/newborn_dataset/pred_w12_by_w1_mse.txt", header = T)
fitted_mat <- read.table("../output/newborn_dataset/pred_w12_by_w1_fit.txt", header = T)

# wilcox.test(CV_result$ER, CV_result$Perm.Lasso, alternative = "less")  ##  p-value < 2.2e-16

methods[6] <- "Composite"

Violin_plot(methods, CV_result, MSE_flag = T, ylab_name = "1 - absolute deviance error",
            label_names = methods[c(1,6,4,2,3,5)])

Fitted_plot(Y, fitted_mat, methods, order = c(1,6,4,2,3,5))




