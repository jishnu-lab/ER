
rm(list=ls())
setwd("~/Documents/Mike/Projects/Application of Essential Regression/code/Code for ER")
source("SupLOVE.R")
setwd("~/Documents/Mike/Projects/Application of Essential Regression/code")
source("K-CV.R")
library(ggplot2)
library(readr)
library(ROCR)

child_data <- read.csv("../dataset/Newborndata/data_child.csv")
# View(child_data)

### Week 1    dims: 56  282
subdata_w1 <- subset(child_data, child_data$Time.point == "Week 1", c("Subject", "Neutrophils"))
head(subdata_w1)

### Week 12   dims: 46 282
subdata_w12 <- subset(child_data, child_data$Time.point == "Week 12")

combined_data <- merge(subdata_w1, subdata_w12, by = "Subject")
# View(combined_data)

raw_X <- combined_data[,-c(1:8, 283)]   # the 282th column is the NA flag

### Impute data via KNN
Impute_K_NN <- function(K, na_col_ind, data_matrix) {
  # na_row_ind <- 7:32
  new_data <- data_matrix
  for (sub_ind in na_col_ind) {
    na_row_ind <- which(is.na(data_matrix[, sub_ind]))
    sub_data <- data_matrix[-na_row_ind, sub_ind]
    K_neighbor <- CompKNN(K, sub_data, na_col_ind, na_row_ind, data_matrix)
    new_data[na_row_ind, sub_ind] <- apply(data_matrix[na_row_ind, K_neighbor], 1, mean)
  }
  return(new_data)
}
CompKNN <- function(K, sub_data, na_col_ind, na_row_ind, data_matrix) {
  candidates <- setdiff(1:ncol(data_matrix), na_col_ind)
  distance <- rep(Inf, length(candidates))
  distance[candidates] <- sapply(candidates, FUN = function(x) {
    mean((data_matrix[-na_row_ind, x] - sub_data) ** 2)
  })
  neighbor_ind <- order(distance)[1:K]
  return(neighbor_ind)
}

data_matrix <- t(raw_X)
na_col_ind <- which(is.na(colSums(data_matrix)))
imputed_data <- Impute_K_NN(5, na_col_ind, data_matrix)

X <- t(imputed_data)    ###  dims:  45 by 274 
Y <- combined_data$Neutrophils.x


################################################################################
#######                       ER on the whole dataset                    #######
################################################################################

# X <- scale(X, T, T)
# Y <- Y - mean(Y)
# 
# set.seed(20200316)
# delta_grid <- seq(0.05, 0.4, 0.01)
# 
# res <- ER(Y, X, delta_grid, pred = F, beta_est = "NULL", rep_CV = 100, verbose = T, merge = F)
# 
# res$K
# selected delta = 0.13

# lbd = 0.5
# delta = 0.13
# 
# res <- ER(Y, X, delta = delta, pred = F, beta_est = "LS", rep_CV = 50, verbose = T,
#           merge = F, CI = T, correction = "Bonferroni", lbd = lbd)
# res$K
# round(res$beta_CIs, 3)
# p_vals <- 2 * pnorm(abs(res$beta) / sqrt(res$beta_var / length(Y)), lower.tail = F)
# 
# clusters <- recoverGroup(res$A)
# 
# sapply(clusters, function(x) as.vector(sort(unlist(x))))
# sapply(res$I_ind, function(x) as.vector(unlist(x)))


# output cluster results
# save_LOVE_result <- function(clusters, delta, lbd, only.pure = F) {
#   K <- length(clusters)
#   if (only.pure)
#     savePath = paste(format(Sys.time(),"../output/Cibersort_dataset/love_%m.%d_nGroup"),K,"delta=",delta,"lbd=",lbd,"_pure",sep='')
#   else
#     savePath = paste(format(Sys.time(),"../output/Cibersort_dataset/love_%m.%d_nGroup"),K,"delta=",delta,"lbd=",lbd, " ", sep='')
#   dir.create(file.path(getwd(),savePath))
# 
#   foreach (k = 1:K) %do% {
#     groupk = sort(unlist(clusters[[k]], use.names = F))
#     name = file.path(getwd(), savePath, paste("group",k,"-size-",length(groupk),
#                                             ".txt", sep=''))
#     write.table(groupk, file=name, sep="\t", col.name=FALSE, row.names=FALSE)
#   }
#   return(savePath)
# }

# save_LOVE_result(clusters, delta, lbd)
# save_LOVE_result(res$I_ind, delta, lbd, T)




################################################################################
######                      Prediction   performance                      ######
################################################################################
X <- scale(X, T, T)

set.seed(20200309)

methods <- c("ER", "PFR", "PLS", "Lasso", "Perm-Lasso")
# methods <- c("PFR", "PLS", "Lasso", "Perm-Lasso")

# result <- replicate(50, KfoldCV(10, Y = Y, X = X, methods = methods,
#                                 pt_option = "perm_X_Y", metric = "mse", standardize = F,
#                                 delta_grid = 0.38, ER_sparse = F))

result <- replicate(50, KfoldCV(10, Y = Y, X = X, methods = methods,
                                pt_option = "perm_Y", metric = "adv", standardize = F,
                                delta_grid = 0.13, ER_sparse = T))

CV_result <- t(sapply(result['MSE',], function(v) v))
fitted_array <- sapply(result['fitted',], function(v) v, simplify = "array")

fitted_mat <- apply(fitted_array, 1:2, mean)

summary(CV_result, digits = 4)


# write.table(CV_result, file = "../output/age_dataset/all_days_CV_result_03_16_mse.txt", row.names = F, col.names = T)
# write.table(fitted_mat, file = "../output/age_dataset/all_days_CV_result_03_16_fit.txt", row.names = F, col.names = T)

CV_result <- read.table("../output/age_dataset/all_days_CV_result_03_16_mse_sparse.txt", header = T)
fitted_mat <- read.table("../output/age_dataset/all_days_CV_result_03_16_fit_sparse.txt", header = T)


CV_result <- read.table("../output/age_dataset/all_days_CV_result_03_16_mse.txt", header = T)
fitted_mat <- read.table("../output/age_dataset/all_days_CV_result_03_16_fit.txt", header = T)


# wilcox.test(CV_result$ER, CV_result$Perm.Lasso, alternative = "less")  ##  p-value < 2.2e-16

# sub_CV_result <- data.frame(sapply(1:ncol(CV_result), FUN = function(k, data) {
#   col_k <- data[,k]
#   col_k[which(col_k <= quantile(col_k, probs = 0.95))]
# }, data = CV_result))
# 
# names(sub_CV_result) <- names(CV_result)



# exploratory result
col_ind <- 1:length(methods)
summary(CV_result[,col_ind], digits = 3)
# making histograms of the MSEs for different methods
data_wide_CV <- data.frame(1 - CV_result[ ,col_ind])
names(data_wide_CV) <- methods
summary(data_wide_CV)
data_long_CV <- reshape(data_wide_CV, direction = "long", timevar = "Method",
                        varying = methods,
                        v.names = "MSE", times = methods)


data_long_CV$Method <- factor(data_long_CV$Method, 
                              levels = c("ER", "Lasso", "PFR", "PLS", "Perm-Lasso"))


ggplot(subset(data_long_CV, !Method %in% c("-Lasso")), aes(Method, MSE, color = Method, fill = Method)) +
  geom_violin(alpha = 0.5, scale = "area") + ylab("root mean square error") +
  xlab("") + theme_classic() + ylim(0.6, max(data_long_CV$MSE))



y_ranges <- range(as.vector(fitted_mat), Y)
par(mfrow=c(2,3))
order_inc <- order(Y)
corr <- c()
for (i in 1:length(methods)) {
  corr[i] <- cor(Y, fitted_mat[,i], method = "spearman")
  # if (i == 1) {
  plot(Y[order_inc], fitted_mat[order_inc,i], ylim = y_ranges, col = i, pch = i, xlim = y_ranges,
       ylab = "predicted Y", xlab = "observed Y", main = paste(methods[i]," (corr = ", round(corr[i],2), ")", sep = ""))
  # } else {
  # points(Y[order_inc], fitted_mat[order_inc,i], col = i, pch = i)
  # }
  # legend("topright", legend = methods, col = 1:length(methods), pch = 1:length(methods))
}
names(corr) <- methods


