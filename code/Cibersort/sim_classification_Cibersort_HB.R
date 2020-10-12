
rm(list=ls())
setwd("~/Documents/Mike/Projects/Application of Essential Regression/code/Code for ER")
source("SupLOVE.R")
setwd("~/Documents/Mike/Projects/Application of Essential Regression/code")
source("K-CV.R")
library(ggplot2)
library(readr)
library(ROCR)

HB_data <- read_delim("../dataset/CIBERsortdata/Cibersort_HB_normal.txt", 
                                  "\t", escape_double = FALSE, col_names = FALSE, 
                                  trim_ws = TRUE)
# View(HB_data)

Y <- HB_data$X1
X <- data.frame(HB_data[,-1])


################################################################################
#######                       ER on the whole dataset                    #######
################################################################################

X <- scale(X, T, T)
Y <- Y - mean(Y)

# set.seed(20200316)
# delta_grid <- seq(0.1, 0.3, 0.01)
# 
# res <- ER(Y, X, delta_grid, pred = F, beta_est = "NULL", rep_CV = 100, verbose = T, merge = F)
# 
# res$K
# selected delta = 0.17

lbd = 0.5
delta = 0.17

res <- ER(Y, X, delta = delta, pred = F, beta_est = "LS", rep_CV = 50, verbose = T,
          merge = F, CI = T, correction = "Bonferroni", lbd = lbd)
res$K
round(res$beta_CIs, 3)
p_vals <- 2 * pnorm(abs(res$beta) / sqrt(res$beta_var / length(Y)), lower.tail = F)

clusters <- recoverGroup(res$A)

sapply(clusters, function(x) as.vector(sort(unlist(x))))
sapply(res$I_ind, function(x) as.vector(unlist(x)))


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

save_LOVE_result(clusters, delta, lbd)
save_LOVE_result(res$I_ind, delta, lbd, T)









################################################################################
######                    Classification performance                      ######
################################################################################
X <- scale(X, T, T)

methods = c("ER", "PFR", "PLS", "Lasso", "Perm-Lasso")


# set.seed(20200316)

# pred_val <- CV_binary(Y, X, methods, seq(0.1, 0.3, 0.01), F)

# write.table(pred_val, file = "../output/Cibersort_dataset/pred_val_binary_HB.txt", row.names = F, col.names = T)


pred_val <- read.table("../output/Cibersort_dataset/pred_val_binary_HB.txt", header = T)
colnames(pred_val) <- methods

true_label <- Y

result <- Compute_metric(pred_val, true_label)

# auc
names(result$auc) <- methods
round(result$auc, 3)
#  ER        PFR        PLS      Lasso Perm-Lasso
# 0.92       0.90       0.93       0.86       0.76


Plot_ROC(result$error, methods, "bottomright")
Plot_ROC(result$error, setdiff(methods, "Perm-Lasso"), "bottomright")








################################################################################
######              Three-way classification performance                  ######
################################################################################
library(e1071)
# standardize the features
X <- scale(X, T, T)

set.seed(20200320)

methods = c("ER", "PFR", "PLS", "Lasso", "Perm-Lasso")
delta_grid <- seq(0.1, 0.3, 0.01)

result <- replicate(50, KfoldCV_tenary(5, Y, X, methods, Y_type = "factor",
                                      pt_option = "perm_Y_before_split",
                                      delta_grid = delta_grid,
                                      ER_sparse = F))

CV_result <- t(sapply(result['metric',], function(v) v))
pred_array <- sapply(result['pred',], function(v) v, simplify = "array")

summary(CV_result, digits = 3)

# write.table(CV_result, file = "../output/Cibersort_dataset/Accuracy_non_sparse.txt", row.names = F, col.names = T)


CV_result <- read.table("../output/Cibersort_dataset/Accuracy_non_sparse.txt", header = T)

# exploratory result
col_ind <- 1:length(methods)
summary(CV_result[,col_ind], digits = 3)
# making histograms of the MSEs for different methods
data_wide_CV <- data.frame(CV_result[ ,col_ind])
names(data_wide_CV) <- methods
summary(data_wide_CV)
data_long_CV <- reshape(data_wide_CV, direction = "long", timevar = "Method",
                        varying = methods,
                        v.names = "MSE", times = methods)
data_long_CV$Method <- factor(data_long_CV$Method, levels = c("ER", "Lasso", "PFR", "PLS", "Perm-Lasso"))

ggplot(data_long_CV, aes(Method, MSE, color = Method, fill = Method)) +
  geom_boxplot(alpha = 0.5) + ylab("Accuracy") + ylim(c(0.25, 1)) + theme_classic()



