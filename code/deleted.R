# # Return the selected optimal turning parameter "K" for the given type of method by using k-folds Cross-Validation
# # args: k: k-fold CV; Y[n], X[n,p], 
# # return: k by length(methods) matrix: MSE for each method
# KfoldCV = function(k, Y, X, methods) {
#   kFolds = splitData(k, Y, X)
#   MSEs = getFoldMSE(kFolds, methods)
#   apply(MSEs, 2, mean)
# }
# 
# 
# # Split the given {@code Y} and {@code X} into k-folds and each fold contains list(Y,X).
# # args: k: k-fold CV
# # args: Y[n] and X[n,p]
# # return: list of k folds, each element contains a list of Y_k[n/k] and X_k[n/k,p]
# splitData = function(k, Y, X) {
#   n = nrow(X)
#   randomIndices = sample(1:n,replace = F)
#   indicesPerGroup = extract(randomIndices,partition(n,k))
#   KfoldsX = KfoldsY = vector("list",length=k)
#   for (i in 1:k) {
#     KfoldsX[[i]] = X[indicesPerGroup[[i]],]
#     KfoldsY[[i]] = Y[indicesPerGroup[[i]]]
#   }
#   return(list(Y=KfoldsY,X=KfoldsX))
# }
# 
# # For each fold, get the training and validation data. Then for each turning parameter K, find the fitted A_hat and 
# # calculate the MSE on the validation set.
# # args: kFolds: the list returned by {@code splitData}
# # return: k by length(K_grid) matrix containing the MSE
# getFoldMSE = function(kFolds, methods) {
#   k = length(kFolds$X)
#   MSEs <- c()
#   for (j in 1:k) {
#     validY = kFolds$Y[[j]]
#     validX = kFolds$X[[j]]
#     trainX = trainY = NULL
#     for (i in setdiff(1:k, j)) {
#       trainX = rbind(trainX, kFolds$X[[i]])
#       trainY = c(trainY, kFolds$Y[[i]])
#     }
#     MSEs_j <- c()
#     for (l in 1:length(methods)) 
#       MSEs_j = c(MSEs_j, calMSE(validY, validX, trainX, trainY, methods[l]))
#     MSEs <- rbind(MSEs, MSEs_j)
#   }
#   return(MSEs)
# }
















#### auc


CV_binary <- function(Y, X, methods = c("Lasso", "PFR", "ER", "PLS"), 
                      N_split = 100) {
  auc <- matrix(NA, length(methods), N_split)
  for (j in 1:N_split) {
    train_ind <- sample(1:nrow(X), floor(nrow(X) * 0.7), replace = F)
    for (i in 1:length(methods)) {
      auc[i,j] <- calAUC(Y[-train_ind], X[-train_ind, ], X[train_ind, ], 
                         Y[train_ind], methods[i])
    }
  }
  rownames(auc) <- methods
  return(auc)
}




set.seed(20191201)
methods = c("Lasso", "PFR", "ER", "PLS", "Perm-ER")
errors <- CV_binary(Y, X, methods = methods, N_split = 100)


# write.table(t(errors), file = "../output/AUC_result.txt", row.names = F, col.names = T)
errors <- read.table("../output/AUC_result.txt", header = T)


data_wide <- data.frame(errors)
names(data_wide) <- methods
data_long <- reshape(data_wide, direction = "long", timevar = "Method",
                        varying = methods,
                        v.names = "AUC", times = methods)

# ggplot(data_long_CV, aes(x = MSE, color = Method, fill = Method)) +
#   geom_histogram(alpha = 0.5, position="identity", bins = 30) +
#   theme(legend.position = "top")


ggplot(data_long, aes(Method, AUC, color = Method, fill = Method)) +
  geom_violin(alpha = 0.5) + theme(legend.position = "top") + 
  ylab("Area under the curve")









