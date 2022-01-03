library(devtools)
install_github('tyler-lovelace1/rCausalMGM')
library(rCausalMGM)
library(tidyverse)
setwd("../Code\ for\ ER/")
source("SupLOVE.R")
setwd("../")
source("K-CV.R")
library(readr)
library(ROCR)
library(e1071)
library(doParallel)
library(foreach)
source("Helper.R")
source("Other_algorithms.R")
setwd("CausalMGM/")
source('StaTS.R')

data <- read.csv('../../output/TB_data/data_mat.csv', header=T)
Y <- factor(round(data[,1]) + 1, levels=c(0,1))
X <- data[,-1]

Z <- read.csv('../../output/TB_data/Z_mat.csv', header=T)
colnames(Z) <- paste('Z', seq(ncol(Z)), sep='')

YZ <- cbind(Y=Y, Z)

stats.select <- StaTS(YZ, 'Y', g=0.1, leaveOneOut=T)
print(stats.select)

## $alpha.hat
## [1] 0.1

## $alphas
## [1] 0.01 0.05 0.10 0.15 0.20 0.25

## $instabilities
## [1] 0.00000000 0.01611111 0.08888889 0.18583333 0.19500000 0.19500000


g <- fciMax(YZ, alpha=stats.select$alpha.hat, fdr=F, verbose=T)
g$edges

print(g)
print(g$edges)
print(g$markov.blankets[['Y']])
edges <- str_split(g$edges, " ")
sif <- data.frame(var1 = edges[[1]][1],
                  edge = edges[[1]][2],
                  var2 = edges[[1]][3])
for (i in 2:length(edges)) {
    sif <- rbind(sif, edges[[i]])
}

sif <- t(t(sif) %>% as.data.frame %>%
         mutate_if(function(x) x[2]=="o-o", function(x) x <- c(x[1], "cc", x[3])))

sif <- t(t(sif) %>% as.data.frame %>%
         mutate_if(function(x) x[2]=="o->", function(x) x <- c(x[1], "ca", x[3])))

sif <- t(t(sif) %>% as.data.frame %>%
         mutate_if(function(x) x[2]=="-->", function(x) x <- c(x[1], "dir", x[3])))

sif <- t(t(sif) %>% as.data.frame %>%
         mutate_if(function(x) x[2]=="<->", function(x) x <- c(x[1], "bidir", x[3])))

rownames(sif) <- NULL
sif

saveGraph(g, '../../output/TB_data/YZ_full_allLatent.txt')
write.table(sif, '../../output/TB_data/YZ_full_allLatent.sif', sep='\t', row.names=F, col.names=F)


############################################################
#####            Run ER on whole dataset              ######
############################################################


data <- read.csv('../../output/TB_data/data_mat.csv', header=T)
Y <- data[,1]
X <- data[,-1]

X <- scale(X, T, T)

set.seed(20200320)

delta <- 0.35
lbd <- 0.5

res <- ER(Y, X, delta = delta, pred = T, beta_est = "LS", rep_CV = 50, verbose = T,
          merge = F, CI = T, correction = NULL, lbd = lbd)
res$K

p_vals <- 2 * pnorm(abs(res$beta) / sqrt(res$beta_var / length(Y)), lower.tail = F)
round(p_vals, 3)

Z_hat <- X %*% res$A %*% solve(crossprod(res$A))
Z_tilde <- Pred_Z_BLP(X, res$A, res$C, res$Gamma, 1:res$K)

############################################################
#####                      K-fold                     ######
############################################################

data <- read.csv('../../output/TB_data/data_mat.csv', header=T)
Y <- data[,1]
X <- data[,-1]

X <- scale(X, T, T)

set.seed(20200320)

methods = c("ER", "CausER", "ER-Lasso-Dantzig", "PLS", "PFR", "Lasso", "Perm-Lasso")
# methods = "ER-Lasso-Dantzig"

delta_grid <- 0.35

no_cores <- 4
registerDoParallel(cores=no_cores) 
cl <- makeCluster(no_cores)

result <- foreach(i = 1:50, .combine = "rbind") %dopar%
  KfoldCV_tenary(10, Y, X, methods, Y_type = "factor", pt_option = "perm_Y_before_split", 
                 delta_grid = delta_grid, alpha=stats.select$alpha.hat, ER_sparse = F)
stopCluster(cl)  

CV_result = t(sapply(result[,'metric'], function(v) v))
# pred_array <- sapply(result[,'pred'], function(v) v, simplify = "array")

write.table(CV_result, file = "../../output/TB_data/CV_result.txt", row.names = F, col.names = T)


methods = c("ER", "CausER", "Composite", "PLS", "PFR", "Lasso", "Perm-Lasso")

CV_result = read.table("../../output/TB_data/CV_result.txt", header = T)

pdf('../../output/TB_data/repcv_violinplot.pdf', height = 5.3, width = 5.7)
Violin_plot(methods[c(2,1,3,4,5,6,7)], CV_result[,c(2,1,3,4,5,6,7)], MSE_flag = F, label_names = methods[c(2,1,3,4,5,6,7)], ylims = c(0.2, 0.9))
dev.off()

methodlist <- c()
for (m in colnames(CV_result))
    methodlist <- c(methodlist, rep(m, 50))
sigtestdf <- data.frame(accuracy = unlist(CV_result, use.names=F),
                        method = methodlist)

## kruskal.test(accuracy ~ method, sigtestdf[sigtestdf$method!='Perm.Lasso',])

## kruskal.test(accuracy ~ method, sigtestdf[sigtestdf$method %in% unique(methodlist)[c(1:2, 6)],])


#### LOOCV 

X <- scale(X, T, T)

set.seed(20191205)

methods = c("ER", "CausER", "ER-Lasso-Dantzig", "PFR", "PLS", "Lasso", "Perm-Lasso")

# result <- KfoldCV_tenary(10, Y, X, methods, Y_type = "factor", pt_option = "perm_Y_before_split", delta_grid = 0.04)
# pred_array <- result$pred


pred_array <- KfoldCV_tenary(nrow(X), Y, X, methods, Y_type = "factor", pt_option = "perm_Y_before_split", alpha=0.1, delta_grid = 0.35, prob_flag=T)

#### ACCURACY ####
## ER            CausER        ER-Lasso-Dantzig    PFR 
## 0.7000000     0.7666667     0.6333333           0.7666667
## PLS            Lasso        Perm-Lasso 
## 0.7666667      0.7666667    0.6333333 

# write.table(pred_array, file = "pred_result_CV10.txt", row.names = F, col.names = F)
write.table(pred_array$pred, file = "pred_result_LOOCV.txt", row.names = F, col.names = methods)
write.table(pred_array$prob, file = "prob_result_LOOCV.txt", row.names = F, col.names = methods)

pred_array <- read.table("pred_result_LOOCV.txt", header=T)
prob_array <- read.table("prob_result_LOOCV.txt", header=T)


methods[3] <- "Composite"

colnames(prob_array) <- methods

true_label <- Y
result <- Compute_metric(prob_array[,c(2,1,3,4,5,6,7)], true_label)

# auc
names(result$auc) <- methods[c(2,1,3,4,5,6,7)]
round(result$auc, 3)

#### AUC ####
## ER      CausER    Composite    PFR        PLS      Lasso      Perm-Lasso 
## 0.715   0.695     0.350        0.860      0.845    0.720      0.550


pdf('../../output/TB_data/loocv_rocplot.pdf', height = 5.3, width = 5.3)
Plot_ROC(result$error, setdiff(methods[c(2,1,3,4,5,6,7)], c("Composite", "Perm-Lasso")), aucs=result$auc, "bottomright")
dev.off()
