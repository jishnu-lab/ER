library(devtools)
install_github('tyler-lovelace1/rCausalMGM')
library(huge)
library(rCausalMGM)
library(tidyverse)
source('StaTS.R')

## FINAL MODEL: ALL LATENT FACTORS

data <- read.delim('../../output/age_dataset/data_mat.csv', sep=',', header=T)
Z <- read.delim('../../output/age_dataset/Z_mat.csv', sep=',', header=T)
colnames(Z) <- paste('Z', seq(ncol(Z)), sep='')
YZ <- cbind(Y=data$Y, Z)

## StARS model selection for GLASSO
set.seed(123)

out.glasso <- huge(as.matrix(YZ), method="glasso" , nlambda=30, lambda.min.ratio=0.1)

out.select <- huge.select(out.glasso, criterion="stars", rep.num=20, stars.thresh=0.1)

ig <- adjMat2Graph(out.select$opt.icov, colnames(YZ))

ig

print(paste("optimal lambda:", out.select$opt.lambda, sep=" "))

## StaTS model selection for FCI-Max

stats.select <- StaTS(YZ, 'Y', lambda=out.select$opt.lambda, g=0.05, leaveOneOut=T)

print(stats.select)

## $alpha.hat
## [1] 0.1

## $alphas
## [1] 0.01 0.05 0.10 0.15 0.20 0.25

## $instabilities
## [1] 0.009520348 0.040457571 0.048336479 0.048336479 0.048336479 0.048336479

## GLASSO with StARS selected lambda = 0.32

## final.glasso <- huge(as.matrix(YZ), method="glasso", lambda=0.32)

## ig <- adjMat2Graph(final.glasso$icov[[1]], colnames(YZ))

print(ig)
saveGraph(ig, '../../output/age_dataset/YZ_full_glasso.txt')

g <- fciMax(YZ, initialGraph=ig, alpha=stats.select$alpha.hat, verbose=T)

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

sif <- t(t(sif) %>% as.data.frame %>% mutate_if(function(x) x[2]=="o-o", function(x) x <- c(x[1], "cc", x[3])))

sif <- t(t(sif) %>% as.data.frame %>% mutate_if(function(x) x[2]=="o->", function(x) x <- c(x[1], "ca", x[3])))

sif <- t(t(sif) %>% as.data.frame %>% mutate_if(function(x) x[2]=="-->", function(x) x <- c(x[1], "dir", x[3])))

sif <- t(t(sif) %>% as.data.frame %>% mutate_if(function(x) x[2]=="<->", function(x) x <- c(x[1], "bidir", x[3])))
sif

saveGraph(g, '../../output/age_dataset/YZ_full_allLatent.txt')
write.table(sif, '../../output/age_dataset/YZ_full_allLatent.sif', sep='\t', row.names=F, col.names=F)

## FINAL MODEL: SIGNIFICANT LATENT FACTORS (CausER)

sig = c(1,2,5,7,9,15,17,19,24,36,41,42,44,45,46,49,52,53)
YZsig <- cbind(Y=data$Y, Z[,sig])

## StARS model selection for GLASSO

## out.glasso <- huge(as.matrix(YZsig), method="glasso" , nlambda=30, lambda.min.ratio=0.1)

## out.select <- huge.select(out.glasso, criterion="stars", rep.num=20, stars.thresh=0.1)

## print(paste("optimal lambda:", out.select$opt.lambda, sep=" "))

## GLASSO with StARS selected lambda = 0.32

final.glasso <- huge(as.matrix(YZsig), method="glasso" , lambda=0.32)

ig <- adjMat2Graph(final.glasso$icov[[1]], colnames(YZsig))

print(ig)
saveGraph(ig, '../../output/age_dataset/YZsig_full_glasso.txt')

g <- fciMax(YZsig, initialGraph=ig, alpha=0.1, verbose=T)

print(g)
print(g$edges)
saveGraph(g, '../../output/age_dataset/YZ_full_CausER.txt')
