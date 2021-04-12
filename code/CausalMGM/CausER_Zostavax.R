library(devtools)
library(huge)
install_github('tyler-lovelace1/rCausalMGM', ref='fci-max')
library(rCausalMGM)

data <- read.delim('../../output/age_dataset/data_mat.csv', sep=',', header=T)
Z <- read.delim('../../output/age_dataset/Z_mat.csv', sep=',', header=T)
YZ <- cbind(Y=data$Y, Z)

out.glasso <- huge(as.matrix(YZ), method="glasso" , nlambda=30, lambda.min.ratio=0.1)

out.select <- huge.select(out.glasso, criterion="stars", rep.num=20, stars.thresh=0.1)

print(paste("optimal lambda:", out.select$opt.lambda, sep=" "))

ig <- adjMat2Graph(out.select$opt.icov, colnames(YZ))

print(ig)
saveGraph(ig, '../../output/age_dataset/YZ_full_glasso.txt')

g <- fciMax(YZ, initialGraph=ig, alpha=0.1, verbose=T)

print(g)
saveGraph(g, '../../output/age_dataset/YZ_full_CausER.txt')
