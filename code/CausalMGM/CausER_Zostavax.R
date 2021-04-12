library(devtools)
library(huge)
install_github('tyler-lovelace1/rCausalMGM', ref='fci-max')
library(rCausalMGM)

## FINAL MODEL: ALL LATENT FACTORS

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
print(g$edges)
saveGraph(g, '../../output/age_dataset/YZ_full_allLatent.txt')

## FINAL MODEL: SIGNIFICANT LATENT FACTORS (CausER)

sig = c(1,2,5,7,9,15,17,19,24,36,41,42,44,45,46,49,52,53)
YZsig <- cbind(Y=data$Y, Z[,sig])

out.glasso <- huge(as.matrix(YZsig), method="glasso" , nlambda=30, lambda.min.ratio=0.1)

out.select <- huge.select(out.glasso, criterion="stars", rep.num=20, stars.thresh=0.1)

print(paste("optimal lambda:", out.select$opt.lambda, sep=" "))

ig <- adjMat2Graph(out.select$opt.icov, colnames(YZsig))

print(ig)
saveGraph(ig, '../../output/age_dataset/YZsig_full_glasso.txt')

g <- fciMax(YZsig, initialGraph=ig, alpha=0.1, verbose=T)

print(g)
print(g$edges)
saveGraph(g, '../../output/age_dataset/YZ_full_CausER.txt')
