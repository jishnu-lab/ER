library(rCausalMGM)
library(huge)

runtimes.undir <- data.frame(CPU.time=c(), Wall.clock.8cpus=c())
runtimes.causal <- data.frame(CPU.time=c(), Wall.clock.8cpus=c())
runtimes.total <- data.frame(CPU.time=c(), Wall.clock.8cpus=c())

#### Age ####

data <- read.csv('../../output/age_dataset/data_mat.csv', header=T)
Y <- data[,1]
X <- data[,-1]

Z <- read.csv('../../output/age_dataset/Z_mat.csv', header=T)
colnames(Z) <- paste('Z', seq(ncol(Z)), sep='')

YZ <- cbind(Y=Y, Z)


## Undirected ##

undir.time <- system.time({ out.glasso <- huge(as.matrix(YZ), method="glasso", lambda=0.31927)
                            ig <- adjMat2Graph(out.glasso$icov[[1]], colnames(YZ))})

runtimes.undir <- rbind(runtimes.undir,
                        data.frame(CPU.time=sum(undir.time[1:2]), Wall.clock.8cpus=undir.time[3]))

## Causal ##

causal.time <- system.time({ g <- fciMax(YZ, initialGraph=ig, alpha=0.05, verbose=F) })

runtimes.causal <- rbind(runtimes.causal,
                         data.frame(CPU.time=sum(causal.time[1:2]), Wall.clock.8cpus=causal.time[3]))



#### Tuberculosis ####

data <- read.csv('../../output/TB_data/data_mat.csv', header=T)
Y <- data[,1]
X <- data[,-1]

Z <- read.csv('../../output/TB_data/Z_mat.csv', header=T)
colnames(Z) <- paste('Z', seq(ncol(Z)), sep='')

YZ <- cbind(Y=Y, Z)


## Undirected ##

## system.time({ out.glasso <- huge(as.matrix(YZ), method="glasso", lambda=0.31927)
##     ig <- adjMat2Graph(out.glasso$icov[[1]], colnames(YZ))})

runtimes.undir <- rbind(runtimes.undir,
                        data.frame(CPU.time=NA, Wall.clock.8cpus=NA))

## Causal ##

causal.time <- system.time({ g <- fciMax(YZ, alpha=0.1, verbose=F) })

runtimes.causal <- rbind(runtimes.causal,
                         data.frame(CPU.time=sum(causal.time[1:2]), Wall.clock.8cpus=causal.time[3]))




#### Term / preterm ####

data <- read.csv('../../output/newborn_dataset/data_mat.csv', header=T)
Y <- data[,1]
X <- data[,-1]

Z <- read.csv('../../output/newborn_dataset/Z_mat.csv', header=T)
colnames(Z) <- paste('Z', seq(ncol(Z)), sep='')

YZ <- cbind(Y=Y, Z)


## Undirected ##

## system.time({ out.glasso <- huge(as.matrix(YZ), method="glasso", lambda=0.31927)
##     ig <- adjMat2Graph(out.glasso$icov[[1]], colnames(YZ))})

causal.time <- runtimes.undir <- rbind(runtimes.undir,
                                       data.frame(CPU.time=NA, Wall.clock.8cpus=NA))

## Causal ##

causal.time <- system.time({ g <- fciMax(YZ, alpha=0.2, verbose=F) })

runtimes.causal <- rbind(runtimes.causal,
                        data.frame(CPU.time=sum(causal.time[1:2]), Wall.clock.8cpus=causal.time[3]))




#### Malaria ####

data <- read.csv('../../output/malaria_dataset/data_mat.csv', header=T)
Y <- factor(data[,1])
X <- data[,-1]

Z <- read.csv('../../output/malaria_dataset/Z_mat.csv', header=T)
colnames(Z) <- paste('Z', seq(ncol(Z)), sep='')

YZ <- cbind(Y=Y, Z)

sig.idx <- c(12,28,34,35,43,57,91,104,116,119,130,138,143,154,161,178,181,188,193,203,211,212,219,231,232,235,241,245,256,269,270,271,283,288,292,295,308,309,335,357,364,376,387,388,390,393,398,452,457,460,468,469,475,479,483,493,494,500,503,508,521,525,528,529,534,536,537,539,541,543,552,562)

YZsig <- cbind(Y=Y, Z[,sig.idx])


## Undirected ##

undir.time <- system.time({ ig <- mgm(YZsig, lambda=0.27) })

runtimes.undir <- rbind(runtimes.undir,
                        data.frame(CPU.time=sum(undir.time[1:2]), Wall.clock.8cpus=undir.time[3]))

## Causal ##

causal.time <- system.time({ g <- fciMax(YZsig, initialGraph=ig, alpha=0.2, verbose=F) })

runtimes.causal <- rbind(runtimes.causal,
                         data.frame(CPU.time=sum(causal.time[1:2]), Wall.clock.8cpus=causal.time[3]))



rownames(runtimes.undir) <- c('age', 'TB', 'newborn', 'malaria')
rownames(runtimes.causal) <- c('age', 'TB', 'newborn', 'malaria')

runtimes.total <- runtimes.undir + runtimes.causal
runtimes.total[rowSums(is.na(runtimes.total))>0,] <- runtimes.causal[rowSums(is.na(runtimes.total))>0,]

runtimes.undir
runtimes.causal
runtimes.total

write.csv(round(runtimes.undir,3), '../../output/undir.runtime.csv', row.names=T)
write.csv(round(runtimes.causal,3), '../../output/causal.runtime.csv', row.names=T)
write.csv(round(runtimes.total,3), '../../output/total.runtime.csv', row.names=T)
