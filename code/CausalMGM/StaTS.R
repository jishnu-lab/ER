library(rCausalMGM)
library(huge)

StaTS <- function(data, target, alphas=c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25), lambda=NULL, g=0.05, ns=20, leaveOneOut=F) {

    theta <- array(0, dim=c(ncol(data), ncol(data), length(alphas)))

    ## poss.mb.set <- c()
    
    if (leaveOneOut) {
        ns = nrow(data)
    }

    for (s in 1:ns) {
        if (leaveOneOut) {
            sub = setdiff(1:ns, s)
        } else {
            sub = sample(1:nrow(data), min(c(round(10*sqrt(nrow(data))), round(0.8*nrow(data)))), replace=F)
        }

        if (is.null(lambda)) {
            ig <- NULL
        } else {
            if (length(lambda) == 3) {
                ig <- mgm(data[sub,], lambda=lambda)
            } else {
                out.glasso <- huge(as.matrix(data[sub,]), method="glasso" , lambda=lambda)
                ig <- adjMat2Graph(out.glasso$icov[[1]], colnames(data[sub,]))

                ## out.glasso <- huge(as.matrix(YZ), method="glasso" , nlambda=30, lambda.min.ratio=0.1)
                ## out.select <- huge.select(out.glasso, criterion="stars", rep.num=20, stars.thresh=0.1)
                ## ig <- adjMat2Graph(out.select$opt.icov, colnames(data))
                ## rm(out.glasso)
                ## rm(out.select)
            }
            ## poss.mb.set <- union(poss.mb.set, ig$markov.blankets[[target]])

            ## for (neighbor in ig$markov.blankets[[target]]) {
            ##     poss.mb.set <- union(poss.mb.set, ig$markov.blankets[[neighbor]])

            ##     ## for (neighbor2 in ig$markov.blankets[[neighbor]]) {
            ##     ##     if (neighbor2 != target) {
            ##     ##         poss.mb.set <- union(poss.mb.set, ig$markov.blankets[[neighbor2]])
            ##     ##     }
            ##     ## }
            ## }
                
        }

        for (a in 1:length(alphas)) {
            graph <- fciMax(data[sub,], initialGraph=ig, alpha=alphas[a], fdr=F)
            ## print(graph$markov.blankets[[target]])
            ## for (i in 1:ncol(data)) {
            i <- which(colnames(data) %in% target)
            jdxs <- which(colnames(data) %in% graph$markov.blankets[[target]])
            for (j in jdxs) {
                theta[i,j,a] <- theta[i,j,a] + 1
            }
            ## }
        }
    }

    ## print(poss.mb.set[poss.mb.set!=target])

    theta <- theta / ns

    instabs <- as.vector(rep(0, length(alphas)))

    for (a in 1:length(alphas)) {
        ## for (i in 1:(ncol(data)-1)) {
        ##     for (j in 1:ncol(data)) {
        i <- which(colnames(data) %in% target)
        ## jdxs <- which(colnames(data) %in% graph$markov.blankets[[target]])
        for (j in 1:ncol(data)) {
            instabs[a] <- instabs[a] + 2 * theta[i,j,a] * (1-theta[i,j,a])
        }
        ##     }
        ## }
    }

    ## print(theta[i,,])
    ## print(rowSums(theta[i,,]))
    ## print(instabs)

    ## if (is.null(lambda)) {
    instabs <- instabs / (ncol(data)-1)
    ## } else {
    ##     instabs <- instabs / (length(poss.mb.set)-1)
    ## }
    

    for (a in 1:length(alphas)) {
        instabs[a] <- max(instabs[1:a])
    }
    
    list('alpha.hat' = alphas[which.min(sapply(instabs-g, abs))], 'alphas' = alphas, 'instabilities' = instabs)
}
