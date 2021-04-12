precisionMat2txt <- function(theta, fname='glasso.txt', nodes=NULL) {
    if (is.null(nodes)) {
        nodes <- colnames(theta)
    }
    
    cat("Graph Nodes:\n", file=fname)
    cat(nodes, file=fname, append=TRUE, sep=",")
    cat("\n\nGraph Edges:\n", file=fname, append=TRUE)

    count <- 1
    for (i in 1:(nrow(theta)-1)) {
        for (j in (i+1):nrow(theta)) {
            if (theta[i,j] != 0) {
                cat(sprintf("%d. %s --- %s\n", count, nodes[i], nodes[j]), file=fname, append=TRUE)
                count <- count + 1
            }
        }
    }
}
