####### This script contains functions to estimate the pure node set ##########
source("Utilities.R")


FindPureNode = function(off_Sigma, delta, Ms, arg_Ms, se_est, merge) {
  # Estimate list of pure node indices for given {@code Sigma} and {@delta}.
  # 
  # Args: 
  #   Sigma: p by p sample covariance matrix.
  #   delta: numerical constant.
  #   Ms: the largest absolute values of each row of Sigma.
  #
  # Returns: 
  #   a list including:
  #     the list of the estimated pure node indices
  #     the vector of the estimated pure node indices
  
  G <- list()

  for (i in 1:nrow(off_Sigma)) {
    row_i <- off_Sigma[i,]
    Si <- FindRowMaxInd(i, Ms[i], arg_Ms[i], row_i, delta, se_est)
    if (length(Si) != 0) {
      pureFlag <- TestPure(row_i, i, Si, Ms, arg_Ms, delta, se_est)
      # pureFlag <- TestPure_new(off_Sigma, i, Si, Ms, arg_Ms, delta, se_est)
      if (pureFlag) {
        if (merge)
          G <- Merge(G, c(Si, i))
        else
          G <- Merge_union(G, c(Si, i))
      }
    }
  }
  return(list(pureInd = G, pureVec = unlist(G)))
}

FindRowMax <- function(Sigma) {
  # Calculate the maximal absolute value for each row of the given matrix.
  #
  # Args:
  #   Sigma: p by p matrix
  # 
  # Returns:
  #   length p vector
  p <- nrow(Sigma)
  M <- arg_M <- rep(0, p)  
  for (i in 1:p) {
    row_i <- Sigma[i,]
    arg_M[i] <- which.max(row_i)
    M[i] <- row_i[arg_M[i]]
  }
  return(list(arg_M = arg_M, M = M))
}

FindRowMaxInd <- function(i, M, arg_M, vector, delta, se_est) {
  # Calculate indices of ith row such that the absolute values of these indices
  # are within 2 * delta from the maximal absolute value {@code M} of this row.
  # 
  # Args:
  #   i: integer denoting for the ith row.
  #   M: the maximal absolute value of the ith row.
  #   vector: the entries of this row.
  #   delta: numerical constant.
  #
  # Returns:
  #   a vector of indices.
  lbd <- delta * se_est[i] * se_est[arg_M] + delta * se_est[i] * se_est 
  indices <- which(M <= lbd + vector)
  return(indices)
}

TestPure <- function(Sigma_row, rowInd, Si, Ms, arg_Ms, delta, se_est) {
  # For given row, check if it is a pure node by iteratively checking the nodes
  # in Si. Return TRUE if the given row corresponds to a pure variable.
  #
  # Args:
  #   Sigma: p by p matrix.
  #   rowInd: integer index.
  #   Si: vector of indices.
  #   Ms: vector of largest absolute values of each rows in Si
  #   delta: numerical constant.
  # Returns:
  #   TRUE or FALSE.
  for (i in 1:length(Si)) {
    j <- Si[i]
    delta_j <- (se_est[rowInd] + se_est[arg_Ms[j]]) * se_est[j] * delta
    if (abs(Sigma_row[j] - Ms[j]) > delta_j)
      return(FALSE)
  }
  return(TRUE)
}

# TestPure_new <- function(off_Sigma, rowInd, Si, Ms, arg_Ms, delta, se_est) {
#   # For given row, check if it is a pure node by iteratively checking the nodes
#   # in Si. Return TRUE if the given row corresponds to a pure variable.
#   #
#   # Args:
#   #   Sigma: p by p matrix.
#   #   rowInd: integer index.
#   #   Si: vector of indices.
#   #   Ms: vector of largest absolute values of each rows in Si
#   #   delta: numerical constant.
#   # Returns:
#   #   TRUE or FALSE.
#   group_ind <- c(Si, rowInd)
#   for (i in 1:length(Si)) {
#     j <- Si[i]
#     delta_j <- (se_est[group_ind][-i] + se_est[arg_Ms[j]]) * se_est[j] * delta
#     if (sum(abs(off_Sigma[group_ind, j][-i] - Ms[j]) > delta_j) > 0)
#       return(FALSE)
#   }
#   return(TRUE)
# }



FindSignPureNode <- function(pureList, Sigma) {
  # Estimate the sign subpartition of pure node sets. If there is an element 
  # of a list is empty, then a empty list will be put in that position
  # 
  # Args: 
  #   pureList: list of pure node indices (Example: list(c(1,2,3),c(4,5,6,7)))
  #   Sigma: p by p sample covariance
  #
  # Returns: 
  #   list of sign subpartition of pure node indices.
  #   Example: list(list(c(1,2),3),list(c(4,7),c(5,6)))
  signPureList <- list()
  for (i in 1:length(pureList)) {
    purei <- sort(pureList[[i]])   ### For simulation purpose only.
    if (length(purei) != 1) {
      firstPure <- purei[1]
      pos <- firstPure
      neg <- c()
      for (j in 2:length(purei)) {
        purej <- purei[j]
        if (Sigma[firstPure, purej] < 0)
          neg <- c(neg, purej)
        else
          pos <- c(pos, purej)
      }
      if (length(neg) == 0)
        neg <- list()
      signPureList[[i]] <- list(pos = pos, neg = neg)
    } else
      signPureList[[i]] <- list(pos = purei, neg = list())
  }
  return(signPureList)
}

Merge <- function(groupList, groupVec) {
  # merge the new group with the previous ones which have common nodes
  if (length(groupList) != 0) {
    for (i in 1:length(groupList)) {
      common_nodes <- intersect(groupList[[i]], groupVec)
      if (length(common_nodes) != 0) {
        groupList[[i]] <- common_nodes
        return(groupList)
      }
    }
  } 
  groupList <- append(groupList, list(groupVec))
  return(groupList)
}

Merge_union <- function(groupList, groupVec) {
  # merge the new group with the previous ones which have common nodes
  if (length(groupList) != 0) {
    common_groups <- sapply(groupList, FUN = function(x, y) {
      length(intersect(x, y))
    }, y = groupVec)
    common_inds <- which(common_groups > 0)
    if (length(common_inds) > 0){
      new_group <- unlist(lapply(common_inds, 
                                        FUN = function(x, y){y[[x]]}, y = groupList))
      remain_group <- lapply(which(common_groups == 0), 
                             FUN = function(x, y){y[[x]]}, y = groupList)
      groupList <- append(remain_group, list(union(groupVec, new_group)))
      return(groupList)
    }
  } 
  groupList <- append(groupList, list(groupVec))
  return(groupList)
}



RecoverAI <- function(estGroupList, p) {
  # Recover the estimated submatrix A_I by given the pure node group.
  # 
  # Args:
  #   estGroupList: list of group indices of the pure node with sign.
  #
  # Returns: 
  #   p by K matrix.
  K <- length(estGroupList)
  A <- matrix(0, p, K)
  for (i in 1:K) {
    groupi <- estGroupList[[i]]
    A[groupi[[1]],i] <- 1
    groupi2 <- groupi[[2]]
    if (length(groupi2) != 0) 
      A[groupi2, i] <- -1
  }
  return(A)
}

EstAI <- function(Sigma, optDelta, se_est, merge) {
  # Use the given {@code optDelta} to calculate the fitted AI, pure variable 
  # indices in list form and vector form. Also return estimated Y and C for 
  # the following Dantzig estimation.
  #
  # Args: 
  #   Sigma: p by p covariance matrix.
  #   optDelta: optDelta to be used.
  #
  # Return: 
  #   the list including:
  #     AI: the p by K estimated AI
  #     pureVec: vector of the indices of estimated pure variables
  #     pureSignInd: list of the indices of estimated pure variables
  off_Sigma <- abs(Sigma)
  diag(off_Sigma) <- 0
  result_Ms <- FindRowMax(off_Sigma)
  Ms <- result_Ms$M
  arg_Ms <- result_Ms$arg_M
  
  resultPure <- FindPureNode(off_Sigma, optDelta, Ms, arg_Ms, se_est, merge)
  estPureIndices <- resultPure$pureInd
  estPureVec <- resultPure$pureVec
  
  estSignPureIndices <- FindSignPureNode(estPureIndices, Sigma)
  AI <- RecoverAI(estSignPureIndices, nrow(off_Sigma))
  
  return(list(AI = AI, pureVec = estPureVec, pureSignInd = estSignPureIndices))
}

