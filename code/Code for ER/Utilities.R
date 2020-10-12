######### This script contains some utilities function to test the performance ########
library(Matrix)
library(gtools)
library(glmnet)

checkElement <- function(element, groupList) {
  # Check if an element is in a list. If it does exist in some group, return the group
  # index in that list and its sublist index. Otherwise, return c(0,0).
  for (i in 1:length(groupList)) {
    for (j in 1:length(groupList[[i]])) {
      if (element %in% groupList[[i]][[j]])
        return(c(i,j))
    }
  }
  return(c(0,0))
}

signPerm <- function(n) {
  # For a given integer K, generate its sign permutation
  # Input: an integer K
  # Return: a 2^K by K matrix containing all sign permutations
  signPerms <- rep(1, n)
  for (i in 1:n) {
    allCombs <- combinations(n, i, 1:n)
    for (j in 1:nrow(allCombs)) {
      thisPerm <- rep(1, n)
      thisPerm[allCombs[j, ]] <- -1
      signPerms <- rbind(signPerms, thisPerm)
    }
  }
  return(signPerms)
}

pureRowInd <- function(A) {
  # For given matrix A, find the row indices which correspond to pure nodes
  #
  # Args: 
  #   A: p by K matrix.
  #
  # Returns: 
  #   vector of pure row indices
  pureVec <- c()
  for (i in 1:ncol(A)) {
    pureVec <- c(pureVec, which(abs(A[ ,i]) == 1))
  }
  return(pureVec)
}

getOptPerm <- function(A, B) {
  # For the given pure node set, find the best column permutation matrix 
  # Aperm of A such that |Aperm - B|_F is minimal for the target B
  #
  # Args: 
  #   A: generic matrix.
  #   B: generic matrix with the same dimension of A.
  #
  # Returns: 
  #   list of column permutation and sign permutation
  K <- ncol(A)
  allPerm <- permutations(K, K, 1:K)
  signPerms <- signPerm(K)
  prevLoss <- norm(A - B, "F")
  optPerm <- vector("list", length <- 2) 
  for (i in 1:nrow(allPerm)) {
    # skip the identity permutation
    permi <- as(as.integer(allPerm[i, ]), "pMatrix")
    permutatedA <- A %*% permi
    for (j in 1:nrow(signPerms)) {
      signPermj <- signPerms[j, ]
      # perform sign permutation along columns of newA
      newA <- t(t(permutatedA) * signPermj)
      currLoss <- norm(newA - B, "F")
      if (currLoss <= prevLoss) {
        optPerm[[1]] <- permi
        optPerm[[2]] <- signPermj
        prevLoss <- currLoss
      }
    }
  }
  return(optPerm)
}

permA <- function(A, B) {
  # Find the optimal permutated matrix of A such that |Aperm-B|_F is minimized.
  #
  # Args: 
  #   A,B: two matrices with the same dimension. 
  # Return: 
  #   Permutated A.
  if (sum(dim(A) == dim(B)) == 2) {
    pureInd <- pureRowInd(B)
    optPerm <- getOptPerm(A[pureInd, ], B[pureInd, ])
    permutatedA <- A %*% optPerm[[1]]
    A <- t(t(permutatedA) * optPerm[[2]])
  }
  return(A)
}

calSPandSN <- function(p, estGroup, trueGroup) {
  # Calculate specificity and sensitivity based on group partition
  #
  # Args:
  #   p: the total number of variables
  #   estGroup(trueGroup): list of variable indices for each group
  #
  # Returns:
  #   Specificity and Sensitivity
  TP <- TN <- FP <- FN <- TPS <- TNS <- FNS <- FPS <- 0
  estTable <- trueTable <- matrix(0, nrow=p, ncol=2)
  for (k in 1:p) {
    estTable[k,] <- checkElement(k,estGroup)
    trueTable[k,] <- checkElement(k,trueGroup)
  }
  for (i in 1:(p-1)) {
    flagiEst <- estTable[i,]
    flagiTrue <- trueTable[i,]
    for (j in (i+1):p) {
      flagjEst <- estTable[j,]
      flagjTrue <- trueTable[j,]
      if (flagiTrue[1] * flagjTrue[1] > 0 && flagjTrue[1] == flagiTrue[1]) {
        if (flagiEst[1] * flagjEst[1] > 0 && flagiEst[1] == flagjEst[1]) {
          TP <- TP + 1
          if (flagiTrue[2] == flagjTrue[2]) {
            if (flagiEst[2] == flagjEst[2]) {
              TPS <- TPS + 1
            } else {
              FNS <- FNS + 1
            }
          } else {
            if (flagiEst[2] == flagjEst[2]) {
              FPS <- FPS + 1
            } else {
              TNS <- TNS + 1
            }
          }
        } else {
          FN <- FN + 1
        }
      } else {
        if (flagiEst[1] * flagjEst[1] > 0 && flagiEst[1] == flagjEst[1]) {
          FP <- FP + 1
        } else {
          TN <- TN + 1
        }
      }
    }
  }
  ifelse (TN + FP == 0, SP <- 0, SP <- TN / (TN + FP))
  ifelse (TP + FN == 0, SN <- 0, SN <- TP / (TP + FN))
  ifelse (FNS + TPS == 0, SNS <- 0, SNS <- TPS / (FNS + TPS))
  ifelse (TNS+FPS == 0, SPS <- 0, SPS <- TNS / (FPS + TNS))
  return(c(SP = SP, SN = SN, SPS = SPS, SNS = SNS))
}

recoverGroup <- function(A) {
  # Recover group structure based given p by K matrix A and perform thresholding
  #
  # Args: 
  #   A: estimated matrix.
  #   thresh: constant > 0.
  # 
  # Returns: 
  #   list of group indices with sign subpartition
  Group <- list()
  for (i in 1:ncol(A)) {
    column <- A[,i]
    posInd <- which(column > 0)
    negInd <- which(column < 0)
    Group[[i]] <- list(pos = posInd, neg = negInd)
  }
  return(Group)
}

calFNRandFPR <- function(estA,trueA) {
  # Calculate the FNR and FPR for given p by K matrix A.
  #
  # Args: 
  #   estA(trueA): two matrices having same dims.
  #
  # Returns: 
  #   a list including:
  #     list of error: FNR and FPR.
  #     errorInd: the row indices where errors occur.
  if (sum(dim(estA) == dim(trueA)) != 2) {
    # If estA and A don't have the same dimension, we return -1.
#    cat("Two matrices don't have the same dimension!\n")
    return(list(error    = c(FPR = -1, FNR = -1, FPSR=-1, FNSR=-1), 
                errorInd = c()))
  }
  errorInd <- c()
  p <- nrow(estA)
  K <- ncol(estA)
  TP <- FP <- TN <- FN <- FPS <- TPS <- TNS <- FNS <- 0
  for (i in 1:p) {
    for (j in 1:K) {
      if (trueA[i,j] != 0) {
        TN <- TN + 1
        if (estA[i,j] == 0) {
          errorInd <- c(errorInd, i) # record the row indices where errors occur
          FN <- FN + 1
        } else {  ### calculate the sign error rate
          if (trueA[i,j] > 0) {
            TPS <- TPS + 1
            if (estA[i,j] < 0) 
               FNS <- FNS + 1
          } else {
            TNS <- TNS + 1
            if (estA[i,j] > 0)
              FPS <- FPS + 1
          }
        }
      } else {
        TP <- TP + 1
        if (estA[i,j] != 0) {
          FP <- FP + 1
          errorInd <- c(errorInd, i)
        }
      }
    }
  }
  if (TPS == 0) 
    FNSR <- 0
  if (TNS == 0) 
    FPSR <- 0
  return(list(error = c(FPR = FP / TP, FNR = FN / TN,FPSR = FPS / TNS, FNSR = FNS / TPS),
              errorInd = unique(errorInd)))
}

calArates <- function(A) {
  # Calculate the scaled matrix l_1 norm and frobenius norm.
  #
  # Args: 
  #   A: matrix
  # 
  # Return: 
  #   vector of three numerical results of c(frob, sup, L_inf).
  frob <- norm(A, "F") / sqrt(ncol(A) * nrow(A))
  l1 <- sum(abs(A)) / nrow(A) / ncol(A)
  return(c(l1 = l1, l2 = frob))
}

singleton <- function(estPureIndices) {
  # Check if there exists an element of the given list has length equal to 1
  # If exists at least one, return TRUE; otherwise return FALSE
  if (length(estPureIndices) == 0)
    return(T)
  else
    ifelse(sum(sapply(estPureIndices, FUN = function(x) {length(x)}) == 1) > 0, T, F)
}

getErrors <- function(estA, trueA, flagPerm = TRUE) {
  estGroup <- recoverGroup(estA)
  trueGroup <- recoverGroup(trueA)
  speAndSen <- calSPandSN(nrow(trueA), estGroup, trueGroup)
  if (flagPerm)
    estAperm <- permA(estA, trueA)
  else
    estAperm <- estA
  fr <- calFNRandFPR(estAperm, trueA)$error
  if (fr[1] == -1) {
    rates <- rep(-1, 3)
  } else {
    rates <- calArates(estAperm - trueA)
  }
  return(list(SPandSN = round(speAndSen,4), FNRandFPR = round(fr,4),rates = round(rates,4)))
}

threshA <- function(A, mu, scale = FALSE) {
  # Threshold the estimated {@code A} based on the given {@code mu}. If {@code scale} is true,
  # then normalize each row of A such that the l-1 norm of each row is not larger than 1.
  scaledA <- A
  for (i in 1:nrow(A)) {
    colInd <- abs(A[i, ]) <= mu
    scaledA[i,colInd] = 0
    if (scale && sum(abs(scaledA[i, ])) > 1)
      scaledA[i, ] <- scaledA[i, ] / sum(abs(scaledA[i, ]))
  }
  return(scaledA)
}

offSum <- function(M, N, weights) {
  # Calculate the sum of squares of the upper off-diagonal elements of two matrices
  # require: M and N have the same dimensions
  tmp <- (M-N) / weights
  tmp <- t(t(tmp) / weights)
  return(sum((tmp[row(tmp) <= (col(tmp) - 1)])^2))
}


Mse <- function(mat_A, mat_B) {mean((mat_A - mat_B) ** 2)}

# Mse_scale <- function(mat_A, mat_B) {mean((mat_A - mat_B) ** 2 / mat_B ** 2)}

Adev <- function(mat_A, mat_B) {mean(abs(mat_A - mat_B))}

Label_Y <- function(vec) {
  label_vec = rep(1, length(vec))
  label_vec[vec <= 4 / 13] = 0
  label_vec[vec >= 10 / 13] = 2
  label_vec
}


Acc <- function(Y_label, res_pred) {
  pred_label <- Label_Y(res_pred)
  mean(Y_label == pred_label)
}


