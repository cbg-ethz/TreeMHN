require(gtools)

##' compare_Theta(true_Theta, pred_Theta, q)
##' This function computes the differences between two Mutual Hazard Networks
##' @param true_Theta A ground truth MHN represented by a square matrix
##' @param pred_Theta An estimated MHN represented by a square matrix
##' @param q A threshold to zero out very small entries in the estimated MHN (Default: 1e-2)
##' @return A vector of performance metrics:
##' \itemize{
##' \item SHD: Structural Hamming Distance between the two matrices;
##' \item TP: True positives in the estimated MHN;
##' \item FP: False positives in the estimated MHN;
##' \item TN: True negatives in the estimated MHN;
##' \item FN: False negatives in the estimated MHN;
##' \item Precision: True positives divided by the total number of edges in the estimated MHN;
##' \item TPR (Recall): True positives divided by the total number of edges in the true MHN;
##' \item FPR_N: False positives divided by the total number of non-edges in the true MHN;
##' \item FPR_P: False positives divided by the total number of edges in the true MHN;
##' \item MSE: Mean squared error between the two matrices.
##' }
##' @author Xiang Ge Luo
##' @export
compare_Theta <- function(true_Theta, pred_Theta, q = 1e-2) {

  if ((nrow(true_Theta) != nrow(pred_Theta)) ||
      (ncol(true_Theta) != ncol(pred_Theta)) ||
      (nrow(true_Theta) != ncol(true_Theta)) ||
      (nrow(pred_Theta) != ncol(pred_Theta))) {
    stop("The dimensions of the two MHN matrices are different!")
  }

  MSE <- mean((true_Theta - pred_Theta)^2)
  n <- length(diag(true_Theta))
  diag(true_Theta) <- 0
  diag(pred_Theta) <- 0
  pred_Theta[pred_Theta > q] <- 1
  pred_Theta[pred_Theta < -q] <- -1
  pred_Theta[(pred_Theta <= q) & (pred_Theta >= -q)] <- 0
  true_Theta[true_Theta > 0] <- 1
  true_Theta[true_Theta < 0] <- -1

  # Number of edges in the estimated Theta
  pred_P <- sum(pred_Theta != 0)

  # Number of edges in the true Theta
  true_P <- sum(true_Theta != 0)

  # Number of non-edges in the true Theta
  true_N <- sum(true_Theta == 0) - n

  # TP, FP, TN, FN, SHD
  TP <- sum((pred_Theta != 0) * (pred_Theta == true_Theta))
  FP <- pred_P - TP
  FN <- sum((pred_Theta == 0) * (true_Theta != 0))
  TN <- sum((pred_Theta == 0) * (true_Theta == 0)) - n
  SHD <- FP + FN

  # Precision
  if ((TP + FP) == 0) {
    Precision <- 0
  } else {
    Precision <- TP / (TP + FP)
  }

  # TPR, FPR_P, FPR_N
  if (true_P == 0) { # true graph is empty
    if (FP >= 0) {
      TPR <- 0
      FPR_P <- 1
    } else {
      TPR <- 1
      FPR_P <- 0
    }
  } else { # true graph is non-empty
    TPR <- TP / true_P
    FPR_P <- FP / true_P
  }

  if (true_N == 0) { # true graph is full
    FPR_N <- 0
  } else { # true graph is not full
    FPR_N <- FP / true_N
  }

  compTheta <- c(SHD,TP,FP,TN,FN,Precision,TPR,FPR_N,FPR_P,MSE)
  names(compTheta) <- c("SHD","TP","FP","TN","FN","Precision","TPR","FPR_N","FPR_P","MSE")
  return(round(compTheta,2))

}

##' Theta_to_pathways(Theta, n_order, prob_only)
##' This function computes the pathway probabilities given a Mutual Hazard Network
##' @param Theta An MHN represented by a square matrix
##' @param n_order Length of the pathways (Default: 4)
##' @param prob_only A Boolean value that determines whether to output only the pathway probabilities
##' or the data frame containing also the pathways (Default: TRUE)
##' @return Pathway probabilities
##' @author Xiang Ge Luo
##' @export
Theta_to_pathways <- function(Theta, n_order = 4, prob_only = TRUE) {

  n <- nrow(Theta)
  if (n < n_order) {
    stop("The number of mutations is smaller than the order. Please check again...")
  }
  pathways <- gtools::permutations(n, n_order)
  temp <- matrix(0, nrow = nrow(pathways), ncol = n_order)
  temp[,1] <- sapply(c(1:n), function (i) rep(exp(Theta[i,i]),
                                              prod(sapply(c(0:(n_order - 2)),
                                                          function (j) (n - 1 - j)))))
  temp[,1] <- temp[,1] / sum(exp(diag(Theta)))
  exp_time <- matrix(0, nrow = nrow(pathways), ncol = n_order)
  exp_time[,1] <- 1 / sum(exp(diag(Theta)))

  if (n_order > 1) {

    # compute lambdas from Theta
    for (i in c(1:nrow(pathways))) {
      for (j in c(2:n_order)) {
        temp[i,j] <- get_lambda(pathways[i,c(1:j)], Theta)
      }
    }

    # compute pathway probabilities
    for (j in c(2:n_order)) {
      nr_pa <- prod(sapply(c(0:(j - 2)), function (i) (n - i)))
      prod_factor <- prod(sapply(c(0:(n_order - j)), function (i) (n - j - i + 1)))
      for (k in c(1:nr_pa)) {
        idx <- c((1 + (k - 1) * prod_factor):(k * prod_factor))
        if (length(idx) > 1) {
          to_sum <- !duplicated(pathways[idx, c(1:j)])
          exp_time[idx, j] <- 1 / sum(temp[idx, j] * to_sum)
          temp[idx, j] <- temp[idx, j] / sum(temp[idx, j] * to_sum)
        } else {
          exp_time[idx, j] <- 1 / sum(temp[idx, j])
          temp[idx, j] <- temp[idx, j] / sum(temp[idx, j])
        }
      }
    }
  }

  prob <- apply(temp,1,prod)
  if (prob_only) {
    return(prob)
  } else {
    df <- data.frame(pathways)
    df$prob <- prob
    df <- cbind(df, exp_time)
    return(df[order(df$prob,decreasing = TRUE),])
  }

}

##' KL_divergence(p, q)
##' This function computes the KL divergence between two probability distributions
##' @param p A probability distribution
##' @param q Another probability distribution
##' @return KL(p || q)
##' @export
KL_divergence <- function(p, q) {
  return(as.numeric(p %*% log(p) - p %*% log(q)))
}

##' trees_to_revolver_W(n, tree_df)
##' This function computes the row-normalized W matrix of the REVOLVER algorithm
##' @param n Number of mutational events
##' @param tree_df A data frame with the following columns:
##' \itemize{
##' \item Tree_ID: IDs of mutation trees, unique for each patient
##' \item Node_ID: IDs of each node in the tree, including the root node (with ID "1"), unique for each node
##' \item Mutation_ID: IDs of each mutational event, the root node has a mutation ID of "0",
##' other mutation IDs can be duplicated in the tree to allow for parallel mutations
##' \item Parent_ID: IDs of the parent node ID. The root node has itself as parent (ID "1").
##' }
##' @return Row-normalized W matrix of the REVOLVER algorithm
##' @author Xiang Ge Luo
##' @export
trees_to_revolver_W <- function(n, tree_df) {

  W <- matrix(0, nrow = n + 1, ncol = n + 1) # This matrix contains also the wild type (GL)

  # Go through all trees and count the edge frequencies
  for (i in c(1:nrow(tree_df))) {
    idx <- which((tree_df$Tree_ID == tree_df$Tree_ID[i]) & (tree_df$Node_ID == tree_df$Parent_ID[i]))
    pa <- tree_df$Mutation_ID[idx]
    W[pa + 1, tree_df$Mutation_ID[i] + 1] <- W[pa + 1, tree_df$Mutation_ID[i] + 1] + 1
  }

  W <- W + 1 # add pseudocounts
  diag(W) <- 0 # no mutation can be the descendant of itself
  W <- W / rowSums(W) # normalize by rows
  return(W)

}

##' get_revolver_pathways(n, tree_df, n_order)
##' This function computes the pathway probabilities inferred from the
##' row-normalized W matrix of the REVOLVER algorithm
##' @param n Number of mutational events
##' @param tree_df A data frame with the following columns:
##' \itemize{
##' \item Tree_ID: IDs of mutation trees, unique for each patient
##' \item Node_ID: IDs of each node in the tree, including the root node (with ID "1"), unique for each node
##' \item Mutation_ID: IDs of each mutational event, the root node has a mutation ID of "0",
##' other mutation IDs can be duplicated in the tree to allow for parallel mutations
##' \item Parent_ID: IDs of the parent node ID. The root node has itself as parent (ID "1").
##' }
##' @param n_order Length of the pathways (Default: 4)
##' @return Pathway probabilities inferred from the row-normalized W matrix of the REVOLVER algorithm
##' @author Xiang Ge Luo
##' @export
get_revolver_pathways <- function(n, tree_df, n_order = 4) {

  W <- trees_to_revolver_W(n, tree_df)
  if (n < n_order) {
    stop("The number of mutations is smaller than the order. Please check again...")
  }
  pathways <- cbind(0, gtools::permutations(n, n_order)) + 1
  probs <- rep(1, nrow(pathways))

  for (i in c(1:nrow(pathways))) {
    for (j in c(2:(n_order + 1))) {
      probs[i] <- probs[i] * W[pathways[i,j-1], pathways[i,j]]
    }
  }

  probs <- probs / sum(probs) # normalize
  return(probs)

}

##' get_hintra_pathways(n, tree_df, n_order)
##' This function computes the pathway probabilities inferred from the
##' row-normalized beta matrix of the HINTRA algorithm
##' @param n Number of mutational events
##' @param tree_df A data frame with the following columns:
##' \itemize{
##' \item Tree_ID: IDs of mutation trees, unique for each patient
##' \item Node_ID: IDs of each node in the tree, including the root node (with ID "1"), unique for each node
##' \item Mutation_ID: IDs of each mutational event, the root node has a mutation ID of "0",
##' other mutation IDs can be duplicated in the tree to allow for parallel mutations
##' \item Parent_ID: IDs of the parent node ID. The root node has itself as parent (ID "1").
##' }
##' @param n_order Length of the pathways (Default: 4)
##' @return Pathway probabilities inferred from the row-normalized beta matrix of the HINTRA algorithm
##' @author Xiang Ge Luo
##' @export
get_hintra_pathways <- function(n, tree_df, n_order = 4) {

  hintra_B_vec <- c("0")
  to_mask <- list()
  for (r in c(1:(n_order - 1))) {
    pathways <- gtools::combinations(n,r)
    temp <- apply(pathways, 1, function(x) paste0(x,collapse = "_"))
    hintra_B_vec <- c(hintra_B_vec,temp)
    to_mask <- append(to_mask, apply(pathways,1,list))
  }

  B <- matrix(0, nrow = length(hintra_B_vec), ncol = n)

  N <- length(unique(tree_df$Tree_ID))
  for (i in c(1:N)) {
    tree <- tree_df[tree_df$Tree_ID == i,]
    for (j in c(2:nrow(tree))) {
      node <- tree$Mutation_ID[j]
      set_P <- c()
      current_pos <- j
      repeat{
        pa_pos <- which(tree$Node_ID == tree$Parent_ID[current_pos])
        pa <- tree$Mutation_ID[pa_pos]
        if (pa == 0) {
          break
        } else {
          set_P <- c(pa, set_P)
          current_pos <- pa_pos
        }
      }
      if (is.null(set_P)) {
        set_P <- 0
      }
      P_idx <- which(hintra_B_vec == paste0(sort(set_P),collapse = "_"))
      B[P_idx,node] <- B[P_idx,node] + 1
    }
  }
  B <- B + 1

  for (i in c(1:length(to_mask))) {
    B[i+1, simplify2array(to_mask[[i]])] <- 0
  }

  B <- B / rowSums(B)

  pathways <- gtools::permutations(n, n_order)
  probs <- B[1, pathways[,1]]

  for (i in c(1:nrow(pathways))) {
    for (j in c(2:n_order)) {
      row_idx <- which(hintra_B_vec == paste0(sort(pathways[i,1:(j-1)]),collapse = "_"))
      probs[i] <- probs[i] * B[row_idx, pathways[i,j]]
    }
  }

  return(probs)

}





