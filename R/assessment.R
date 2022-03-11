
##' @name compare_Theta
##' @title Compute the differences between two Mutual Hazard Networks
##' @description This function computes the differences between two Mutual Hazard Networks.
##' @param true_Theta A ground truth MHN represented by a square matrix.
##' @param pred_Theta An estimated MHN represented by a square matrix.
##' @param q A threshold to zero out very small entries in the estimated MHN (Default: 1e-2).
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

##' @name Theta_to_pathways
##' @title Compute the pathway probabilities given a Mutual Hazard Network
##' @description This function computes the pathway probabilities given a Mutual Hazard Network.
##' @param Theta An MHN represented by a square matrix.
##' @param n_order Length of the pathways (Default: 4).
##' @param prob_only A Boolean value that determines whether to output only the pathway probabilities
##' or the data frame containing also the pathways (Default: TRUE).
##' @return Pathway probabilities
##' @author Xiang Ge Luo
##' @import gtools
##' @export
Theta_to_pathways <- function(Theta, n_order = 4, prob_only = TRUE) {

  n <- nrow(Theta)
  if (n < n_order) {
    stop("The number of mutations is smaller than the order. Please check again...")
  }
  pathways <- permutations(n, n_order)
  temp <- matrix(0, nrow = nrow(pathways), ncol = n_order)
  if (n_order == 1) {
    temp[,1] <- sapply(c(1:n), function (i) exp(Theta[i,i]))
  } else {
    temp[,1] <- sapply(c(1:n), function (i) rep(exp(Theta[i,i]),
                                                prod(sapply(c(0:(n_order - 2)),
                                                            function (j) (n - 1 - j)))))
  }
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

##' @name KL_divergence
##' @title Compute the KL divergence between two probability distributions
##' @description This function computes the KL divergence between two probability distributions.
##' @param p A probability distribution.
##' @param q Another probability distribution.
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
##' @return Row-normalized W matrix of the REVOLVER algorithm.
##' @author Xiang Ge Luo
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

##' @name get_revolver_pathways
##' @title Compute REVOLVER pathway probabilities
##' @description This function computes the pathway probabilities inferred from the
##' row-normalized W matrix of the REVOLVER algorithm
##' @param n Number of mutational events.
##' @param tree_df A data frame with the following columns:
##' \itemize{
##' \item Tree_ID: IDs of mutation trees, unique for each patient
##' \item Node_ID: IDs of each node in the tree, including the root node (with ID "1"), unique for each node
##' \item Mutation_ID: IDs of each mutational event, the root node has a mutation ID of "0",
##' other mutation IDs can be duplicated in the tree to allow for parallel mutations
##' \item Parent_ID: IDs of the parent node ID. The root node has itself as parent (ID "1").
##' }
##' @param n_order Length of the pathways (Default: 4).
##' @return A list containing the matrix W and the pathway probabilities 
##' inferred from the row-normalized W matrix of the REVOLVER algorithm.
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
  
  res <- list(W = W, probs = probs)
  return(res)

}

##' @name get_hintra_pathways
##' @title Compute HINTRA pathway probabilities
##' @description This function computes the pathway probabilities inferred from the
##' row-normalized beta matrix of the HINTRA algorithm.
##' @param n Number of mutational events.
##' @param tree_df A data frame with the following columns:
##' \itemize{
##' \item Tree_ID: IDs of mutation trees, unique for each patient
##' \item Node_ID: IDs of each node in the tree, including the root node (with ID "1"), unique for each node
##' \item Mutation_ID: IDs of each mutational event, the root node has a mutation ID of "0",
##' other mutation IDs can be duplicated in the tree to allow for parallel mutations
##' \item Parent_ID: IDs of the parent node ID. The root node has itself as parent (ID "1").
##' }
##' @param n_order Length of the pathways (Default: 4).
##' @return A list containing the matrix beta and the pathway probabilities 
##' inferred from the row-normalized beta matrix of the HINTRA algorithm.
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

  res <- list(B = B, probs = probs)
  return(res)

}

get_children <- function(n, pathway) {
  
  nr_ch <- n - length(pathway)
  to_add <- setdiff(c(1:n), pathway) 
  pathways <- vector("list", nr_ch)
  if (nr_ch > 0) {
    for (i in c(1:nr_ch)) {
      pathways[[i]] <- c(pathway, to_add[i])
    }
  }
  return(pathways)
  
}

Theta_to_pathways_w_sampling <- function(Theta, top_M = 10, lambda_s = 1) {
  
  n <- nrow(Theta)
  pathways <- vector("list", top_M)
  probs <- rep(0, top_M)
  
  current_pathways <- as.list(c(1:n))
  while (length(current_pathways) > 0) {
    
    next_pathways <- list()
    for (p in current_pathways) {
      p_prob <- 1
      for (i in c(1:length(p))) {
        pp <- p[c(1:i)]
        num <- TreeMHN:::get_lambda(pp, Theta)
        denom_set <- get_children(n, pp[-length(pp)])
        denom <- lambda_s + sum(sapply(denom_set, function (l) TreeMHN:::get_lambda(l, Theta)))
        p_prob <- p_prob * num / denom
      }
      # times the probability of the pathway stopping at the sampling event
      p_ch <- get_children(n, p)
      if (length(p_ch) > 0) {
        p_prob <- p_prob * lambda_s / (lambda_s + sum(sapply(p_ch, function (l) TreeMHN:::get_lambda(l, Theta)))) 
      }
      if (p_prob > min(probs)) {
        to_replace <- which.min(probs)
        pathways[[to_replace]] <- p
        probs[to_replace] <- p_prob
        next_pathways <- append(next_pathways, p_ch)
      }
    }
    current_pathways <- next_pathways
    
  }
  
  probs <- probs / (1 - lambda_s / (lambda_s + sum(exp(diag(Theta)))))
  pathways <- sapply(pathways, function (x) paste(x, collapse = "_"))
  df <- data.frame(pathways, probs) %>% arrange(desc(probs))
  return(df)
  
}

pathway_helper <- function(tree, index = 1, pathway = character(0)) {
  
  ch_set <- tree$children[[index]]
  ch_set <- ch_set[tree$in_tree[ch_set]]
  
  if (index != 1) {
    if (length(pathway) == 0) {
      pathway <- as.character(tree$nodes[index])
    } else {
      pathway <- paste(pathway, tree$nodes[index], sep = "_")
    }
  }
  
  if (length(ch_set) == 0) {
    return(pathway)
  } else {
    return(sapply(ch_set, function(ch) pathway_helper(tree, ch, pathway)))
  }
  
}

##' @import dplyr
get_observed_pathways <- function(tree_obj) {
  
  pathways <- c()
  
  trees <- tree_obj$trees
  
  for (tree in trees) {
    pathways <- c(pathways, unlist(pathway_helper(tree)))
  }
  
  df <- data.frame(pathways) %>% 
    group_by_all() %>%
    summarise(count = n()) %>%
    mutate(probs = count / sum(count)) %>%
    arrange(desc(count))
  
  return(df)
  
}

##' @import dplyr
get_next_mutations <- function(tree, Theta, mutations = NULL) {
  
  n <- nrow(Theta)
  
  if (is.null(mutations)) {
    mutations <- as.character(seq(1,n))
  } else {
    if (length(mutations) != n) {
      stop("The number of mutations doesn't match matrix dimension. Please check again...")
    } else if (length(unique(mutations)) != n) {
      stop("Mutation names must be unique. Please check again...")
    }
  }
  
  next_mutations <- which(tree$in_tree == FALSE)
  nr_next_mutations <- length(next_mutations)
  next_lambdas <- rep(0, nr_next_mutations)
  pathways <- rep("", nr_next_mutations)
  last_mutations <- rep("", nr_next_mutations)
  
  for (i in c(1:nr_next_mutations)) {
    
    pos <- next_mutations[i]
    node <- get_pathway(tree$nodes, pos, tree$parents)
    next_lambdas[i] <- get_lambda(node, Theta)
    pathways[i] <- paste(c("Root", mutations[node]), collapse = "->")
    last_mutations[i] <- mutations[tree$nodes[pos]]
    
  }
  
  df <- data.frame(pathways, last_mutations, next_lambdas) %>%
    mutate(probs = next_lambdas / sum(next_lambdas)) %>%
    select(!next_lambdas) %>%
    arrange(desc(probs)) %>%
    mutate(rank = seq(1, nr_next_mutations))
  
  return(df)
  
}


