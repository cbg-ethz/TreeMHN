require(Matrix)

##' initialize_Theta(n, N, trees, lambda_s)
##' A function that initializes the diagonal entries of a MHN based on the
##' relative frequencies of each event starting from root.
##' @param n Number of events in the MHNs
##' @param N Sample size
##' @param trees Mutation tree structures
##' @param lambda_s Rate of the sampling event
##' @return An n-by-n matrix representing the initialized MHN
initialize_Theta <- function(n, N, trees, lambda_s) {
  initial_Theta <- matrix(0,nrow = n, ncol = n)

  for (l in c(1:N)) {
    tree <- trees[[l]]
    in_tree_nodes <- tree$nodes[tree$in_tree]
    in_tree_pa <- tree$parents[tree$in_tree]
    for (k in c(2:sum(tree$in_tree))) {
      j <- in_tree_nodes[k]
      if (in_tree_pa[k] == 1) {
        initial_Theta[j,j] <- initial_Theta[j,j] + 1
      }
    }
  }

  epsilon <- runif(n, 1e-4, 1e-3) # a vector of random small numbers
  diag(initial_Theta) <- diag(initial_Theta) / N

  temp <- (diag(initial_Theta) + epsilon) / (1- diag(initial_Theta) + epsilon) * lambda_s

  diag(initial_Theta) <- log(temp)
  off_diagonals <- which(initial_Theta == 0)

  # also randomize the off diagonal entries with small numbers
  initial_Theta[off_diagonals] <- runif(length(off_diagonals), 1e-4, 1e-3) *
    sample(c(-1,1),length(off_diagonals),replace = TRUE)

  return(initial_Theta)
}

##' learn_MHN(tree_obj, gamma, lambda_s, Theta_init, M, iterations, to_mask, use_EM, verbose)
##' This function learns a Mutual Hazard Network from a set of mutation trees
##' in the format of a TreeMHN object
##' @param tree_obj A TreeMHN object
##' @param gamma Penalization parameter in the objective function (Default: 0.5)
##' @param lambda_s Sampling rate (Default: 1)
##' @param Theta_init Initial value of the MHN provided to the optimization procedure (Default: NULL)
##' @param M Number of Monte Carlo samples to be drawn (Default: 100)
##' @param iterations Number of iterations for the EM/MCEM algorithm (Default: 1000)
##' @param to_mask An integer vector of indices by column, which is used to mask the
##' off-diagonal entries of an MHN (Default: an empty vector)
##' @param use_EM A boolean value to determine whether the EM/MCEM algorithm is used (Default: FALSE)
##' @param verbose A boolean value to determine whether optimization steps are printed (Default: FALSE)
##' @return A Mutual Hazard Network Theta
##' @author Xiang Ge Luo
##' @export
learn_MHN <- function(tree_obj, gamma = 0.5, lambda_s = 1, Theta_init = NULL,
                      M = 100, iterations = 1000, to_mask = integer(0),
                      use_EM = FALSE, verbose = FALSE) {

  n <- tree_obj$n
  N <- tree_obj$N
  trees <- tree_obj$trees
  N_patients <- tree_obj$N_patients ## New
  weights <- tree_obj$weights ## New
  smallest_tree_size <- min(sapply(trees, function (tree) sum(tree$in_tree) - 1))

  # initialize Theta
  if (is.null(Theta_init)) {
    cat("Initializing Theta...\n")
    Theta <- initialize_Theta(n, N, trees, lambda_s)
  } else if (nrow(Theta_init) != n) {
    stop("The dimension of the provided MHN is different from n. Please check again...")
  } else {
    Theta <- Theta_init
  }

  round <- 0
  reltol <- Inf
  ll <- -1e10
  cat("Checking whether MCEM is needed...\n")
  MC_flags <- get_MC_flags(N, n, trees)

  if (any(MC_flags) || use_EM) { # EM/MCEM

    cat("Running hybrid EM/MCEM...\n")
    while ((round < iterations) && (reltol > 1e-6)) {

      # E-step
      timed_trees <- get_timed_trees(n, N, trees, Theta, lambda_s, M, MC_flags)

      # M-step
      optim_res <- optim(Theta, full_MHN_objective, full_MHN_grad, timed_trees,
                         gamma, n, N, lambda_s, to_mask, weights, N_patients, smallest_tree_size,
                         method = "L-BFGS-B",
                         lower = -10, upper = 10,
                         control = list(fnscale = -1,
                                        trace = 0,
                                        maxit = 500,
                                        factr = 1e10))
      Theta <- optim_res$par
      round <- round + 1
      reltol <- abs(optim_res$value - ll) / abs(ll)
      if (M < 1000) {
        M <- M + 20
      }
      ll <- optim_res$value
      if ((round %% 10 == 0) && verbose) {
        cat("EM iteration round: ", round, "\n")
        cat("log-likelihood: ", ll, "\n")
      }
    }

  } else { # MLE

    cat("Running MLE...\n")
    obj_grad_help <- obj_grad_helper(n, N, trees, Theta, lambda_s)
    optim_res <- optim(Theta, obs_MHN_objective, obs_MHN_grad, n, N, lambda_s,
                       trees, gamma, obj_grad_help, to_mask, weights, N_patients, smallest_tree_size,
                       method = "L-BFGS-B",
                       lower = -10, upper = 10,
                       control = list(fnscale = -1,
                                      trace = verbose,
                                      maxit = 500,
                                      factr = 1e10))

    Theta <- optim_res$par

  }

  # Set masked elements to zero
  if (length(to_mask) > 0) {
    Theta[to_mask] <- 0
  }

  return(Theta)

}

##' get_MC_flags(N, trees)
##' A function that determines whether Monte Carlo sampling is needed for each tree
##' @param N Sample size
##' @param n Number of events
##' @param trees Mutation tree structures
##' @return A boolean vector indicating whether Monte Carlo sampling is used for each tree
get_MC_flags <- function(N, n, trees) {

  MC_flags <- logical(N)
  for (i in c(1:N)) {
    tree <- trees[[i]]

    if ((n <= 10) & (sum(tree$in_tree) > (n+1))) {
      MC_flags[i] <- TRUE
    } else if ((n > 10) & (sum(tree$in_tree) > (n/2+1))) {
      MC_flags[i] <- TRUE
    } else {
      p <- build_poset(tree)
      comp_geno <- compatible_genotypes(p)
      if (nrow(comp_geno) < 500) { # trees with fewer than 500 subtrees
        MC_flags[i] <- FALSE
      } else {
        MC_flags[i] <- TRUE
      }
    }
  }
  return(MC_flags)

}

##' get_timed_trees(n, N, trees, Theta, lambda_s, M, MC_flags)
##' A function that computes the expected waiting times of all edges in a set of trees
##' given a mutual hazard network
##' @param n Number of events in the MHNs
##' @param N Sample size
##' @param trees Mutation tree structures
##' @param Theta An n-by-n matrix representing the mutual hazard network
##' @param lambda_s Rate of the sampling event
##' @param M Number of Monte Carlo copies
##' @param MC_flags A boolean vector indicating whether Monte Carlo sampling is used for each tree
##' @return Timed trees
get_timed_trees <- function(n, N, trees, Theta, lambda_s, M, MC_flags) {

  timed_trees <- list()

  # Loop through all trees
  for (i in c(1:N)) {

    tree <- trees[[i]]

    if (MC_flags[i]) {  # Monte Carlo E-step for large trees

      time_diffs <- tree_importance_sampling(Theta, tree$nodes, tree$children,
                                             tree$in_tree, n, M, lambda_s)

    } else { # Exact E-step for small trees

      p <- build_poset(tree)
      comp_geno <- compatible_genotypes(p)
      node_labels <- get_node_labels(p)
      time_diffs <- tree_E_step(Theta, n, lambda_s,
                                tree$nodes, tree$children, tree$in_tree,
                                comp_geno, node_labels)

    }

    tree$time_diffs <- time_diffs
    timed_trees <- append(timed_trees, list(tree))

  }

  return(timed_trees)

}

full_MHN_objective <- function(Theta, trees, gamma, n, N, lambda_s, to_mask, 
                               weights, N_patients, smallest_tree_size = 1) {

  score <- full_MHN_objective_(Theta, trees, gamma, n, N, lambda_s, to_mask, weights)
  exp_Theta <- exp(matrix(Theta, nrow = n))
  if (smallest_tree_size == 1) {
    lambdas <- diag(exp_Theta)
    score <- score - N_patients * log(1 - prob_empty_tree(lambdas, lambda_s))
  }
  if (smallest_tree_size == 2) {
    lambdas <- diag(exp_Theta)
    score <- score - N_patients * log(1 - prob_empty_tree(lambdas, lambda_s) - prob_one_tree(n, exp_Theta, lambda_s))
  }

  return(score)
}

full_MHN_grad <- function(Theta, trees, gamma, n, N, lambda_s, to_mask, 
                          weights, N_patients, smallest_tree_size = 1) {

  gd <- full_MHN_grad_(Theta, trees, gamma, n, N, lambda_s, to_mask, weights)
  exp_Theta <- exp(matrix(Theta, nrow = n))
  if (smallest_tree_size == 1) {
    lambdas <- diag(exp_Theta)
    p_empty <- prob_empty_tree(lambdas, lambda_s)
    diag(gd) <- diag(gd) - N_patients * p_empty^2 / lambda_s / (1 - p_empty + 1e-10) * lambdas
  }
  if (smallest_tree_size == 2) {
    lambdas <- diag(exp_Theta)
    p_empty <- prob_empty_tree(lambdas, lambda_s)
    p_one <- prob_one_tree(n, exp_Theta, lambda_s)
    gd <- gd - N_patients * (diag(p_empty^2 / lambda_s * lambdas) - dp_one(n, exp_Theta, lambda_s)) / (1 - p_empty - p_one + 1e-10)
  }

  return(gd)
}

obs_MHN_objective <- function(Theta, n, N, lambda_s, trees, gamma, obj_grad_help, to_mask, 
                              weights, N_patients, smallest_tree_size = 1) {

  score <- obs_MHN_objective_(Theta, n, N, lambda_s, trees, gamma, obj_grad_help, to_mask, weights)
  exp_Theta <- exp(matrix(Theta, nrow = n))
  if (smallest_tree_size == 1) {
    lambdas <- diag(exp_Theta)
    score <- score - N_patients * log(1 - prob_empty_tree(lambdas, lambda_s))
  }
  if (smallest_tree_size == 2) {
    lambdas <- diag(exp_Theta)
    score <- score - N_patients * log(1 - prob_empty_tree(lambdas, lambda_s) - prob_one_tree(n, exp_Theta, lambda_s))
  }

  return(score)
}

obs_MHN_grad <- function(Theta, n, N, lambda_s, trees, gamma, obj_grad_help, to_mask, 
                         weights, N_patients, smallest_tree_size = 1) {

  gd <- obs_MHN_grad_(Theta, n, N, lambda_s, trees, gamma, obj_grad_help, to_mask, weights)
  exp_Theta <- exp(matrix(Theta, nrow = n))
  if (smallest_tree_size == 1) {
    lambdas <- diag(exp_Theta)
    p_empty <- prob_empty_tree(lambdas, lambda_s)
    diag(gd) <- diag(gd) - N_patients * p_empty^2 / lambda_s / (1 - p_empty + 1e-10) * lambdas
  }
  if (smallest_tree_size == 2) {
    lambdas <- diag(exp_Theta)
    p_empty <- prob_empty_tree(lambdas, lambda_s)
    p_one <- prob_one_tree(n, exp_Theta, lambda_s)
    gd <- gd - N_patients * (diag(p_empty^2 / lambda_s * lambdas) - dp_one(n, exp_Theta, lambda_s)) / (1 - p_empty - p_one + 1e-10)
  }

  return(gd)
}


prob_empty_tree <- function(lambdas, lambda_s) {
  return(lambda_s / (lambda_s + sum(lambdas)))
}

prob_one_tree <- function(n, exp_Theta, lambda_s) {

  p <- 0
  sum_lambdas <- sum(diag(exp_Theta))
  for (i in c(1:n)) {
    temp_denom <- sum(sapply(setdiff(c(1:n),i), function (j) exp_Theta[j,j] * (1 + exp_Theta[j,i])))
    p <- p + exp_Theta[i,i] / (lambda_s + sum_lambdas) * lambda_s / (lambda_s + temp_denom)
  }
  return(p)

}

dp_one <- function(n, exp_Theta, lambda_s) {

  sum_lambdas <- sum(diag(exp_Theta))
  sum_all <- lambda_s + sum_lambdas
  dp <- matrix(0, nrow = n, ncol = n)

  for (i in c(1:n)) {
    temp_denom_i <- sum(sapply(setdiff(c(1:n),i), function (j) exp_Theta[j,j] * (1 + exp_Theta[j,i])))
    dp[-i,i] <- - exp_Theta[i,i] / sum_all * lambda_s * diag(exp_Theta)[-i] / (lambda_s + temp_denom_i)^2

    dp[i,i] <- lambda_s / (lambda_s + temp_denom_i) * (sum_all^(-1) - exp_Theta[i,i] * sum_all^(-2))
    for (k in setdiff(c(1:n),i)) {
      temp_denom_k <- sum(sapply(setdiff(c(1:n),k), function (j) exp_Theta[j,j] * (1 + exp_Theta[j,k])))
      dp[i,i] <- dp[i,i] - exp_Theta[k,k] * lambda_s *
        (sum_all^(-2) * (lambda_s + temp_denom_k)^(-1) + sum_all^(-1) * (1 + exp_Theta[i,k]) * (lambda_s + temp_denom_k)^(-2))
    }
  }

  return(dp)
}

