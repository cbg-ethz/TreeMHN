require(Matrix)

##' initialize_Theta(n, N, trees, lambda_s)
##' A function that initializes the diagonal entries of a MHN based on the
##' relative frequencies of each event starting from root.
##' @param n Number of events in the MHNs
##' @param N Sample size
##' @param trees Mutation tree structures
##' @param lambda_s Rate of the sampling event
##' @return An n-by-n matrix representing the initialized MHN
##' @importFrom stats runif
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

##' @name learn_MHN
##' @title Learn a Mutual Hazard Network from a set of mutation trees
##' @description This function learns a Mutual Hazard Network from a set of mutation trees
##' in the format of a TreeMHN object.
##' @param tree_obj A TreeMHN object
##' @param gamma Penalization parameter in the objective function (Default: 0.5).
##' @param lambda_s Sampling rate (Default: 1).
##' @param Theta_init Initial value of the MHN provided to the optimization procedure (Default: NULL).
##' @param M Number of Monte Carlo samples to be drawn (Default: 100).
##' @param iterations Number of iterations for the EM/MCEM algorithm (Default: 500).
##' @param to_mask An integer vector of indices by column, which is used to mask the
##' off-diagonal entries of an MHN (Default: an empty vector).
##' @param use_EM A boolean value to determine whether the EM/MCEM algorithm is used (Default: FALSE).
##' @param verbose A boolean value to determine whether optimization steps are printed (Default: FALSE).
##' @param MC_threshold A threshold on the maximum number of subtrees of a given tree,
##' above which Monte Carlo sampling will be used (Default: 500). 
##' @param increment_M The step size to increment the number of Monte Carlo samples (Default: 0).
##' @param increment_M_bound The upper bound on the number of Monte Carlo samples (Default: 500).
##' @return A Mutual Hazard Network Theta
##' @author Xiang Ge Luo
##' @importFrom stats optim
##' @export
learn_MHN <- function(tree_obj, gamma = 0.5, lambda_s = 1, Theta_init = NULL,
                      M = 100, iterations = 500, to_mask = integer(0),
                      use_EM = FALSE, verbose = FALSE, MC_threshold = 500,
                      increment_M = 0, increment_M_bound = 500) {

  n <- tree_obj$n
  N <- tree_obj$N
  trees <- tree_obj$trees
  N_patients <- tree_obj$N_patients ## New
  weights <- tree_obj$weights ## New
  smallest_tree_size <- min(sapply(trees, function (tree) sum(tree$in_tree) - 1))

  # initialize Theta
  if (is.null(Theta_init)) {
    if (verbose) {
      cat("Initializing Theta...\n")
    }
    Theta <- initialize_Theta(n, N, trees, lambda_s)
  } else if (nrow(Theta_init) != n) {
    stop("The dimension of the provided MHN is different from n. Please check again...")
  } else {
    Theta <- Theta_init
  }

  round <- 1
  reltol <- Inf
  ll <- -1e10
  if (verbose) {
    cat("Checking whether MCEM is needed...\n")
  }
  MC_flags <- get_MC_flags(N, n, trees, MC_threshold)

  if (any(MC_flags) || use_EM) { # EM/MCEM
    
    # Re-order trees such that trees with exact inference go first
    nr_exact <- N - sum(MC_flags)
    timed_trees <- list()
    comp_geno_vec <- list()
    node_labels_vec <- list()
    for (i in order(MC_flags)) {
      tree <- trees[[i]]
      tree$time_diffs <- numeric(0) # Initialize the time difference vector
      if (!MC_flags[i]) { # if exact inference is needed
        p <- build_poset(tree)
        comp_geno_vec <- append(comp_geno_vec, list(compatible_genotypes(p)))
        node_labels_vec <- append(node_labels_vec, list(get_node_labels(p)))
      }
      timed_trees <- append(timed_trees, list(tree))
    }
    rm(tree)
    rm(trees)
    rm(MC_flags)
    
    if (verbose) {
      
      cat("Running hybrid EM/MCEM...\n")
      
      while ((round < iterations) && (reltol > 1e-6)) {
        cat("EM iteration round: ", round, "\n")
        
        # E-step
        cat("E-step...\n")
        start_time <- Sys.time()
        update_timed_trees(n, N, timed_trees, Theta, lambda_s, M,
                           comp_geno_vec, node_labels_vec, nr_exact)
        print(Sys.time() - start_time)
        
        # M-step
        cat("M-step...\n")
        start_time <- Sys.time()
        optim_res <- optim(Theta, full_MHN_objective, full_MHN_grad, timed_trees,
                           gamma, n, N, lambda_s, to_mask, weights, N_patients, smallest_tree_size,
                           method = "L-BFGS-B",
                           lower = -10, upper = 10,
                           control = list(fnscale = -1,
                                          trace = 0,
                                          maxit = 100,
                                          factr = 1e10))
        Theta <- optim_res$par
        round <- round + 1
        reltol <- abs(optim_res$value - ll) / abs(ll)
        if (M < increment_M_bound) {
          M <- M + increment_M
        }
        ll <- optim_res$value
        print(Sys.time() - start_time)
        cat("Complete-data log-likelihood: ", ll, "\n")
      }
      
    } else {
      
      while ((round < iterations) && (reltol > 1e-6)) {
        
        # E-step
        update_timed_trees(n, N, timed_trees, Theta, lambda_s, M,
                           comp_geno_vec, node_labels_vec, nr_exact)
        
        # M-step
        optim_res <- optim(Theta, full_MHN_objective, full_MHN_grad, timed_trees,
                           gamma, n, N, lambda_s, to_mask, weights, N_patients, smallest_tree_size,
                           method = "L-BFGS-B",
                           lower = -10, upper = 10,
                           control = list(fnscale = -1,
                                          trace = 0,
                                          maxit = 100,
                                          factr = 1e10))
        Theta <- optim_res$par
        round <- round + 1
        reltol <- abs(optim_res$value - ll) / abs(ll)
        if (M < increment_M_bound) {
          M <- M + increment_M
        }
        ll <- optim_res$value
      }
      
    }


  } else { # MLE
    
    if (verbose) {
      cat("Running MLE...\n")
    }

    obj_grad_help <- obj_grad_helper(n, N, trees, Theta, lambda_s)
    tr_mat_vec <- obj_grad_help$tr_mat_vec
    comp_geno_vec <- obj_grad_help$comp_geno_vec
    node_labels_vec <- obj_grad_help$node_labels_vec
    log_prob_vec <- obj_grad_help$log_prob_vec
    rm(obj_grad_help)
    optim_res <- optim(Theta, obs_MHN_objective, obs_MHN_grad, n, N, lambda_s, trees, 
                       gamma, tr_mat_vec, log_prob_vec, comp_geno_vec, node_labels_vec,
                       to_mask, weights, N_patients, smallest_tree_size,
                       method = "L-BFGS-B",
                       lower = -10, upper = 10,
                       control = list(fnscale = -1,
                                      trace = verbose,
                                      maxit = iterations,
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
##' @param MC_threshold A threshold on the maximum number of subtrees of a given tree,
##' above which Monte Carlo sampling will be used (Default: 500). 
##' @return A boolean vector indicating whether Monte Carlo sampling is used for each tree
get_MC_flags <- function(N, n, trees, MC_threshold) {

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
      if (nrow(comp_geno) < MC_threshold) {
        MC_flags[i] <- FALSE
      } else {
        MC_flags[i] <- TRUE
      }
    }
  }
  return(MC_flags)

}

