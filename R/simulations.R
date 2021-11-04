
##' random_Theta(n, sparsity, exclusive_ratio)
##' This function randomly generates an MHN with log-transformed parameters.
##' @param n Number of mutational events
##' @param sparsity Percentage of off diagonal elements with a value of zero in the MHN. (Default: 0.5)
##' @param exclusive_ratio Percentage of non-zero elements that are negative,
##' meaning that one mutation is inhibiting another mutation (Default: 0.9)
##' @return An n-by-n matrix representing the MHN
##' @author Rudolf Schill, Xiang Ge Luo
##' @references Schill et al. (2020). Modelling cancer progression using Mutual Hazard Networks. Bioinformatics, 36(1), 241–249.
##' @export
random_Theta <- function(n, sparsity = 0.5, exclusive_ratio = 0.5, sim_setup = 1) {
  Theta  <- matrix(0, nrow=n, ncol=n)

  diag(Theta)  <- NA
  nonZeroSize <- (n^2 - n) * (1 - sparsity)
  nonZeros <- sample(which(!is.na(Theta)), size=nonZeroSize)
  exclusive <- sample(nonZeros, floor(nonZeroSize * exclusive_ratio))
  cooccur <- setdiff(nonZeros, exclusive)

  # Parameters (can change)
  switch(as.character(sim_setup),
         "1" = {
           # Setup 1: hyper parameters estimated from real data
           diag(Theta) <- sample(c(runif(ceiling(n/2), -3, -1), runif(floor(n/2), -6, -3)))
           Theta[exclusive] <- - rgamma(length(exclusive), 4, 2.5)
           Theta[cooccur] <- rgamma(length(cooccur), 4, 2.5)
         },
         "2" = {
           # Setup 2:
           diag(Theta) <- sort(rnorm(100,-2,1))[c(1:n)]
           temp <- rnorm(length(nonZeros),0,3)
           Theta[nonZeros] <- sapply(temp, function (x) ifelse(abs(x) > 0.1, x, runif(1,0.1,1) * sign(x)))
         })

  return(round(Theta, 2))
}

##' generate_trees(n, N, lambda_s, Theta, sparsity, exclusive_ratio, non_empty)
##' This function generate a set of trees based on an MHN, which is either
##' provided by the user or generated randomly within the function.
##' @param n Number of mutational events (Default: 10)
##' @param N Number of trees (Default: 100)
##' @param lambda_s The rate of the sampling event. (Default: 1)
##' @param Theta A mutual hazard network. (Default: NULL for the function to generate)
##' @param sparsity Percentage of off diagonal elements with a value of zero in the MHN. (Default: 0.5)
##' @param exclusive_ratio Percentage of non-zero elements that are negative,
##' meaning that one mutation is inhibiting another mutation. (Default: 0.9)
##' @param non_empty A flag for the function to generate non-empty trees only. (Default: TRUE)
##' @param mutations A list of mutation names, which must be unique values. (Default: NULL)
##' @param perturb A flag for the function to perturb the generated tree structures (Default: FALSE)
##' @param epsilon Noise level for perturbing the trees
##' @return A treeMHN object with the true MHN for generating the trees
##' @author Xiang Ge Luo
##' @export
generate_trees <- function(n = 10, N = 100, lambda_s = 1, Theta = NULL,
                           sparsity = 0.5, exclusive_ratio = 0.5, non_empty = TRUE,
                           mutations = NULL, perturb = FALSE, epsilon = 0.05, sim_setup = 1) {

  if (is.null(Theta)) {
    Theta <- random_Theta(n, sparsity, exclusive_ratio, sim_setup)
  } else if (nrow(Theta) != n) {
    stop("The dimension of the provided MHN is different from n. Please check again...")
  }

  if (is.null(mutations)) {
    mutations <- as.character(seq(1,n))
  } else {
    if (length(mutations) != n) {
      stop("The list of mutations has length different from n. Please check again...")
    } else if (length(unique(mutations)) != n) {
      stop("Mutation names must be unique. Please check again...")
    }
  }

  trees <- list()
  all_mutation_check <- FALSE

  if (non_empty) { # empty trees are not allowed
    while ((length(trees) < N) || !all_mutation_check) {
      tree <- generate_tree_MHN(n, Theta, lambda_s)
      if (length(tree) != 0) {
        # check if the tree is empty
        if (sum(tree$in_tree) <= 1) {
          next
        }
        if (length(trees) >= N) {
          trees[[sample(seq(1,N),1)]] <- NULL
        }
        trees <- append(trees, list(tree))
        all_mutation_check <- check_mutations(n, trees)
      }
    }
  } else { # empty trees are allowed
    while ((length(trees) < N) || !all_mutation_check) {
      tree <- generate_tree_MHN(n, Theta, lambda_s)
      if (length(tree) != 0) {
        if (length(trees) >= N) {
          trees[[sample(seq(1,N),1)]] <- NULL
        }
        trees <- append(trees, list(tree))
        all_mutation_check <- check_mutations(n, trees)
      }
    }
  }

  # Add tree IDs
  for (i in c(1:N)) {
    trees[[i]]$tree_ID <- i
  }

  tree_obj <- list("n" = n,
                   "N" = N,
                   "mutations" = mutations,
                   "tree_labels" = as.character(c(1:N)),
                   "trees" = trees,
                   "Theta" = Theta)

  tree_obj <- output_tree_df(tree_obj)
  res <- tree_df_to_trees(n, tree_obj$tree_df)
  tree_obj$tree_df <- res$tree_df
  tree_obj$trees <- res$trees

  # Perturb tree
  if (perturb) {
    res_perturb <- perturb_trees(n, tree_obj$tree_df, epsilon)
    tree_obj$tree_df <- res_perturb$tree_df
    tree_obj$trees <- res_perturb$trees
    tree_obj$epsilon <- epsilon
  }

  return(tree_obj)

}

generate_tree_MHN <- function(n, Theta, lambda_s) {

  t_s <- rexp(1, lambda_s) # sampling time
  current_pos <- c(1) # current positions being visited
  nodes <- c(0) # nodes in the tree T or A(T)
  path_times <- c(0) # waiting times for the pathways
  time_diffs <- c(0) # waiting time differences
  parents <- c(1) # list of positions for the parents
  children <- list(c()) # list of children
  in_tree <- c(TRUE) # whether the node is in T [TRUE] or A(T) [FALSE]
  genotypes <- rep(0,n)

  repeat{

    node_count <- length(current_pos)
    if (node_count == 0) {
      break
    }

    # Check if the tree is too large. Trees with size greater than the number of
    # mutational events (or half if n > 10) are typically not observed in practice.
    if ((n <= 10) && (sum(in_tree) > (n+1))) {
      return(list())
    } else if ((n > 10) && (sum(in_tree) > (n/2 + 1))) {
      return(list())
    }

    temp_pos <- c()
    for (i in c(1:node_count)) {

      path <- get_pathway(nodes, current_pos[i], parents)
      path <- setdiff(path, 0)

      for (j in setdiff(c(1:n), path)) {
        nodes <- c(nodes, j)
        parents <- c(parents, current_pos[i])
        children <- append(children, list(c()))
        children[[current_pos[i]]] <- c(children[[current_pos[i]]], length(nodes))

        if (length(path) == 0) {
          lambda_j <- exp(Theta[j,j])
        } else {
          lambda_j <- exp(Theta[j,j] + sum(sapply(path, function(k) Theta[j,k])))
        }
        time_j_diff <- rexp(1, lambda_j)
        time_diffs <- c(time_diffs, time_j_diff)
        time_j <- time_j_diff + path_times[current_pos[i]]
        path_times <- c(path_times, time_j)


        if (time_j < t_s) {
          in_tree <- c(in_tree, TRUE)
          temp_pos <- c(temp_pos, length(nodes))
          temp_genotype <- rep(0,n)
          temp_genotype[c(path,j)] <- 1
          genotypes <- rbind(genotypes,temp_genotype)
        } else {
          in_tree <- c(in_tree, FALSE)
        }
      }
    }

    current_pos <- temp_pos
  }

  tree <- list("nodes" = nodes,
               "parents" = parents,
               "children" = children,
               "in_tree" = in_tree,
               "genotypes" = genotypes)
  return(tree)

}

get_pathway <- function(nodes, pos, parents) {

  path <- c(nodes[pos])
  repeat {
    pa <- nodes[parents[pos]]
    if (pa == 0) {
      break
    }
    path <- c(pa, path)
    pos <- parents[pos]
  }
  return(path)

}


perturb_trees <- function(n, tree_df, epsilon = 0.05) {

  unique_tree_IDs <- unique(tree_df$Tree_ID)
  new_tree_df <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(new_tree_df) <- c("Tree_ID","Node_ID","Mutation_ID","Parent_ID")

  for (i in c(1:length(unique_tree_IDs))) {

    one_tree_df <- sort_one_tree_df(tree_df[tree_df$Tree_ID == unique_tree_IDs[i],])
    new_tree_df <- rbind(new_tree_df, perturb_tree(n, one_tree_df, epsilon))

  }

  res <- tree_df_to_trees(n, new_tree_df)
  return(res)

}

perturb_tree <- function(n, one_tree_df, epsilon = 0.05) {

  nr_nodes <- nrow(one_tree_df)
  tree_id <- one_tree_df$Tree_ID[1]
  node_to_remove <- c()

  # loop through the nodes except the root
  for (i in c(2:nr_nodes)) {

    node_id <- one_tree_df$Node_ID[i]
    parent <- one_tree_df$Parent_ID[i]
    children <- one_tree_df$Node_ID[which(one_tree_df$Parent_ID == node_id)]
    nr_ch <- length(children)

    switch(as.character(nr_ch),
           "0" = { #leaf node

             if (runif(1) < epsilon) { #perturb with an error rate of epsilon
               switch(as.character(sample.int(4,1)), #randomly select one of the three cases
                      "1" = { #add a parent
                        one_tree_df$Parent_ID[i] <- nrow(one_tree_df) + 1
                        one_tree_df <- rbind(one_tree_df, data.frame(Tree_ID = tree_id,
                                                                     Node_ID = nrow(one_tree_df) + 1,
                                                                     Mutation_ID = sample.int(n,1),
                                                                     Parent_ID = parent))
                      },
                      "2" = { #add a child
                        one_tree_df <- rbind(one_tree_df, data.frame(Tree_ID = tree_id,
                                                                     Node_ID = nrow(one_tree_df) + 1,
                                                                     Mutation_ID = sample.int(n,1),
                                                                     Parent_ID = node_id))
                      },
                      "3" = { #prune and reattach to grandparent/sibling
                        grandparent <- one_tree_df$Parent_ID[which(one_tree_df$Node_ID == parent)]
                        siblings <- setdiff(one_tree_df$Node_ID[one_tree_df$Parent_ID == parent], node_id)
                        to_attach <- c(grandparent, siblings)
                        if (length(to_attach) == 1) {
                          one_tree_df$Parent_ID[i] <- to_attach
                        } else {
                          one_tree_df$Parent_ID[i] <- sample(to_attach, 1)
                        }
                      },
                      "4" = { #be removed if not the only node in the tree
                        if ((nrow(one_tree_df) - length(node_to_remove)) > 2) {
                          node_to_remove <- c(node_to_remove, i)
                          one_tree_df$Node_ID[i] <- -1
                          one_tree_df$Parent_ID[i] <- -1
                        }
                      })
             }

           },
           "1" = { #interior node with one child

             if (runif(1) < epsilon) {
               switch(as.character(sample.int(4,1)),
                      "1" = { #add a parent
                        one_tree_df$Parent_ID[i] <- nrow(one_tree_df) + 1
                        one_tree_df <- rbind(one_tree_df, data.frame(Tree_ID = tree_id,
                                                                     Node_ID = nrow(one_tree_df) + 1,
                                                                     Mutation_ID = sample.int(n,1),
                                                                     Parent_ID = parent))
                      },
                      "2" = { #add a child
                        one_tree_df <- rbind(one_tree_df, data.frame(Tree_ID = tree_id,
                                                                     Node_ID = nrow(one_tree_df) + 1,
                                                                     Mutation_ID = sample.int(n,1),
                                                                     Parent_ID = node_id))
                      },
                      "3" = { #switch order with the child
                        one_tree_df$Parent_ID[which(one_tree_df$Node_ID  == children)] <- parent
                        one_tree_df$Parent_ID[i] <- children
                      },
                      "4" = { #be removed if not the only node in the tree
                        if ((nrow(one_tree_df) - length(node_to_remove)) > 2) {
                          node_to_remove <- c(node_to_remove, i)
                          one_tree_df$Parent_ID[which(one_tree_df$Node_ID  == children)] <- parent
                          one_tree_df$Node_ID[i] <- -1
                          one_tree_df$Parent_ID[i] <- -1
                        }
                      })
             }

           },
           { #interior node with more than one child

             if (runif(1) < epsilon) {
               switch(as.character(sample.int(2,1)),
                      "1" = { #add a parent
                        one_tree_df$Parent_ID[i] <- nrow(one_tree_df) + 1
                        one_tree_df <- rbind(one_tree_df, data.frame(Tree_ID = tree_id,
                                                                     Node_ID = nrow(one_tree_df) + 1,
                                                                     Mutation_ID = sample.int(n,1),
                                                                     Parent_ID = parent))
                      },
                      "2" = { #add a child
                        one_tree_df <- rbind(one_tree_df, data.frame(Tree_ID = tree_id,
                                                                     Node_ID = nrow(one_tree_df) + 1,
                                                                     Mutation_ID = sample.int(n,1),
                                                                     Parent_ID = node_id))
                      })
             }

           })

  }

  if (length(node_to_remove) > 0) {
    one_tree_df <- one_tree_df[-node_to_remove,]
  }

  return(one_tree_df)

}

