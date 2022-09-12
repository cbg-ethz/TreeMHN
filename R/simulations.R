
##' random_Theta(n, sparsity, exclusive_ratio)
##' This function randomly generates an MHN with log-transformed parameters.
##' @param n Number of mutational events
##' @param sparsity Percentage of off diagonal elements with a value of zero in the MHN. (Default: 0.5)
##' @param exclusive_ratio Percentage of non-zero elements that are negative,
##' meaning that one mutation is inhibiting another mutation (Default: 0.9)
##' @return An n-by-n matrix representing the MHN
##' @author Rudolf Schill, Xiang Ge Luo
##' @importFrom stats rgamma runif
##' @references Schill et al. (2020). Modelling cancer progression using Mutual Hazard Networks. Bioinformatics, 36(1), 241–249.
random_Theta <- function(n, sparsity = 0.5, exclusive_ratio = 0.5) {
  Theta  <- matrix(0, nrow=n, ncol=n)

  diag(Theta)  <- NA
  nonZeroSize <- (n^2 - n) * (1 - sparsity)
  nonZeros <- sample(which(!is.na(Theta)), size=nonZeroSize)
  exclusive <- sample(nonZeros, floor(nonZeroSize * exclusive_ratio))
  cooccur <- setdiff(nonZeros, exclusive)

  # Parameters (can change)
  diag(Theta) <- sample(c(runif(ceiling(n/2), -3, -1), runif(floor(n/2), -6, -3)))
  Theta[exclusive] <- - rgamma(length(exclusive), 4, 2.5)
  Theta[cooccur] <- rgamma(length(cooccur), 4, 2.5)

  return(round(Theta, 2))
}


##' @name generate_trees
##' @title Generate a set of trees based on a Mutual Hazard Network
##' @description This function generate a set of trees based on an MHN, which is either
##' provided by the user or generated randomly within the function.
##' @param n Number of mutational events (Default: 10).
##' @param N Number of trees (Default: 100).
##' @param lambda_s The rate of the sampling event (Default: 1).
##' @param Theta A mutual hazard network (Default: NULL for the function to generate).
##' @param sparsity Percentage of off diagonal elements with a value of zero in the MHN (Default: 0.5).
##' @param exclusive_ratio Percentage of non-zero elements that are negative,
##' meaning that one mutation is inhibiting another mutation (Default: 0.9).
##' @param smallest_tree_size The smallest size of the tree (Default: 1). 
##' The tree only containing the root node has a size of zero.
##' @param largest_tree_size The largest size of the tree. By default, the tree size
##' is no bigger than the number of mutations n (or half if n > 10).
##' @param mutations A list of mutation names, which must be unique values (Default: NULL).
##' @param perturb A flag for the function to perturb the generated tree structures (Default: FALSE).
##' @param epsilon Noise level for perturbing the trees (Default: 0.05).
##' @return A treeMHN object with the true MHN for generating the trees
##' @author Xiang Ge Luo
##' @export
generate_trees <- function(n = 10, N = 100, lambda_s = 1, Theta = NULL,
                           sparsity = 0.5, exclusive_ratio = 0.5, 
                           smallest_tree_size = 1, largest_tree_size = NULL,
                           mutations = NULL, perturb = FALSE, epsilon = 0.05) {
  
  if (is.null(Theta)) {
    Theta <- random_Theta(n, sparsity, exclusive_ratio)
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
  
  if (is.null(largest_tree_size)) {
    if (n <= 10) {
      largest_tree_size <- n
    } else {
      largest_tree_size <- floor(n / 2)
    }
  }
  
  tree_dfs <- list()
  tree_id <- 1
  while (length(tree_dfs) < N) {
    tree_df <- generate_tree_MHN(n, Theta, lambda_s, largest_tree_size)
    if (!is.null(tree_df)) {
      if (nrow(tree_df) - 1 >= smallest_tree_size) {
        tree_df$Patient_ID <- rep(tree_id, nrow(tree_df))
        tree_df$Tree_ID <- rep(tree_id, nrow(tree_df))
        tree_dfs <- append(tree_dfs, list(tree_df))
        tree_id <- tree_id + 1
      }
    }
  }
  tree_df <- do.call(rbind, tree_dfs) %>% 
    select(Patient_ID, Tree_ID, Node_ID, Mutation_ID, Parent_ID)
  
  if (perturb) {
    tree_df <- perturb_trees(n, tree_df, epsilon)
  }
  
  tree_obj <- input_tree_df(n = n,
                            tree_df = tree_df,
                            patients = as.character(c(1:N)), 
                            tree_labels = as.character(c(1:N)),
                            mutations = mutations)
  
  tree_obj$Theta <- Theta
  tree_obj$epsilon <- epsilon
  tree_obj$perturb <- perturb
  tree_obj$lambda_s <- lambda_s
  
  return(tree_obj)
  
}


##' @import dplyr
##' @importFrom stats rexp
generate_tree_MHN <- function(n, Theta, lambda_s, largest_tree_size) {
  
  # sampling time
  t_s <- rexp(1, lambda_s) 
  
  # initialize dataframe to store the tree
  tree_df <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(tree_df) <- c("Node_ID","Mutation_ID","Parent_ID", "Occurrence_Time")
  tree_df <- tree_df %>% bind_rows(c(Node_ID = 1, 
                                     Mutation_ID = 0, 
                                     Parent_ID = 1, 
                                     Occurrence_Time = 0))
  
  # current positions being visited
  current_nodes <- c(1)
  
  # recursively expand the tree
  while (length(current_nodes) > 0) {
    
    # check tree size
    if (nrow(tree_df) - 1 > largest_tree_size) {
      tree_df <- NULL
      break
    }
    
    # initialize next positions
    next_nodes <- c()
    
    # loop through all existing nodes
    for (i in c(1:length(current_nodes))) {
      
      node_i <- current_nodes[i]
      idx_i <- match(node_i, tree_df$Node_ID)
      pathway_i <- get_pathway_tree_df(tree_df, idx_i)
      t_i <- tree_df$Occurrence_Time[idx_i]
      
      # loop through all next nodes
      for (j in setdiff(c(1:n), pathway_i)) {
        
        # get rate
        if (length(pathway_i) == 0) {
          lambda_j <- exp(Theta[j, j])
        } else {
          lambda_j <- exp(Theta[j, j] + sum(sapply(pathway_i, function (k) Theta[j, k])))
        }
        
        # draw new waiting time
        t_j <- t_i + rexp(1, lambda_j)
        if (t_j < t_s) {
          # expand the tree by one node
          node_j <- nrow(tree_df) + 1
          tree_df <- tree_df %>% bind_rows(c(Node_ID = node_j, 
                                             Mutation_ID = j, 
                                             Parent_ID = node_i, 
                                             Occurrence_Time = t_j))
          # add this node to next positions
          next_nodes <- c(next_nodes, node_j)
        }
      }
    }
    
    current_nodes <- next_nodes
    
  }
  
  if (!is.null(tree_df)) {
    tree_df <- tree_df %>% select(Node_ID, Mutation_ID, Parent_ID)
  }
  
  return(tree_df)
  
}


perturb_trees <- function(n, tree_df, epsilon = 0.05) {

  unique_tree_IDs <- unique(tree_df$Tree_ID)
  new_tree_df <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(new_tree_df) <- c("Patient_ID", "Tree_ID","Node_ID","Mutation_ID","Parent_ID")

  for (i in c(1:length(unique_tree_IDs))) {

    one_tree_df <- sort_one_tree_df(tree_df[tree_df$Tree_ID == unique_tree_IDs[i],])
    new_tree_df <- rbind(new_tree_df, perturb_tree(n, one_tree_df, epsilon))

  }
  
  return(new_tree_df)

}

##' @importFrom stats runif
perturb_tree <- function(n, one_tree_df, epsilon = 0.05) {

  nr_nodes <- nrow(one_tree_df)
  patient_id <- one_tree_df$Patient_ID[1]
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
                        one_tree_df <- rbind(one_tree_df, data.frame(Patient_ID = patient_id,
                                                                     Tree_ID = tree_id,
                                                                     Node_ID = nrow(one_tree_df) + 1,
                                                                     Mutation_ID = sample.int(n,1),
                                                                     Parent_ID = parent))
                      },
                      "2" = { #add a child
                        one_tree_df <- rbind(one_tree_df, data.frame(Patient_ID = patient_id,
                                                                     Tree_ID = tree_id,
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
                        one_tree_df <- rbind(one_tree_df, data.frame(Patient_ID = patient_id,
                                                                     Tree_ID = tree_id,
                                                                     Node_ID = nrow(one_tree_df) + 1,
                                                                     Mutation_ID = sample.int(n,1),
                                                                     Parent_ID = parent))
                      },
                      "2" = { #add a child
                        one_tree_df <- rbind(one_tree_df, data.frame(Patient_ID = patient_id,
                                                                     Tree_ID = tree_id,
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
                        one_tree_df <- rbind(one_tree_df, data.frame(Patient_ID = patient_id,
                                                                     Tree_ID = tree_id,
                                                                     Node_ID = nrow(one_tree_df) + 1,
                                                                     Mutation_ID = sample.int(n,1),
                                                                     Parent_ID = parent))
                      },
                      "2" = { #add a child
                        one_tree_df <- rbind(one_tree_df, data.frame(Patient_ID = patient_id,
                                                                     Tree_ID = tree_id,
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

