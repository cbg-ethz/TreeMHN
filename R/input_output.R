
##' input_tree_df(n, tree_df, mutations, tree_labels)
##' This function processes a data frame of mutation trees and output a TreeMHN object
##' @param n Number of mutational events
##' @param tree_df A data frame with the following columns:
##' \itemize{
##' \item Patient_ID: IDs of patients, unique for each patient
##' \item Tree_ID: IDs of mutation trees, unique for each tree
##' \item Node_ID: IDs of each node in the tree, including the root node (with ID "1"), unique for each node
##' \item Mutation_ID: IDs of each mutational event, the root node has a mutation ID of "0",
##' other mutation IDs can be duplicated in the tree to allow for parallel mutations
##' \item Parent_ID: IDs of the parent node ID. The root node has itself as parent (ID "1").
##' }
##' @param patients A list of patient labels, which must be unique values.
##' If no labels are given, then the patient IDs will be used.
##' @param tree_labels A list of tree labels, which must be unique values.
##' If no labels are given, then the tree IDs will be used.
##' @param mutations A list of mutation names, which must be unique values.
##' If no names are given, then the mutation IDs will be used.
##' @param weights Weights of the trees. If no values are given, weights are
##' assigned equally to the trees such that each patient has a weight of 1.
##' @return A TreeMHN object
##' @author Xiang Ge Luo
##' @export
input_tree_df <- function(n, tree_df, patients = NULL, tree_labels = NULL, 
                          mutations = NULL, weights = NULL) {

  # Check input format
  if (length(setdiff(colnames(tree_df), c("Patient_ID","Tree_ID","Node_ID","Mutation_ID","Parent_ID"))) != 0) {
    stop("Please check the input format of the data frame...")
  }
  
  # Check patient labels
  N_patients <- length(unique(tree_df$Patient_ID))
  if (is.null(patients)) {
    patients <- as.character(seq(1, N_patients))
  } else {
    if (length(unique(patients)) != N_patients) {
      stop("The list of patient labels has length different from the number of patients. Please check again...")
    }
  }
  
  # Check tree labels
  N <- length(unique(tree_df$Tree_ID))
  if (is.null(tree_labels)) {
    tree_labels <- as.character(seq(1, N))
  } else {
    if (length(unique(tree_labels)) != N) {
      stop("The list of tree labels has length different from the number of trees. Please check again...")
    }
  }
  
  # Check mutation names
  if (is.null(mutations)) {
    mutations <- as.character(seq(1, n))
  } else {
    if (length(mutations) != n) {
      stop("The list of mutations has length different from n. Please check again...")
    } else if (length(unique(mutations)) != n) {
      stop("Mutation names must be unique. Please check again...")
    }
  }
  
  # Check tree weights
  if (is.null(weights)) {
    weights <- tree_df %>% 
      select(Patient_ID, Tree_ID) %>% 
      distinct(Tree_ID, .keep_all = TRUE) %>% 
      group_by(Patient_ID) %>%
      mutate(temp = n()) %>%
      mutate(weights = 1 / temp) %>%
      ungroup(Patient_ID) %>% 
      .$weights
  } else {
    patient_level_weights <- tree_df %>% 
      select(Patient_ID, Tree_ID) %>% 
      distinct(Tree_ID, .keep_all = TRUE) %>%
      mutate(tree_weight = weights[Tree_ID]) %>% 
      group_by(Patient_ID) %>%
      summarise(tree_weight = sum(tree_weight)) %>%
      .$tree_weight
    
    idx <- which(patient_level_weights != 1)
    if (length(idx) != 0) {
      stop(paste("The tree weights of patients with ID", 
                 paste(idx, collapse = ", "), 
                 "do not sum to 1. Please check again..."))
    }
  }

  # Convert data frame into tree format
  res <- tree_df_to_trees(n, tree_df)

  tree_obj <- list("n" = n,
                   "N" = N,
                   "N_patients" = N_patients,
                   "mutations" = mutations,
                   "tree_labels" = tree_labels,
                   "patients" = patients,
                   "tree_df" = res$tree_df,
                   "trees" = res$trees,
                   "weights" = weights)
  return(tree_obj)
}

##' tree_df_to_trees(n, tree_df)
##' This function processes a data frame of mutation trees
##' and output a list of processed mutation trees
##' @param n Number of mutational events
##' @param tree_df A data frame with the following columns:
##' \itemize{
##' \item Tree_ID: IDs of mutation trees, unique for each patient
##' \item Node_ID: IDs of each node in the tree, including the root node (with ID "1"), unique for each node
##' \item Mutation_ID: IDs of each mutational event, the root node has a mutation ID of "0",
##' other mutation IDs can be duplicated in the tree to allow for parallel mutations
##' \item Parent_ID: IDs of the parent node ID. The root node has itself as parent (ID "1").
##' }
##' @return A named list of a tree data frame and the corresponding trees
tree_df_to_trees <- function(n, tree_df) {
  # Convert data frame into tree format
  trees <- list()
  unique_tree_IDs <- unique(tree_df$Tree_ID)
  for (i in c(1:length(unique_tree_IDs))) {

    # Ensure that the node IDs are in ascending order
    one_tree_df <- sort_one_tree_df(tree_df[tree_df$Tree_ID == unique_tree_IDs[i],])
    tree_df[tree_df$Tree_ID == unique_tree_IDs[i],] <- one_tree_df

    # Construct the tree T
    tree_to_add <- list()
    tree_to_add$patient_ID <- one_tree_df$Patient_ID[1] ## New
    tree_to_add$tree_ID <- unique_tree_IDs[i]
    tree_to_add$nodes <- one_tree_df$Mutation_ID
    tree_to_add$parents <- one_tree_df$Parent_ID
    tree_to_add$children <- replicate(nrow(one_tree_df), list(integer(0)))
    for (j in c(1:nrow(one_tree_df))) {
      temp <- which(one_tree_df$Parent_ID == j)
      if (length(temp) != 0) {
        tree_to_add$children[[j]] <- temp
      }
    }
    tree_to_add$children[[1]] <- tree_to_add$children[[1]][-1]
    if (length(tree_to_add$children[[1]]) == 0) {
      stop(paste("Tree with ID",unique_tree_IDs[i],"does not contain edges from the root. Please check again..."))
    }
    tree_to_add$in_tree <- rep(TRUE, nrow(one_tree_df))

    tree_to_add$genotypes <- tree_df_to_genotypes(n, one_tree_df)

    trees <- append(trees, list(tree_to_add))
  }

  res <- list()
  res$tree_df <- tree_df

  # Construct the augmented trees A(T)
  res$trees <- get_augmented_trees(n, trees)
  return(res)
}

##' tree_df_to_trees(one_tree_df)
##' A helper function that sorts the node IDs such that they are in ascending order
##' @param n Number of mutational events
##' @param one_tree_df A tree data frame for a single patient/tumor
##' @return A sorted tree data frame
sort_one_tree_df <- function(one_tree_df) {
  one_tree_df$Parent_ID <- sapply(c(1:nrow(one_tree_df)),
                                  function (j) which(one_tree_df$Node_ID == one_tree_df$Parent_ID[j]))
  one_tree_df$Node_ID <- seq(1, nrow(one_tree_df))
  return(one_tree_df)
}

##' output_tree_df(tree_obj)
##' A helper function that construct a tree data frame from a list of mutation trees
##' @param tree_obj A TreeMHN object
##' @return A TreeMHN object with the tree data frame
output_tree_df <- function(tree_obj) {

  tree_df <- data.frame(matrix(ncol = 5, nrow = 0)) ## New

  for (tree in tree_obj$trees) {
    tree_df <- rbind(tree_df, output_one_tree_df(tree))
  }

  colnames(tree_df) <- c("Patient_ID","Tree_ID","Node_ID","Mutation_ID","Parent_ID") ## New
  tree_obj$tree_df <- tree_df
  return(tree_obj)

}

output_one_tree_df <- function(tree) {
  nr_nodes <- sum(tree$in_tree)
  patient_ID_temp <- rep(tree$patient_ID, nr_nodes) ## New
  tree_ID_temp <- rep(tree$tree_ID, nr_nodes)
  node_ID_temp <- seq(1,length(tree$nodes))[tree$in_tree]
  mutation_ID_temp <- tree$nodes[tree$in_tree]
  parent_ID_temp <- tree$parents[tree$in_tree]
  temp_df <- cbind(patient_ID_temp,tree_ID_temp, node_ID_temp, mutation_ID_temp, parent_ID_temp) ## New
  colnames(temp_df) <- c("Patient_ID","Tree_ID","Node_ID","Mutation_ID","Parent_ID") ## New
  return(sort_one_tree_df(data.frame(temp_df)))
}

##' tree_df_to_genotypes(n, tree_df)
##' This function processes a data frame of mutation trees
##' and output a list of subclonal genotypes.
##' Note that the wild type (a vector of zeros) is also included.
##' @param n Number of mutational events
##' @param tree_df A data frame with the following columns:
##' \itemize{
##' \item Tree_ID: IDs of mutation trees, unique for each patient
##' \item Node_ID: IDs of each node in the tree, including the root node (with ID "1"), unique for each node
##' \item Mutation_ID: IDs of each mutational event, the root node has a mutation ID of "0",
##' other mutation IDs can be duplicated in the tree to allow for parallel mutations
##' \item Parent_ID: IDs of the parent node ID. The root node has itself as parent (ID "1").
##' }
##' @return A matrix of subclonal genotypes
tree_df_to_genotypes <- function(n, tree_df) {

  genotypes <- matrix(0, nrow = nrow(tree_df), ncol = n)
  check_change <- -1
  repeat {

    if (sum(genotypes) == check_change) {
      break
    } else {
      check_change <- sum(genotypes)
    }

    for (i in c(2:nrow(tree_df))) {
      pa_idx <- which(tree_df$Node_ID == tree_df$Parent_ID[i])
      genotypes[i,] <- genotypes[pa_idx,]
      genotypes[i, tree_df$Mutation_ID[i]] <- 1
    }

  }
  return(genotypes)

}
