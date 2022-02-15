##' @import Matrix
obj_grad_helper <- function(n, N, trees, Theta, lambda_s) {
  tr_mat_vec <- list()
  comp_geno_vec <- list()
  node_labels_vec <- list()
  log_prob_vec <- c()

  for (i in c(1:N)) {
    tree <- trees[[i]]
    p <- build_poset(tree)
    comp_geno <- compatible_genotypes(p)
    comp_geno_vec <- append(comp_geno_vec, list(comp_geno))
    node_labels <- get_node_labels(p)
    node_labels_vec <- append(node_labels_vec, list(node_labels))
    tr_mat <- build_tr_mat(n, Theta, comp_geno, node_labels, lambda_s)
    tr_mat_vec <- append(tr_mat_vec, list(tr_mat))
    log_prob <- compute_obs_ll(tr_mat, lambda_s)
    log_prob_vec <- c(log_prob_vec, log_prob)
  }

  res <- list("tr_mat_vec" = tr_mat_vec,
              "comp_geno_vec" = comp_geno_vec,
              "node_labels_vec" = node_labels_vec,
              "log_prob_vec" = log_prob_vec)
  return(res)
}

compatible_genotypes <- function(poset) {
  if (is.null(poset)) {
    return(matrix(0))
  } else {
    nr_nodes <- ncol(poset)
    genotypes <- matrix(rep(0, nr_nodes), ncol = nr_nodes)
    for (j in c(1:nr_nodes)) {
      temp <- genotypes
      pa_i <- which(poset[,j] == 1)
      if (length(pa_i) == 0) {
        for (k in c(1:nrow(temp))) {
          g <- matrix(temp[k,],ncol = nr_nodes)
          g[1,j] <- 1
          genotypes <- rbind(genotypes,g)
        }
      } else {
        for (k in c(1:nrow(temp))) {
          g <- matrix(temp[k,],ncol = nr_nodes)
          if (g[1,pa_i] == 1) {
            g[1,j] <- 1
            genotypes <- rbind(genotypes,g)
          }
        }
      }
    }
    if (nr_nodes > 1) {
      return(genotypes[order(apply(genotypes,1,sum)),])
    } else {
      return(genotypes)
    }
  }
}

get_node_labels <- function(poset) {
  if (is.null(poset)) {
    return(NULL)
  } else {
    mut <- as.integer(rownames(poset))
    nr_nodes <- nrow(poset)
    node_labels <- list()
    for (i in c(1:nr_nodes)) {
      label <- c(mut[i])
      j <- i
      repeat{
        pa_i <- which(poset[,j] == 1)
        if (length(pa_i) == 0) {
          break
        } else {
          label <- c(mut[pa_i], label)
          j <- pa_i
        }
      }
      node_labels <- c(node_labels, list(label))
    }
    return(node_labels)
  }
}

##' @import ggm
build_poset <- function(tree) {
  in_tree_idx <- which(tree$in_tree)[-1]
  nodes_in_tree <- tree$nodes[in_tree_idx]
  nr_nodes <- length(in_tree_idx)
  if (nr_nodes == 0) {
    return(NULL)
  } else {
    poset <- matrix(0, nrow = nr_nodes, ncol = nr_nodes)
    for (i in c(1:nr_nodes)) {
      pa_idx <- tree$parents[in_tree_idx[i]]
      if (pa_idx > 1) {
        pa_i <- which(in_tree_idx == pa_idx);
        poset[pa_i, i] <- 1
      }
    }
    dimnames(poset) <- list(c(nodes_in_tree),c(nodes_in_tree))
    if (nr_nodes > 1) {
      return(ggm::topSort(poset))
    } else {
      return(poset)
    }
  }
}

check_mutations <- function(n, trees) {
  all_muts <- c(1:n)
  for (tree in trees) {
    all_muts <- setdiff(all_muts, tree$nodes[tree$in_tree])
    if (length(all_muts) == 0) {
      return(TRUE)
    }
  }
  return(FALSE)
}


