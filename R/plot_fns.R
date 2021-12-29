require(DiagrammeR)
require(ggplot2)

##' @name plot_tree
##' @title Plot a tree
##' @description This function plots a tree in named list format.
##' @param tree A tree in named list format.
##' @param mutations A list of mutation names corresponding to the mutation IDs.
##' @param tree_label The title of the tree (Default: NULL).
##' @author Xiang Ge Luo
##' @export
plot_tree <- function(tree, mutations, tree_label = NULL) {

  # graphviz dot language
  graph_dot <- "
  digraph g {
  labelloc='t';
  fontname='Helvetica bold';
  fontsize=28;
  "

  # Add title
  if (is.null(tree_label)) {
    graph_dot <- paste(graph_dot, "label = 'Tree", tree$tree_ID, "';")
  } else {
    graph_dot <- paste(graph_dot, "label = '", tree_label, "';")
  }

  # Add nodes
  nr_nodes <- sum(tree$in_tree)
  node_idx <- c(1:length(tree$nodes))[tree$in_tree]
  node_labels <- c("Root", mutations[tree$nodes[tree$in_tree][-1]])
  for (i in c(1:nr_nodes)) {
    graph_dot <- paste(graph_dot, node_idx[i], "[label = '", node_labels[i], "'];")
  }

  # Add edges
  parents <- tree$parents[tree$in_tree]
  for (i in c(2:nr_nodes)) {
    graph_dot <- paste(graph_dot, parents[i], "->", node_idx[i], ";")
  }

  graph_dot <- paste(graph_dot,"}")
  g <- grViz(graph_dot)
  return(g)

}

##' @name plot_pathways
##' @title Plot most probable pathways of a given length computed from a Mutual Hazard Network
##' @description This function takes a Mutual Hazard Network as input and plot the most probable
##' mutational pathways of a given length.
##' @param Theta An n-by-n matrix representing a Mutual Hazard Network
##' @param mutations A list of mutation names, which must be unique values.
##' If no names are given, then the mutation IDs will be used.
##' @param n_order The length of the pathways (Default: 4).
##' @param top_M Number of most probable pathways to plot (Default: 10).
##' @param log2 A boolean flag indicating whether the x axis is scaled by log2.
##' @author Xiang Ge Luo
##' @export
plot_pathways <- function(Theta, mutations = NULL, n_order = 4, top_M = 10, log2 = TRUE) {

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
  pathways <- Theta_to_pathways(Theta, n_order = n_order, prob_only = FALSE)
  pathways <- pathways[c(1:top_M),]
  pathways[,c((n_order+2):(2*n_order+1))] <- t(apply(pathways[,c((n_order+2):(2*n_order+1))],1,cumsum))
  pathways$prob <- pathways$prob*100
  waiting_time <- c()
  probability <- c()
  labels <- c()
  for (i in c(1:n_order)) {
    waiting_time <- c(waiting_time, pathways[,(i + n_order + 1)])
    probability <- c(probability, paste0(round(pathways$prob,3),"%"))
    labels <- c(labels, mutations[pathways[,i]])
  }
  df <- data.frame(waiting_time,probability,labels)

  g <- ggplot(df, aes(x = waiting_time, y = factor(probability), label = labels)) +
    geom_label(aes(fill = factor(labels))) +
    xlab("Expected waiting time relative to the sampling rate") + ylab("Probability") +
    theme_classic() + guides(fill="none") +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size=16),
          axis.text = element_text(size=14),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16))

  if (log2) {
    g <- g + scale_x_continuous(trans='log2')
  }

  return(g)

}

##' @name next_mutation
##' @title Plot the next most probable mutational events
##' @description Given a particular tree and a Mutual Hazard Network, this function 
##' finds the next most probable mutational events.
##' @param n Number of mutational events.
##' @param tree A tree in named list format.
##' @param mutations A list of mutation names, which must be unique values.
##' If no names are given, then the mutation IDs will be used.
##' @param tree_label The title of the tree (Default: NULL).
##' @param top_M Number of most probable mutational events to plot (Default: 1).
##' @export
next_mutation <- function(n, tree, Theta, mutations = NULL, tree_label = NULL, top_M = 1) {

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
  next_lambdas <- rep(0, length(next_mutations))

  for (i in c(1:length(next_mutations))) {

    pos <- next_mutations[i]
    node <- get_pathway(tree$nodes,pos,tree$parents)
    next_lambdas[i] <- get_lambda(node, Theta)

  }

  probs <- next_lambdas / sum(next_lambdas)

  cat("Top",top_M,"most probable mutational events that will happen next:\n")
  top_M_idx <- order(probs, decreasing = TRUE)[c(1:top_M)]
  temp <- tree
  for (i in c(1:top_M)) {

    idx <- top_M_idx[i]
    most_likely_pos <- next_mutations[idx]
    node <- get_pathway(tree$nodes,most_likely_pos,tree$parents)
    node <- mutations[node]
    cat("The next most probable node:", paste(c("Root",node),collapse = "->"), "\n")
    cat("Probability:", round(probs[idx]*100,3), "%\n")
    temp$in_tree[most_likely_pos] <- TRUE

  }

  temp <- output_one_tree_df(temp)
  new_tree <- tree_df_to_trees(n, temp)$trees[[1]]

  g <- plot_tree(new_tree,mutations,tree_label)
  return(g)

}
