
##' @name plot_tree_list
##' @title Plot a tree in named list format
##' @description This function plots a tree in named list format.
##' @param tree A tree in named list format.
##' @param mutations A list of mutation names corresponding to the mutation IDs.
##' @param tree_label The title of the tree (Default: NULL).
##' @author Xiang Ge Luo
##' @import DiagrammeR
##' @export
plot_tree_list <- function(tree, mutations, tree_label = NULL) {
  
  # graphviz dot language
  graph_dot <- "
  digraph g {
  labelloc='t';
  fontname='Arial';
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
    graph_dot <- paste(graph_dot, 
                       node_idx[i], 
                       "[label = '", node_labels[i], 
                       "', fontname='Arial'];")
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

##' @name plot_tree_df
##' @title Plot a tree in data frame format
##' @description This function plots a tree in data frame format.
##' @param tree_df A tree in data frame format.
##' @param mutations A list of mutation names corresponding to the mutation IDs.
##' @param tree_label The title of the tree (Default: NULL).
##' @param main_node_color The color of the nodes in the main tree structure (Default: pale turquoise).
##' @param next_node_color The color of the nodes in the augmented tree structure (Default: thistle).
##' @author Xiang Ge Luo
##' @import DiagrammeR
##' @export
plot_tree_df <- function(tree_df, mutations, tree_label = NULL, main_node_color = "paleturquoise3", next_node_color = "thistle") {
  
  # graphviz dot language
  graph_dot <- "
  digraph g {
  labelloc='t';
  fontname='Arial';
  fontsize=28;
  "
  
  # Add title
  if (is.null(tree_label)) {
    graph_dot <- paste(graph_dot, "label = 'Tree", tree_df$Tree_ID[1], "';")
  } else {
    graph_dot <- paste(graph_dot, "label = '", tree_label, "';")
  }
  
  # Add nodes
  node_labels <- c("Root", mutations[tree_df$Mutation_ID[-1]])
  if ("Existing" %in% colnames(tree_df)) {
    for (i in c(1:nrow(tree_df))) {
      if (tree_df$Existing[i]) {
        graph_dot <- paste(graph_dot, 
                           tree_df$Node_ID[i], 
                           "[label = '", 
                           node_labels[i], 
                           "', fontname='Arial', style=filled, color=",
                           main_node_color,
                           "];")
      } else {
        if ("Prob" %in% colnames(tree_df)) {
          graph_dot <- paste(graph_dot, 
                             tree_df$Node_ID[i], 
                             "[label = '", node_labels[i], "\n", tree_df$Prob[i], "%",
                             "', fontname='Arial', style=filled, color=",
                             next_node_color,
                             "];")
        } else {
          graph_dot <- paste(graph_dot, 
                             tree_df$Node_ID[i], 
                             "[label = '", node_labels[i],
                             "', fontname='Arial', style=filled, color=",
                             next_node_color,
                             "];")
        }
        
      }
    }
  } else {
    for (i in c(1:nrow(tree_df))) {
      graph_dot <- paste(graph_dot, 
                         tree_df$Node_ID[i], 
                         "[label = '", 
                         node_labels[i], 
                         "', fontname='Arial', style=filled, color=",
                         main_node_color,
                         "];")
    }
  }
  
  # Add edges
  for (i in c(2:nrow(tree_df))) {
    graph_dot <- paste(graph_dot, 
                       tree_df$Parent_ID[i], 
                       "->", 
                       tree_df$Node_ID[i], 
                       ";")
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
##' @import ggplot2
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
  if(n_order > 1) {
    pathways[,c((n_order+2):(2*n_order+1))] <- t(apply(pathways[,c((n_order+2):(2*n_order+1))],1,cumsum))
  }
  pathways$prob <- pathways$prob*100
  waiting_time <- c()
  probability <- c()
  labels <- c()
  for (i in c(1:n_order)) {
    waiting_time <- c(waiting_time, pathways[,(i + n_order + 1)])
    probability <- c(probability, round(pathways$prob, 3))
    labels <- c(labels, mutations[pathways[,i]])
  }
  df <- data.frame(waiting_time,probability,labels)
  
  g <- ggplot(df, aes(x = waiting_time, y = factor(probability, ordered = TRUE), label = labels)) +
    geom_label(aes(fill = factor(labels))) +
    xlab("Expected waiting time relative to the sampling rate") + ylab("Probability") +
    theme_classic() + guides(fill="none") +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size=16),
          axis.text = element_text(size=14),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16)) +
    scale_y_discrete(labels=sapply(sort(pathways$prob), function(x) paste0(round(x, 3), "%")))
  
  if (log2) {
    g <- g + scale_x_continuous(trans='log2')
  }
  
  return(g)
  
}

##' @name plot_next_mutations
##' @title Plot the next most probable mutational events
##' @description Given a particular tree and a Mutual Hazard Network, this function
##' finds the next most probable mutational events.
##' @param n Number of mutational events.
##' @param tree_df A tree in data frame format.
##' @param Theta An n-by-n matrix representing a Mutual Hazard Network
##' @param mutations A list of mutation names, which must be unique values.
##' If no names are given, then the mutation IDs will be used.
##' @param tree_label The title of the tree (Default: NULL).
##' @param top_M Number of most probable mutational events to plot (Default: 5).
##' @export
plot_next_mutations <- function(n, tree_df, Theta,
                                mutations = NULL, tree_label = NULL, top_M = 5) {
  
  if (is.null(mutations)) {
    mutations <- as.character(seq(1,n))
  } else {
    if (length(mutations) != n) {
      stop("The number of mutations doesn't match matrix dimension. Please check again...")
    } else if (length(unique(mutations)) != n) {
      stop("Mutation names must be unique. Please check again...")
    }
  }
  
  tree_df$Existing <- rep(1, nrow(tree_df))
  tree_df$Prob <- rep(0, nrow(tree_df))
  
  next_mutations <- c()
  next_lambdas <- c()
  next_parents <- c()
  for (i in c(1:nrow(tree_df))) {
    pathway_i <- get_pathway_tree_df(tree_df, i)
    siblings <- setdiff(tree_df$Mutation_ID[tree_df$Parent_ID == tree_df$Node_ID[i]], c(0))
    for (j in setdiff(c(1:n), c(pathway_i, siblings))) {
      next_mutations <- c(next_mutations, j)
      next_lambdas <- c(next_lambdas, get_lambda(c(pathway_i, j), Theta))
      next_parents <- c(next_parents, tree_df$Node_ID[i])
    }
  }
  
  cat("Top", top_M, "most probable mutational events that will happen next:\n")
  probs <- next_lambdas / sum(next_lambdas)
  llr <- log(probs * length(probs))
  top_M_idx <- order(probs, decreasing = TRUE)[c(1:top_M)]
  for (i in c(1:top_M)) {
    idx <- top_M_idx[i]
    tree_df <- rbind(tree_df, data.frame(Patient_ID = tree_df$Patient_ID[1],
                                         Tree_ID = tree_df$Tree_ID[1],
                                         Node_ID = nrow(tree_df) + 1,
                                         Mutation_ID = next_mutations[idx],
                                         Parent_ID = next_parents[idx],
                                         Existing = 0,
                                         Prob = round(probs[idx]*100, 3)))
    node <- get_pathway_tree_df(tree_df, nrow(tree_df))
    node <- mutations[node]
    cat("The next most probable node:", paste(c("Root", node),collapse = "->"), "\n")
    cat("Probability:", round(probs[idx]*100, 3), "%\n")
    cat("Log ratio model vs random:", round(llr[idx], 3), "\n")
  }
  
  g <- plot_tree_df(tree_df, mutations, tree_label)
  return(g)
  
}


##' @name plot_pathways_w_sampling
##' @title Plot most probable pathways before a sampling event computed from a Mutual Hazard Network
##' @description This function takes a Mutual Hazard Network as input and plot the most probable
##' mutational pathways before a sampling event.
##' @param Theta An n-by-n matrix representing a Mutual Hazard Network.
##' @param mutations A list of mutation names, which must be unique values.
##' If no names are given, then the mutation IDs will be used.
##' @param top_M Number of most probable pathways to plot (Default: 10).
##' @param lambda_s Sampling rate (Default: 1)
##' @param mutation_colors A named vector with the color codes for all mutations (Default: NULL)
##' @param prob_digits Number of digits to show for the probabilities (Default: 2)
##' @param full A boolean flag indicating whether the full MHN is used.
##' If false, then only mutations with non-zero off-diagonal entries will be used (default: FALSE).
##' @author Xiang Ge Luo
##' @import ggplot2
##' @importFrom stats rnorm
##' @export
plot_pathways_w_sampling <- function(Theta, mutations, top_M = 10, lambda_s = 1, 
                                     mutation_colors = NULL, prob_digits = 2, full = FALSE) {
  
  if (full) {
    to_use_mat <- Theta
  } else {
    to_use <- sapply(c(1:nrow(Theta)), function (i) any(Theta[i,-i] != 0) || any(Theta[-i,i] != 0))
    to_use_mat <- Theta[to_use, to_use]
    mutations <- mutations[to_use]
    mutation_colors <- mutation_colors[to_use]
  }
  
  n <- nrow(to_use_mat)
  
  pathway_df <- Theta_to_pathways_w_sampling(to_use_mat, top_M, lambda_s)
  pathway_df$probs <- pathway_df$probs*100 + rnorm(top_M, sd = 0.0001) # Avoid overlapping labels
  
  waiting_time <- c()
  probability <- c()
  labels <- c()
  probability2 <- c()
  line_start <- c()
  line_end <- c()
  
  for (i in c(1:top_M)) {
    p <- as.integer(unlist(strsplit(pathway_df$pathways[i], "_")))
    labels <- c(labels, "Root", mutations[p])
    probability <- c(probability, rep(pathway_df$probs[i], length(p) + 1))
    probability2 <- c(probability2, rep(pathway_df$probs[i], length(p)))
    times <- rep(0, length(p))
    for (j in c(1:length(p))) {
      pp <- p[c(1:j)]
      denom_set <- get_children(n, pp[-length(pp)])
      denom <- lambda_s + sum(sapply(denom_set, function (l) get_lambda(l, to_use_mat)))
      time_diff <- 1 / denom
      if (j == 1) {
        times[j] <- time_diff
        line_start <- c(line_start, 0.1)
        line_end <- c(line_end, time_diff - nchar(mutations[p[j]]) / 50)
      } else {
        times[j] <- times[j - 1] + time_diff
        line_start <- c(line_start, times[j - 1] + nchar(mutations[p[j - 1]]) / 50)
        line_end <- c(line_end, times[j - 1] + time_diff - nchar(mutations[p[j]]) / 50)
      }
    }
    waiting_time <- c(waiting_time, 0, times)
  }
  
  df <- data.frame(waiting_time, probability, labels)
  df2 <- data.frame(probability2, line_start, line_end)
  
  g <- ggplot(df, aes(x = waiting_time, y = factor(probability, ordered = TRUE), label = labels)) +
    geom_label(aes(fill = factor(labels))) +
    xlab("Expected waiting time relative to the sampling rate") + ylab("Probability") +
    theme_classic() + guides(fill="none") +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size=16),
          axis.text = element_text(size=14),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16)) +
    geom_segment(data = df2, 
                 mapping = aes(x = line_start, 
                               xend = line_end,
                               y = factor(probability2, ordered = TRUE),
                               yend = factor(probability2, ordered = TRUE),
                               label = NULL),
                 arrow = arrow(length = unit(0.1, "cm"))) +
    scale_y_discrete(labels=sapply(sort(pathway_df$probs), function(x) paste0(round(x, prob_digits), "%"))) +
    xlim(c(0, max(df$waiting_time)*1.1))
  
  if (!is.null(mutation_colors)) {
    mutation_colors <- c("#E0E0E0", mutation_colors)
    names(mutation_colors) <- c("Root", mutations)
    g <- g + scale_fill_manual(values = mutation_colors)
  }
  
  return(g)
}


##' @name plot_observed_pathways
##' @title Plot most frequent observed pathways from a cohort of mutation trees
##' @description This function takes a cohort of mutation trees and the estimated
##' Mutual Hazard Network as input and plot the most frequent observed pathways.
##' @param tree_obj A TreeMHN object.
##' @param Theta An n-by-n matrix representing a Mutual Hazard Network.
##' @param top_M Number of most frequent pathways to plot (Default: 10).
##' @param lambda_s Sampling rate (Default: 1)
##' @param mutation_colors A named vector with the color codes for all mutations (Default: NULL)
##' @param prob_digits Number of digits to show for the probabilities (Default: 2)
##' @param at_least_twice A boolean flag indicating whether only trajectories that
##' appear at least twice in the dataset are shown (Default: TRUE). If true, this flag
##' will overwrite the top_M argument.
##' @author Xiang Ge Luo
##' @importFrom stats rnorm
##' @import dplyr
##' @export
plot_observed_pathways <- function(tree_obj, Theta, top_M = 10, lambda_s = 1, 
                                   mutation_colors = NULL, prob_digits = 2, at_least_twice = TRUE) {
  
  n <- nrow(Theta)
  mutations <- tree_obj$mutations
  pathway_df <- get_observed_pathways(tree_obj)
  
  if (at_least_twice) {
    pathway_df <- pathway_df %>% filter(count >= 2)
    top_M <- nrow(pathway_df)
  } else {
    pathway_df <- pathway_df[c(1:top_M),]
  }
  
  pathway_df$probs <- pathway_df$probs*100 + rnorm(top_M, sd = 0.0001) # Avoid overlapping labels
  
  waiting_time <- c()
  probability <- c()
  labels <- c()
  probability2 <- c()
  line_start <- c()
  line_end <- c()
  
  for (i in c(1:top_M)) {
    p <- as.integer(unlist(strsplit(pathway_df$pathways[i], "_")))
    labels <- c(labels, "Root", mutations[p])
    probability <- c(probability, rep(pathway_df$probs[i], length(p) + 1))
    probability2 <- c(probability2, rep(pathway_df$probs[i], length(p)))
    times <- rep(0, length(p))
    for (j in c(1:length(p))) {
      pp <- p[c(1:j)]
      denom_set <- get_children(n, pp[-length(pp)])
      denom <- lambda_s + sum(sapply(denom_set, function (l) get_lambda(l, Theta)))
      time_diff <- 1 / denom
      if (j == 1) {
        times[j] <- time_diff
        line_start <- c(line_start, 0.1)
        line_end <- c(line_end, time_diff - nchar(mutations[p[j]]) / 50)
      } else {
        times[j] <- times[j - 1] + time_diff
        line_start <- c(line_start, times[j - 1] + nchar(mutations[p[j - 1]]) / 50)
        line_end <- c(line_end, times[j - 1] + time_diff - nchar(mutations[p[j]]) / 50)
      }
    }
    waiting_time <- c(waiting_time, 0, times)
  }
  
  df <- data.frame(waiting_time, probability, labels)
  df2 <- data.frame(probability2, line_start, line_end)
  
  g <- ggplot(df, aes(x = waiting_time, y = factor(probability, ordered = TRUE), label = labels)) +
    geom_label(aes(fill = factor(labels))) +
    xlab("Expected waiting time relative to the sampling rate") + ylab("Relative frequency") +
    theme_classic() + guides(fill="none") +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_text(size=16),
          axis.text = element_text(size=14),
          legend.text = element_text(size=16),
          legend.title = element_text(size=16)) +
    geom_segment(data = df2, 
                 mapping = aes(x = line_start, 
                               xend = line_end,
                               y = factor(probability2, ordered = TRUE),
                               yend = factor(probability2, ordered = TRUE),
                               label = NULL),
                 arrow = arrow(length = unit(0.1, "cm"))) +
    scale_y_discrete(labels=sapply(sort(pathway_df$probs), function(x) paste0(round(x, prob_digits), "%"))) +
    xlim(c(0, max(df$waiting_time)*1.1))
  
  if (!is.null(mutation_colors)) {
    mutation_colors <- c("#E0E0E0", mutation_colors)
    names(mutation_colors) <- c("Root", mutations)
    g <- g + scale_fill_manual(values = mutation_colors)
  }
  
  return(g)
}

##' @name plot_Theta
##' @title Plot a Mutual Hazard Network
##' @description This function plots a Mutual Hazard Network by separating the 
##' diagonal entries from the off-diagonal entries.
##' @param Theta An n-by-n matrix representing a Mutual Hazard Network.
##' @param mutations A list of mutation names, which must be unique values.
##' If no names are given, then the mutation IDs will be used.
##' @param full A boolean flag indicating whether the full MHN is shown (default: TRUE).
##' If not, then only mutations with non-zero off-diagonal entries will be shown.
##' @param sort_diag A boolean flag indicating whether the rows and columns of the MHN are
##' sorted based on the decreasing order of the diagonal entries (default: TRUE).
##' @param to_show A vector indicating which rows and columns of the MHN are shown (default: NULL).
##' @import ggplot2
##' @importFrom gridExtra grid.arrange
##' @importFrom reshape2 melt
##' @importFrom grDevices colors
##' @importFrom ggpubr annotate_figure
##' @export
plot_Theta <- function(Theta, mutations = NULL, full = TRUE, sort_diag = TRUE, to_show = NULL) {
  
  n <- nrow(Theta)
  
  if (is.null(mutations)) {
    mutations <- sapply(c(1:n), function(x) paste0("V", x))
  } else {
    if (length(mutations) != n) {
      stop("The number of mutations doesn't match matrix dimension. Please check again...")
    } else if (length(unique(mutations)) != n) {
      stop("Mutation names must be unique. Please check again...")
    }
  }
  
  dimnames(Theta) <- list(mutations, mutations)
  
  if (sort_diag) {
    idx <- order(diag(Theta), decreasing = TRUE)
  } else {
    idx <- c(1:n)
  }
  
  if (!is.null(to_show)) {
    idx <- to_show
  }
  
  if (full) {
    
    Theta_no_diag <- Theta
    diag(Theta_no_diag) <- 0
    Theta_diag <- matrix(diag(Theta)[idx], ncol = 1)
    rownames(Theta_diag) <- mutations[idx]
    colnames(Theta_diag) <- c("temp")
    
    g1 <- melt(Theta_diag) %>%
      ggplot(aes(x = Var2, y = Var1)) + 
      geom_raster(aes(fill=value)) +
      scale_fill_gradient(low = "white", high = colors()[77], 
                          labels = function(x) sprintf("%.1f", round(x, 1))) +
      scale_y_discrete(limits=rev) +
      theme(axis.text.x=element_blank(),
            axis.ticks=element_blank(),
            legend.title = element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
    
    g1 <- annotate_figure(g1, bottom = text_grob("Baseline rates"))
    
    g2 <- melt(Theta_no_diag[idx, idx]) %>%
      ggplot(aes(x = Var2, y = Var1)) + 
      geom_raster(aes(fill=value)) +
      scale_fill_gradient2(labels = function(x) sprintf("%.1f", round(x, 1))) +
      scale_y_discrete(limits=rev) + scale_x_discrete(position = "top") +
      theme(axis.text.x=element_text(angle=45, hjust = 0.2, vjust = 0.1),
            axis.ticks=element_blank(),
            legend.title = element_blank()) +
      xlab("ancestor") + ylab("descendant")
    
    g2 <- annotate_figure(g2, bottom = text_grob("Pattern of exclusivity and co-occurrence"))
    
    G <- grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2,6))
    
    
  } else {
    
    if (is.null(to_show)) {
      to_show <- sapply(c(1:n), function (i) any(Theta[i,-i] != 0) || any(Theta[-i,i] != 0))
      to_show_mat <- Theta[to_show, to_show]
      dimnames(to_show_mat) <- list(mutations[to_show], mutations[to_show])
      to_show_diag <- matrix(diag(Theta)[to_show], ncol = 1)
      rownames(to_show_diag) <- mutations[to_show]
      colnames(to_show_diag) <- c("temp")
      diag(to_show_mat) <- 0
      
      if (sort_diag) {
        to_show_idx <- order(to_show_diag, decreasing = TRUE)
        to_show_diag <- as.matrix(to_show_diag[to_show_idx,])
        
        to_show_no_diag <- melt(to_show_mat[to_show_idx, to_show_idx])
      } else {
        diag(to_show_mat) <- 0
        to_show_no_diag <- melt(to_show_mat)
      }
      
    } else {
      to_show_mat <- Theta[to_show, to_show]
      dimnames(to_show_mat) <- list(mutations[to_show], mutations[to_show])
      to_show_diag <- matrix(diag(Theta)[to_show], ncol = 1)
      rownames(to_show_diag) <- mutations[to_show]
      colnames(to_show_diag) <- c("temp")
      diag(to_show_mat) <- 0
      to_show_no_diag <- melt(to_show_mat)
    }
    
    g1 <- melt(to_show_diag) %>%
      ggplot(aes(x = Var2, y = Var1)) + 
      geom_raster(aes(fill=value)) +
      scale_fill_gradient(low = "white", high = colors()[77], 
                          labels = function(x) sprintf("%.1f", round(x, 1))) +
      scale_y_discrete(limits=rev) +
      theme(axis.text.x=element_blank(),
            axis.ticks=element_blank(),
            legend.title = element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
    
    g1 <- annotate_figure(g1, bottom = text_grob("Baseline rates"))
    
    g2 <- ggplot(to_show_no_diag, aes(x = Var2, y = Var1)) + 
      geom_raster(aes(fill=value)) + 
      geom_label(data = to_show_no_diag %>% mutate(text = ifelse(Var1 == Var2, as.character(Var2), NA)) %>% filter(!is.na(text)),
                 mapping = aes(label = text), 
                 size = 2.5,
                 label.padding = unit(0.1, "lines")) +
      scale_fill_gradient2(labels = function(x) sprintf("%.1f", round(x, 1))) +
      scale_y_discrete(limits=rev) + scale_x_discrete(position = "top") +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            legend.title = element_blank()) +
      xlab("ancestor") + ylab("descendant")
    
    # + theme(plot.background = element_rect(color = "white", fill = "lightgrey"),
    #         legend.background = element_rect(fill = "lightgrey"))
    
    g2 <- annotate_figure(g2, bottom = text_grob("Pattern of exclusivity and co-occurrence"))
    
    G <- grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(1,4))
  }
  
  return(G)
  
}