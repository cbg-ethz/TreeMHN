library(ggplot2)
library(ggpubr)
library(TreeMHN)
library(dplyr)
library(forcats)

###################### Precision - Recall curves ##################
plot_precision_recall <- function(n, N, s = 0.5, er = 0.5, 
                                  mc = 50, m = 300, gamma = 1, norder = 4, 
                                  ncores = 125, nsubsampling = 200, thr = 0.95,
                                  v = 1, mytitle = NULL, top_mutations = TRUE) {
  
  dir_name <- sprintf("n%d_N%d_s%.1f_er%.1f_MC%d_M%d_gamma%.2f_norder%d_cores%d_subsampling%d_threshold%.2f_v%d",
                      n, N, s, er, mc, m, gamma, norder, ncores, nsubsampling, thr, v)
  dir_name <- paste0("./structure_learning_w_SS/", dir_name)
  files <- list.files(path = dir_name,
                      full.names = TRUE,
                      recursive = FALSE)
  
  df <- data.frame(matrix(nrow = 0,ncol = 12))
  
  for (f in files) {
    load(f)
    true_Theta <- res$tree_obj$Theta
    if (top_mutations) { # only plot for the top half of mutations ranked by baseline rates
      idx <- order(diag(true_Theta),decreasing = TRUE)[c(1:(res$tree_obj$n /2))]
    } else {
      idx <- order(diag(true_Theta),decreasing = TRUE)
    }
    
    for (i in c(1:length(res$gamma_list))) {
      df <- rbind(df, c("TreeMHN",res$gamma_list[i],
                        compare_Theta(true_Theta[idx,idx],res$pred_Thetas[i,idx,idx])))
      df <- rbind(df, c("TreeMHN_w_SS",res$gamma_list[i],
                        compare_Theta(true_Theta[idx,idx],res$pred_Thetas_w_SS[i,idx,idx])))
      df <- rbind(df, c("MHN",res$gamma_list[i],
                        compare_Theta(true_Theta[idx,idx],res$pred_Thetas_genotype[i,idx,idx])))
      if (res$tree_obj$n < 20) {
        df <- rbind(df, c("MHN_w_SS",res$gamma_list[i],
                          compare_Theta(true_Theta[idx,idx],res$pred_Thetas_genotype_w_SS[i,idx,idx])))
      }
      
    }
  }
  
  df[,c(2:ncol(df))] <- lapply(df[,c(2:ncol(df))],as.numeric)
  
  colnames(df) <- c("Method","gamma","SHD","TP","FP","TN","FN","Precision","TPR","FPR_N","FPR_P","MSE")
  
  df_mean <- with(df, aggregate(. ~ Method + gamma, df, mean))
  
  if (res$tree_obj$n < 20) {
    myPalette <- c("#FF8000","#FFB266","#9933FF","#CC99FF")
  } else {
    myPalette <- c("#FF8000","#9933FF","#CC99FF")
  }
  
  random_precision <- 0.5 * (1 - s)
  g <- ggplot(df_mean) + xlim(c(0,1)) + ylim(c(0,1)) + xlab("Recall") + labs(title = mytitle) +
    geom_path(mapping = aes(x = TPR, y = Precision, group = Method, colour = Method),size = 1) +
    geom_segment(aes(x = 0, xend = 0.5, y = random_precision, yend = random_precision), color = "darkgrey", lty= "dashed") +
    geom_segment(aes(x = 0.5, xend = 0.5, y = 0, yend = random_precision), color = "darkgrey", lty= "dashed") +
    scale_color_manual(labels = c("Genotype MHN", "Genotype MHN (stability selection)",
                                  "TreeMHN","TreeMHN (stability selection)"),
                       values = myPalette) +
    theme(plot.title = element_text(size=10),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.text=element_text(size=12))
  return(g)
  
}

# Plot precision recall curves (only top half of mutations ranked by baseline rate)
g1 <- plot_precision_recall(n = 10, N = 100, mytitle = "n = 10, N = 100", top_mutations = TRUE)
g2 <- plot_precision_recall(n = 10, N = 200, mytitle = "n = 10, N = 200", top_mutations = TRUE)
g3 <- plot_precision_recall(n = 10, N = 300, mytitle = "n = 10, N = 300", top_mutations = TRUE)
g4 <- plot_precision_recall(n = 10, N = 500, mytitle = "n = 10, N = 500", top_mutations = TRUE)
g5 <- plot_precision_recall(n = 15, N = 100, mytitle = "n = 15, N = 100", top_mutations = TRUE)
g6 <- plot_precision_recall(n = 15, N = 200, mytitle = "n = 15, N = 200", top_mutations = TRUE)
g7 <- plot_precision_recall(n = 15, N = 300, mytitle = "n = 15, N = 300", top_mutations = TRUE)
g8 <- plot_precision_recall(n = 15, N = 500, mytitle = "n = 15, N = 500", top_mutations = TRUE)
g9 <- plot_precision_recall(n = 20, N = 100, mytitle = "n = 20, N = 100", top_mutations = TRUE)
g10 <- plot_precision_recall(n = 20, N = 200, mytitle = "n = 20, N = 200", top_mutations = TRUE)
g11 <- plot_precision_recall(n = 20, N = 300, mytitle = "n = 20, N = 300", top_mutations = TRUE)
g12 <- plot_precision_recall(n = 20, N = 500, mytitle = "n = 20, N = 500", top_mutations = TRUE)

fig <- ggarrange(g1, g2, g3, g4,
                 g5, g6, g7, g8,
                 g9, g10, g11, g12,
                 nrow = 3, ncol = 4, common.legend = TRUE, legend = "top")
annotate_figure(fig, left = "Precision", bottom = "Recall")

# Plot precision recall curves (all mutations)
g1 <- plot_precision_recall(n = 10, N = 100, mytitle = "n = 10, N = 100", top_mutations = FALSE)
g2 <- plot_precision_recall(n = 10, N = 200, mytitle = "n = 10, N = 200", top_mutations = FALSE)
g3 <- plot_precision_recall(n = 10, N = 300, mytitle = "n = 10, N = 300", top_mutations = FALSE)
g4 <- plot_precision_recall(n = 10, N = 500, mytitle = "n = 10, N = 500", top_mutations = FALSE)
g5 <- plot_precision_recall(n = 15, N = 100, mytitle = "n = 15, N = 100", top_mutations = FALSE)
g6 <- plot_precision_recall(n = 15, N = 200, mytitle = "n = 15, N = 200", top_mutations = FALSE)
g7 <- plot_precision_recall(n = 15, N = 300, mytitle = "n = 15, N = 300", top_mutations = FALSE)
g8 <- plot_precision_recall(n = 15, N = 500, mytitle = "n = 15, N = 500", top_mutations = FALSE)
g9 <- plot_precision_recall(n = 20, N = 100, mytitle = "n = 20, N = 100", top_mutations = FALSE)
g10 <- plot_precision_recall(n = 20, N = 200, mytitle = "n = 20, N = 200", top_mutations = FALSE)
g11 <- plot_precision_recall(n = 20, N = 300, mytitle = "n = 20, N = 300", top_mutations = FALSE)
g12 <- plot_precision_recall(n = 20, N = 500, mytitle = "n = 20, N = 500", top_mutations = FALSE)

fig <- ggarrange(g1, g2, g3, g4,
                 g5, g6, g7, g8,
                 g9, g10, g11, g12,
                 nrow = 3, ncol = 4, common.legend = TRUE, legend = "top")
annotate_figure(fig, left = "Precision", bottom = "Recall")


##################### Pathway KL/Runtime/Memory #######################

get_pathwayKL_runtime_memory <- function(method = NULL, n = 10, N = 500, s = 0.5, 
                                         er = 0.5, mc = 50, m = 300, gamma = 1, 
                                         norder = 4, ncores = 100) {
  
  df <- data.frame(matrix(nrow = 0, ncol = 4))
  
  dir_name <- sprintf("n%d_N%d_s%.1f_er%.1f_MC%d_M%d_gamma%.2f_norder%d_%d",
                      n, N, s, er, mc, m, gamma, norder, ncores)
  dir_name <- paste0("./pathKL_runtime_memory/", dir_name)
  files <- list.files(path = dir_name,
                      pattern = method,
                      full.names = TRUE,
                      recursive = FALSE)
  
  for (f in files) {
    load(f)
    if (method == "Others") {
      for (x in all_res) {
        df <- rbind(df, c(x$Pathway_KL_REVOLVER, 0, 0, "REVOLVER"))
        df <- rbind(df, c(x$Pathway_KL_HINTRA, 0, 0, "HINTRA"))
      }
    } else {
      df <- rbind(df, t(sapply(all_res, function(x) c(x$Pathway_KL, x$Runtime, x$Memory))))
    }
  }
  
  if (method != "Others") {
    df <- cbind(df, method)
  }

  df <- cbind(seq(1:nrow(df)), df, n, N)
  
  colnames(df) <- c("Sample", "Pathway_KL", "Runtime", "Memory", "Method", "n", "N")
  df$Pathway_KL <- as.numeric(df$Pathway_KL)
  df$Runtime <- as.numeric(df$Runtime)
  df$Memory <- as.numeric(df$Memory)
  return(df)
  
}

plot_pathway_KL <- function(n) {
  
  all_info_df <- rbind(get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 100),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 100),
                       get_pathwayKL_runtime_memory("Others", n, 100),
                       get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 200),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 200),
                       get_pathwayKL_runtime_memory("Others", n, 200),
                       get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 300),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 300),
                       get_pathwayKL_runtime_memory("Others", n, 300),
                       get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 500),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 500),
                       get_pathwayKL_runtime_memory("Others", n, 500))
  
  if (n <= 20) {
    all_info_df <- rbind(all_info_df,
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 100),
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 200),
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 300),
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 500))
  }
  
  # Take the minimum of "TreeMHN_MLE" and "TreeMHN_MCEM" to get "TreeMHN" pathway KL
  # (They are actually quite similar though...)
  temp <- all_info_df %>% 
    filter(Method %in% c("TreeMHN_MLE", "TreeMHN_MCEM")) %>% 
    group_by(Sample, n, N) %>% 
    summarise(Pathway_KL = min(Pathway_KL)) %>%
    mutate(Runtime = 0, Memory = 0, Method = "TreeMHN") %>%
    arrange(Sample, Pathway_KL, Runtime, Memory, Method, n, N)
  
  all_info_df <- bind_rows(all_info_df, temp) %>% select(!Sample)
  
  if (n <= 20) {
    g <- all_info_df %>%
      filter(Method %in% c("TreeMHN", 
                           "Baseline_MHN",
                           "REVOLVER",
                           "HINTRA"),
             n == n) %>%
      mutate(Method = fct_relevel(Method, 
                                  "TreeMHN", 
                                  "Baseline_MHN",
                                  "REVOLVER",
                                  "HINTRA")) %>%
      ggplot(mapping = aes(x = factor(N), y = Pathway_KL, fill = Method)) +
      geom_boxplot() + xlab("") + ylab("") + ylim(c(0, 5)) +
      scale_fill_manual(values = c("#CC99FF",
                                   "#FFB266",
                                   "#CCFF99",
                                   "#99CCFF"), 
                        labels = c("TreeMHN", 
                                   "Genotype MHN",
                                   "REVOLVER",
                                   "HINTRA")) +
      labs(title = paste0("n = ", n)) +
      theme(legend.text=element_text(size=12)) +
      ylab("Trajectory KL divergence") +
      xlab("Sample Size N")
  } else {
    g <- all_info_df %>%
      filter(Method %in% c("TreeMHN",
                           "REVOLVER",
                           "HINTRA"),
             n == n) %>%
      mutate(Method = fct_relevel(Method, 
                                  "TreeMHN",
                                  "REVOLVER",
                                  "HINTRA")) %>%
      ggplot(mapping = aes(x = factor(N), y = Pathway_KL, fill = Method)) +
      geom_boxplot() + xlab("") + ylab("") + ylim(c(0, 5)) +
      scale_fill_manual(values = c("#CC99FF",
                                   "#CCFF99",
                                   "#99CCFF"), 
                        labels = c("TreeMHN",
                                   "REVOLVER",
                                   "HINTRA")) +
      labs(title = paste0("n = ", n)) +
      theme(legend.text=element_text(size=12)) +
      ylab("Trajectory KL divergence") +
      xlab("Sample Size N")
  }
  
  return(g)
  
}

plot_runtime <- function(n) {
  
  all_info_df <- rbind(get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 100),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 100),
                       get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 200),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 200),
                       get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 300),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 300),
                       get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 500),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 500))
  
  if (n <= 20) {
    all_info_df <- rbind(all_info_df,
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 100),
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 200),
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 300),
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 500))
  }
  
  
  g <- all_info_df %>%
    filter(Method %in% c("TreeMHN_MLE",
                         "TreeMHN_MCEM",
                         "Baseline_MHN"),
           n == n) %>%
    mutate(Method = fct_relevel(Method, 
                                "TreeMHN_MLE",
                                "TreeMHN_MCEM",
                                "Baseline_MHN")) %>%
    ggplot(mapping = aes(x = factor(N), y = Runtime, fill = Method)) +
    geom_boxplot() + xlab("") + ylab("") +
    scale_fill_manual(values = c("#FF66B2", "#CC99FF", "#FFB266"),
                      labels = c("TreeMHN (MLE)",
                                 "TreeMHN (MC-EM)",
                                 "Genotype MHN")) +
    scale_y_log10(limits = c(1, 1e5)) +
    labs(title = paste0("n = ", n)) +
    theme(legend.text=element_text(size=12)) +
    ylab("Runtime (seconds)") +
    xlab("Sample Size N")
  
  return(g)
  
}


plot_memory <- function(n) {
  
  all_info_df <- rbind(get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 100),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 100),
                       get_pathwayKL_runtime_memory("Others", n, 100),
                       get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 200),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 200),
                       get_pathwayKL_runtime_memory("Others", n, 200),
                       get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 300),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 300),
                       get_pathwayKL_runtime_memory("Others", n, 300),
                       get_pathwayKL_runtime_memory("TreeMHN_MLE", n, 500),
                       get_pathwayKL_runtime_memory("TreeMHN_MCEM", n, 500),
                       get_pathwayKL_runtime_memory("Others", n, 500))
  
  if (n <= 20) {
    all_info_df <- rbind(all_info_df,
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 100),
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 200),
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 300),
                         get_pathwayKL_runtime_memory("Baseline_MHN", n, 500))
  }
  
  g <- all_info_df %>%
    filter(Method %in% c("TreeMHN_MLE",
                         "TreeMHN_MCEM",
                         "Baseline_MHN"),
           n == n) %>%
    mutate(Method = fct_relevel(Method, 
                                "TreeMHN_MLE",
                                "TreeMHN_MCEM",
                                "Baseline_MHN")) %>%
    ggplot(mapping = aes(x = factor(N), y = Memory, fill = Method)) +
    geom_boxplot() + xlab("") + ylab("") +
    scale_fill_manual(values = c("#FF66B2", "#CC99FF", "#FFB266"),
                      labels = c("TreeMHN (MLE)",
                                 "TreeMHN (MC-EM)",
                                 "Genotype MHN")) +
    scale_y_log10(limits = c(1, 1e7)) +
    labs(title = paste0("n = ", n)) +
    theme(legend.text=element_text(size=12)) +
    ylab("Memory (MiB)") +
    xlab("Sample Size N")
  
  return(g)
  
}


# Plot pathway KL
g1 <- plot_pathway_KL(10)
g2 <- plot_pathway_KL(15)
g3 <- plot_pathway_KL(20)
g4 <- plot_pathway_KL(30)

g <- ggarrange(g1 + rremove("xlab") + rremove("ylab"),
               g2 + rremove("xlab") + rremove("ylab"), 
               g3 + rremove("xlab") + rremove("ylab"), 
               g4 + rremove("xlab") + rremove("ylab"), 
               nrow = 1, ncol = 4, common.legend = TRUE, legend = "top")
annotate_figure(g, left = text_grob("Trajectory KL divergence", size= 12, rot = 90), 
                bottom = text_grob("Sample size N",size=12))

# Plot runtime
g1 <- plot_runtime(10)
g2 <- plot_runtime(15)
g3 <- plot_runtime(20)
g4 <- plot_runtime(30)
fig <- ggarrange(g1 + rremove("xlab") + rremove("ylab"), 
                 g2 + rremove("xlab") + rremove("ylab"), 
                 g3 + rremove("xlab") + rremove("ylab"), 
                 g4 + rremove("xlab") + rremove("ylab"), 
                 ncol = 4, nrow = 1, common.legend = TRUE, legend = "top")
G1 <- annotate_figure(fig, left = "Runtime (seconds)")

# Plot memory
g1 <- plot_memory(10)
g2 <- plot_memory(15)
g3 <- plot_memory(20)
g4 <- plot_memory(30)
fig <- ggarrange(g1 + rremove("xlab") + rremove("ylab"), 
                 g2 + rremove("xlab") + rremove("ylab"), 
                 g3 + rremove("xlab") + rremove("ylab"), 
                 g4 + rremove("xlab") + rremove("ylab"), 
                 ncol = 4, nrow = 1, common.legend = TRUE, legend = "none")
G2 <- annotate_figure(fig, left = "Memory (MiB)", bottom = "Sample size N")

# Putting together
ggarrange(G1, G2, ncol = 1, nrow = 2)
