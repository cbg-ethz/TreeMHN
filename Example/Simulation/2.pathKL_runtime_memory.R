# Make sure you are in the "TreeMHN/Example/Simulation" folder
source("0.setup.R")

# Options for you to run this script from command line
option_list <- list( 
  make_option(c("-n", "--mutations"), action = "store", type = "integer",
              default = 10, help = "Total number of mutations [default %default]"),
  make_option(c("-N", "--samples"), action = "store", type = "integer",
              default = 500, help = "Total number of patient samples [default %default]"),
  make_option(c("-s", "--sparsity"), action = "store", type = "double",
              default = 0.5, help = "Network sparsity (proportion of zero entries) [default %default]"),
  make_option(c("--er"), action = "store", type = "double",
              default = 0.5, help = "Proportion of negative off-diagonal entries [default %default]"),
  make_option(c("--mc"), action = "store", type = "integer",
              default = 50, help = "A threshold on the maximum number of subtrees of a given tree, above which Monte Carlo sampling will be used [default %default]"),
  make_option(c("-m", "--mcsamples"), action = "store", type = "integer",
              default = 300, help = "Number of Monte Carlo samples to be drawn [default %default]"),
  make_option(c("-g", "--gamma"), action = "store", type = "double", default = 1.0,
              help = "Penalization level gamma [default %default]"),
  make_option(c("--norder"), action = "store", type = "integer", default = 4,
              help = "Length of the pathways [default %default]"),
  make_option(c("--iter"), action = "store", type = "integer", default = 500,
              help = "Number of EM iterations [default %default]"),
  make_option(c("--method"), default = "TreeMHN",
              help = "Inference method: \n
              1. TreeMHN_MLE [default];\n
              2. TreeMHN_MCEM;\n
              3. Baseline_MHN;\n
              4. Others."),
  make_option(c("--ncores"), action = "store", type = "integer", default = 100,
              help = "Number of cores to use [default %default]"),
  make_option(c("--version"), action = "store", type = "integer", default = 1,
              help = "Version of simulated data [default %default]"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false", 
              dest="verbose", help="Print little output")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

# Function to get one simulation run
run_once <- function(n = 10, N = 500, gamma = 1, s = 0.5, er = 0.5, mc = 50, 
                     m = 300, norder = 4, iter = 500, i = 1, method = "TreeMHN_MLE", 
                     v = 1) {
  
  # Generate trees
  set.seed((i + N)*n + v) # seed for the trees
  tree_obj <- generate_trees(n = n, N = N, sparsity = s, exclusive_ratio = er)
  true_Theta <- tree_obj$Theta # true MHN
  trees <- tree_obj$trees # extract trees
  true_p <- Theta_to_pathways(true_Theta)
  
  # Run methods
  switch(method,
         "TreeMHN_MLE" = {
           
           cat("TreeMHN (MLE) starts\n")
           start_time <- Sys.time()
           p <- profmem({Theta_TreeMHN_MLE <- learn_MHN(tree_obj, gamma = gamma, MC_threshold = 10000)})
           end_time <- Sys.time()
           
           res <- list(Theta_TreeMHN_MLE = Theta_TreeMHN_MLE,
                       Pathway_KL = KL_divergence(true_p, Theta_to_pathways(Theta_TreeMHN_MLE)),
                       Runtime = difftime(end_time, start_time, units = "secs"),
                       Memory = sum(p$bytes, na.rm = TRUE) / 1024^2)
           
           cat("Memory used: ", res$Memory, "MB\n")
           cat("TreeMHN (MLE) ends\n")
           
         },
         "TreeMHN_MCEM" = {
           
           cat("TreeMHN (MCEM) starts\n")
           start_time <- Sys.time()
           p <- profmem({Theta_TreeMHN_MCEM <- learn_MHN(tree_obj, gamma = gamma, use_EM = TRUE, 
                                                         MC_threshold = mc, M = m, iterations = iter)})
           end_time <- Sys.time()
           
           res <- list(Theta_TreeMHN_MCEM = Theta_TreeMHN_MCEM,
                       Pathway_KL = KL_divergence(true_p, Theta_to_pathways(Theta_TreeMHN_MCEM)),
                       Runtime = difftime(end_time, start_time, units = "secs"),
                       Memory = sum(p$bytes, na.rm = TRUE) / 1024^2)
           
           cat("Memory used: ", res$Memory, "MB\n")
           cat("TreeMHN (MCEM) ends\n")
           
         },
         "Baseline_MHN" = {
           
           cat("Computing pD...\n")
           pD <- genotypes_to_pD(trees)
           cat("Baseline MHN starts\n")
           start_time <- Sys.time()
           p <- profmem({Theta_baseline_MHN <- Learn.MHN(pD, lambda = gamma/N)})
           end_time <- Sys.time()
           
           res <- list(Theta_baseline_MHN = Theta_baseline_MHN,
                       Pathway_KL = KL_divergence(true_p, Theta_to_pathways(Theta_baseline_MHN)),
                       Runtime = difftime(end_time, start_time, units = "secs"),
                       Memory = sum(p$bytes, na.rm = TRUE) / 1024^2)
           
           cat("Memory used: ", res$Memory, "MB\n")
           cat("Baseline MHN ends\n")
           
         },
         "Others" = {
           # REVOLVER
           cat("REVOLVER starts\n")
           Pathway_KL_REVOLVER <- KL_divergence(true_p, 
                                                get_revolver_pathways(n, 
                                                                      tree_obj$tree_df, 
                                                                      n_order = norder)$probs)
           cat("REVOLVER ends\n")
           
           # HINTRA
           cat("HINTRA starts\n")
           Pathway_KL_HINTRA <- KL_divergence(true_p, 
                                              get_hintra_pathways(n, 
                                                                  tree_obj$tree_df, 
                                                                  n_order = norder)$probs)
           cat("HINTRA ends\n")
           
           res <- list(Pathway_KL_REVOLVER = Pathway_KL_REVOLVER, 
                       Pathway_KL_HINTRA = Pathway_KL_HINTRA)
           
         })
  
  return(res)
  
}

all_res <- mclapply(c(1:opt$ncores), 
                    function(i) run_once(n = opt$mutations, N = opt$samples, 
                                         gamma = opt$gamma, s = opt$sparsity, 
                                         er = opt$er, mc = opt$mc, 
                                         m = opt$mcsamples, norder = opt$norder,
                                         iter = opt$iter, i = i, method = opt$method, 
                                         v = opt$version), 
                    mc.cores = opt$ncores,
                    mc.preschedule = FALSE)

if (!file.exists("./pathKL_runtime_memory")) {
  dir.create("./pathKL_runtime_memory")
}
setwd("./pathKL_runtime_memory/")

dir_name <- sprintf("n%d_N%d_s%.1f_er%.1f_MC%d_M%d_gamma%.2f_norder%d_%d",
                   opt$mutations, opt$samples, opt$sparsity, opt$er, opt$mc, opt$mcsamples,
                   opt$gamma, opt$norder, opt$ncores)
if (!file.exists(dir_name)) {
  dir.create(dir_name)
}
setwd(dir_name)

fname <- sprintf("%s_v%d.RData", opt$method, opt$version)
save(all_res, file = fname)


