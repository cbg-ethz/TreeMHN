# Make sure you are in the "TreeMHN/Example/Simulation" folder
source("0.setup.R")

# Functions for stability selection
subsample_once <- function(subsample_size, tree_obj, gamma) {
  subsample_patients <- sample(tree_obj$N_patients, subsample_size)
  subsample_trees <- list()
  subsample_weights <- c()
  for (i in c(1:tree_obj$N)) {
    tree <- tree_obj$trees[[i]]
    if (tree$patient_ID %in% subsample_patients) {
      subsample_trees <- append(subsample_trees, list(tree))
      subsample_weights <- c(subsample_weights, tree_obj$weights[i])
    }
  }
  tree_obj$trees <- subsample_trees
  tree_obj$N <- length(subsample_trees)
  tree_obj$N_patients <- subsample_size
  tree_obj$weights <- subsample_weights
  
  pred_Theta <- learn_MHN(tree_obj, gamma)
  pD <- genotypes_to_pD(subsample_trees)
  pred_Theta_genotype <- Learn.MHN(pD, lambda = gamma/subsample_size)
  return(list(pred_Theta, pred_Theta_genotype))
}

get_mask <- function(n, SS_res, idx, thr = 0.95) {
  to_mask <- sapply(SS_res, function (x) as.vector(x[[idx]]))
  to_mask[abs(to_mask) > 1e-3] <- 1
  to_mask[to_mask != 1] <- 0
  to_mask <- which(matrix(apply(to_mask,1,mean) > thr, nrow = n) == 0)
  return(to_mask)
}

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
  make_option(c("--norder"), action = "store", type = "integer", default = 4,
              help = "Length of the pathways [default %default]"),
  make_option(c("--iter"), action = "store", type = "integer", default = 500,
              help = "Number of EM iterations [default %default]"),
  make_option(c("-i", "--run"), action = "store", type = "integer", default = 1,
              help = "Simulation run ID [default %default]"),
  make_option(c("-g", "--gamma"), action = "store", type = "double", default = 1.0,
              help = "Penalization level gamma for the final output [default %default]"),
  make_option(c("--gamma_list"), default = "c(0.05, 0.1, 0.5, 1, 1.5, 2, 2.5, 3)",
              help = "List of penalty levels to test [default %default]"),
  make_option(c("--SS_threshold"), action = "store", type = "double", default = 0.95,
              help = "Threshold for stability selection [default %default]"),
  make_option(c("--ncores"), action = "store", type = "integer", default = 100,
              help = "Number of cores to use [default %default]"),
  make_option(c("--nsubsampling"), action = "store", type = "integer", default = 200,
              help = "Number of times to draw subsamples in stability selection [default %default]"),
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
run_once <- function(n = 10, N = 500, final_gamma = 1, s = 0.5, er = 0.5, mc = 50,
                     m = 300, norder = 4, iter = 500, i = 1, nsubsampling = 200,
                     v = 1, gamma_list = "c(0.05, 0.1, 0.5, 1, 1.5, 2, 2.5, 3)",
                     ncores = 100, thr = 0.95) {

  # Generate trees
  set.seed((i + N)*n + v) # seed for the trees
  tree_obj <- generate_trees(n = n, N = N, sparsity = s, exclusive_ratio = er)
  true_Theta <- tree_obj$Theta # true MHN
  trees <- tree_obj$trees # extract trees
  gamma_list <- eval(parse(text = gamma_list)) # penalty levels
  nr_gammas <- length(gamma_list) # number of penalty levels
  subsample_size <- floor(N/2)

  # list to store items
  res <- list()
  res$tree_obj <- tree_obj
  res$pD <- genotypes_to_pD(trees)
  res$gamma_list <- gamma_list

  # run through all penalty values
  pred_Thetas <- array(0, dim = c(nr_gammas, n, n))
  pred_Thetas_genotype <- array(0, dim = c(nr_gammas, n, n))
  pred_Thetas_w_SS <- array(0, dim = c(nr_gammas, n, n))
  pred_Thetas_genotype_w_SS <- array(0, dim = c(nr_gammas, n, n))

  for (i in c(1:nr_gammas)) {

    gamma <- res$gamma_list[i]

    ## without stability selection
    pred_Thetas[i,,] <- learn_MHN(tree_obj, gamma = gamma) #TreeMHN wo SS
    pred_Thetas_genotype[i,,] <- Learn.MHN(res$pD, lambda = gamma/N) #MHN wo SS

    ## with stability selection
    SS_res <- mclapply(c(1:nsubsampling),
                       function(i) subsample_once(subsample_size, tree_obj, gamma),
                       mc.cores = ncores,
                       mc.preschedule = FALSE)

    TreeMHN_to_mask <- get_mask(n, SS_res, 1, thr)
    pred_Thetas_w_SS[i,,] <- learn_MHN(tree_obj, gamma = final_gamma, to_mask = TreeMHN_to_mask)
    MHN_to_mask <- get_mask(n, SS_res, 2, thr)
    pred_Thetas_genotype_w_SS[i,,] <- Learn.MHN(res$pD, lambda = final_gamma/N, to_mask = MHN_to_mask)

  }

  res$pred_Thetas <- pred_Thetas
  res$pred_Thetas_genotype <- pred_Thetas_genotype
  res$pred_Thetas_w_SS <- pred_Thetas_w_SS
  res$pred_Thetas_genotype_w_SS <- pred_Thetas_genotype_w_SS

  return(res)

}

res <- run_once(n = opt$mutations, N = opt$samples,
                final_gamma = opt$gamma, s = opt$sparsity,
                er = opt$er, mc = opt$mc,
                m = opt$mcsamples, norder = opt$norder,
                iter = opt$iter, i = opt$run,
                nsubsampling = opt$nsubsampling, v = opt$version,
                gamma_list = opt$gamma_list, ncores = opt$ncores)

if (!file.exists("./structure_learning_w_SS")) {
  dir.create("./structure_learning_w_SS")
}
setwd("./structure_learning_w_SS/")

dir_name <- sprintf("n%d_N%d_s%.1f_er%.1f_MC%d_M%d_gamma%.2f_norder%d_cores%d_subsampling%d_threshold%.2f_v%d",
                    opt$mutations, opt$samples, opt$sparsity, opt$er, opt$mc, opt$mcsamples,
                    opt$gamma, opt$norder, opt$ncores, opt$nsubsampling, opt$SS_threshold, opt$version)

if (!file.exists(dir_name)) {
  dir.create(dir_name)
}
setwd(dir_name)

fname <- sprintf("run%d.RData", opt$run)
save(res, file = fname)