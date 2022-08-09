############ Load Packages #################

library(TreeMHN)
library(devtools)
library(parallel)
library(profmem)
library(optparse)

############ Functions for MHNs (Schill et al., 2020) #################

## Original code for genotype MHNs 
## Check if downloaded:
if (!file.exists("MHN-master")) {
  download.file(url = "https://github.com/RudiSchill/MHN/archive/refs/heads/master.zip",
                destfile = "MHN-master.zip")
  unzip(zipfile = "MHN-master.zip")
}

setwd("./MHN-master/")
source("UtilityFunctions.R")
source("ModelConstruction.R")
source("Likelihood.R") # may need "CC = /usr/local/opt/llvm/bin/clang -fopenmp" in your "~/.R/Makevars"
source("RegularizedOptimization.R")
setwd(dir = "..")


## The following two functions are modified to account for the case that
## no observed genotypes are wild-type genotypes (all zeros)
# Regularized Score
Score.Reg <- function(Theta, pD, lambda, to_mask = integer(0)){
  
  #Reshape parameters as matrix, after internal handling by BFGS as vectors.
  n <- sqrt(length(Theta))
  Theta <- matrix(Theta,nrow=n,ncol=n)
  
  if (length(to_mask) != 0) {
    Theta[to_mask] <- 0
  }
  
  # Score(Theta,pD) - lambda*L1(Theta)
  # Condition on non-empty genotypes
  p_empty <- 1 / (1 + sum(exp(diag(Theta)))) # equivalent: p_empty <- Generate.pTh(Theta)[1]
  Score(Theta,pD) - lambda*L1(Theta) - log(1 - p_empty)
}

# Regularized Gradient
Grad.Reg  <- function(Theta, pD, lambda, to_mask = integer(0)){
  n <- sqrt(length(Theta))
  Theta <- matrix(Theta,nrow=n,ncol=n)
  if (length(to_mask) != 0) {
    Theta[to_mask] <- 0
  }
  
  # Grad(Theta,pD) - lambda*L1_(Theta)
  # Condition on non-empty genotypes
  p_empty <- 1 / (1 + sum(exp(diag(Theta)))) # equivalent: p_empty <- Generate.pTh(Theta)[1]
  gd <- Grad(Theta,pD) - lambda*L1_(Theta) - diag(p_empty^2 / (1 - p_empty) * exp(diag(Theta)))
  if (length(to_mask) != 0) {
    gd[to_mask] <- 0
  }
  return(gd)
  
}

# This function is modified to add bound constraints on the optimization to avoid numerical errors!
Learn.MHN <- function(pD, init=NULL, lambda=0 ,maxit=5000, trace=0, reltol=1e-07, round=TRUE, to_mask = integer(0)){
  n <- log(length(pD), base=2)
  
  #Initialize the parameters from the independence model
  if(is.null(init)){
    init <- Learn.Indep(pD)
  }
  
  for (i in 1:n) {
    if (init[i,i] == Inf) {
      init[i,i] <- 100
    } else if (init[i,i] == -Inf) {
      init[i,i] <- -100
    }
  }
  
  opt <- optim(init, fn=Score.Reg, gr=Grad.Reg, pD, lambda, to_mask,
               method = "L-BFGS-B", lower = -20, upper = 20,
               control=list(fnscale=-1,trace=trace,maxit=maxit,factr = 1e11))
  
  Theta <- matrix(opt$par,nrow=n,ncol=n)
  
  if(round){
    Theta <- round(Theta,2)
  }
  
  return(Theta)
}

# This function is to convert subclonal genotypes to input vectors for genotype MHN
genotypes_to_pD <- function(trees) {
  N <- length(trees)
  pD <- 0
  for (i in c(1:N)) {
    pD <- pD + Data.to.pD(trees[[i]]$genotypes)
  }
  pD <- pD / N
  pD[1] <- 0
  pD <- pD / sum(pD)
  return(pD)
}
