This document demonstrates how to use the `TreeMHN` package. The package
is created based on the article “Joint inference of repeated
evolutionary trajectories and patterns of clonal exclusivity or
co-occurrence from tumor mutation trees”.

# 1 Load required packages

``` r
library(TreeMHN)
library(parallel)
library(ggplot2)
```

# 2 Simulated data

## 2.1 Generate trees

The function `generate_trees` can generate a random Mutual Hazard
Network *Θ* and a set of mutation trees from *Θ* according to the tree
generating process introduced in the paper.

``` r
set.seed(6666)
n <- 10 # number of events
N <- 500 # number of samples
lambda_s <- 1 # sampling event rate
gamma <- 0.5 # penalty parameter
sparsity <- 0.5 # sparsity of the MHN
exclusive_ratio <- 0.5 # proportion of inhibiting edges

# This function will generate a random MHN along with a collection of mutation trees
tree_obj <- generate_trees(n = n, N = N, lambda_s = lambda_s, sparsity = sparsity, 
                           exclusive_ratio = exclusive_ratio)
true_Theta <- tree_obj$Theta # true MHN
trees <- tree_obj$trees # extract trees
```

We can plot one of the trees using the `plot_tree_list` function:

``` r
tree <- trees[[66]]
plot_tree_list(tree, tree_obj$mutations)
```

![](Demo_files/figure-gfm/display_one_tree-1.png)<!-- -->

## 2.2 Learn MHN from the generated trees

The function `learn_MHN` takes a `TreeMHN` object and learns an MHN *Θ̂*.

-   Without stability selection:

``` r
pred_Theta <- learn_MHN(tree_obj, gamma = gamma, lambda_s = lambda_s, verbose = TRUE)
```

    ## Initializing Theta...
    ## Checking whether MCEM is needed...
    ## Running MLE...
    ## iter   10 value 3534.182513
    ## iter   20 value 3520.105848
    ## iter   30 value 3515.849267
    ## iter   40 value 3514.837437
    ## iter   50 value 3514.370185
    ## final  value 3514.281048 
    ## converged

-   With stability selection ([Meinshausen and
    Bühlmann (2010)](https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9868.2010.00740.x)):

``` r
# Function to subsample half of the trees and learn one MHN
subsample_size <- floor(N/2)
subsample_once <- function(subsample_size, tree_obj, gamma) {
  subsample_trees <- sample(tree_obj$trees, subsample_size)
  tree_obj$trees <- subsample_trees
  tree_obj$N <- subsample_size
  tree_obj$N_patients <- subsample_size
  tree_obj$weights <- rep(1, subsample_size)
  pred_Theta <- learn_MHN(tree_obj, gamma)
  return(pred_Theta)
}

# Function to get a vector of non-selected (masked) elements 
get_mask <- function(n, SS_res, thr = 0.95) {
  to_mask <- sapply(SS_res, function (x) as.vector(x))
  to_mask[abs(to_mask) > 1e-2] <- 1
  to_mask[to_mask != 1] <- 0
  to_mask <- which(matrix(apply(to_mask,1,mean) > thr,nrow = n) == 0)
  return(to_mask)
}
```

It is recommended to run the stability selection procedure using the
`parallel` package.

``` r
SS_res <- mclapply(c(1:1000), 
                   function(i) subsample_once(subsample_size, tree_obj, gamma), 
                   mc.cores = detectCores())
TreeMHN_to_mask <- get_mask(n, SS_res, 0.99)
pred_Theta_w_SS <- learn_MHN(tree_obj, gamma = gamma, to_mask = TreeMHN_to_mask)
```

## 2.3 Performance assessment

We first plot the true MHN and the two estimated MHNs by ordering the
entries based on the true baseline rates. At a regularization level
*γ* = 0.5, most entries at the top left corner are recovered. Without
stability selection, there are many false non-zero entries due to
overfitting, including some very small values at the top left corner.
With stability selection, we see that those entries are removed, along
with some true positive entries at the lower left corner.

``` r
idx <- order(diag(true_Theta),decreasing = TRUE)
top_idx <- idx[c(1:(n/2))]
```

![](Demo_files/figure-gfm/theta_plots2-1.png)<!-- -->

We can compute the precision and recall (= true positive rate) based on
the off-diagonal differences using function `compare_Theta`. (See the
Supplementary Material for more details.)

| *Θ* (left) *Θ̂* (top) | *i*  *j* | *i* → *j* | *i* ⊣ *j* |
|----------------------|----------|-----------|-----------|
| *i*  *j*             | TN       | FP        | FP        |
| *i* → *j*            | FN       | TP        | FP        |
| *i* ⊣ *j*            | FN       | FP        | TP        |

-   Without stability selection:

``` r
compare_Theta(true_Theta, pred_Theta)
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     34.00     33.00     23.00     23.00     11.00      0.59      0.73      0.51 
    ##     FPR_P       MSE 
    ##      0.51      0.77

-   With stability selection:

``` r
compare_Theta(true_Theta, pred_Theta_w_SS)
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     30.00     16.00      1.00     44.00     29.00      0.94      0.36      0.02 
    ##     FPR_P       MSE 
    ##      0.02      0.89

If we focus on the first half of the events with higher baseline rates,
we can see an increase in recall/TPR.

-   Without stability selection:

``` r
compare_Theta(true_Theta[top_idx,top_idx], pred_Theta[top_idx,top_idx])
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##      7.00     10.00      7.00      3.00      0.00      0.59      1.00      0.70 
    ##     FPR_P       MSE 
    ##      0.70      0.11

-   With stability selection:

``` r
compare_Theta(true_Theta[top_idx,top_idx], pred_Theta_w_SS[top_idx,top_idx])
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##      1.00      9.00      0.00     10.00      1.00      1.00      0.90      0.00 
    ##     FPR_P       MSE 
    ##      0.00      0.12

# 3 Real data

## 3.1 Input dataset

Here we use the tree dataset from [Morita et
al. (2020)](https://www.nature.com/articles/s41467-020-19119-8).

``` r
# note that we summarize the mutations at the gene level
load("AML_tree_obj.RData")
plot_tree_df(AML$tree_df[AML$tree_df$Tree_ID == match("AML-38_AML-38-001", AML$tree_labels),],
             AML$mutations, "AML-38_AML-38-001")
```

![](Demo_files/figure-gfm/read_AML-1.png)<!-- -->

To use another dataset, please make sure it is in dataframe format with
five columns:

-   `Patient_ID`: IDs of patients, unique for each patient;

-   `Tree_ID`: IDs of mutation trees, unique within each patient;

-   `Node_ID`: IDs of each node in the tree, including the root node
    (with ID “1”), unique for each node;

-   `Mutation_ID`: IDs of each mutational event. The root node has a
    mutation ID of “0”, and other mutation IDs can be duplicated in the
    tree to allow for parallel mutations;

-   `Parent_ID`: IDs of the parent node ID. The root node has itself as
    parent (ID “1”).

For example,

``` r
head(AML$tree_df)
```

    ##   Patient_ID Tree_ID Node_ID Mutation_ID Parent_ID
    ## 1          1       1       1           0         1
    ## 2          1       1       2           4         1
    ## 3          1       1       3           2         2
    ## 4          1       1       4           1         3
    ## 5          1       1       5           5         3
    ## 6          1       1       6           6         3

To convert a dataframe to an `TreeMHN` object, use the `input_tree_df`
function. For example,

``` r
# not run
# input_tree_df(n = AML$n, tree_df = AML$tree_df, mutations = AML$mutations, tree_labels = AML$tree_labels)
```

## 3.2 Learn the MHN

To ensure enough precision, we run stability selection with *γ* = 0.1
and a threshold of 95% and obtain a vector of non-selected elements over
1000 subsamples. Again, we recommend to run the code using the
`parallel` package on a cluster.

``` r
RNGkind("L'Ecuyer-CMRG")
gamma <- 0.1
subsample_size <- floor(AML$N / 2)
SS_res <- mclapply(c(1:1000), 
                   function(i) subsample_once(subsample_size, AML, gamma), 
                   mc.cores = detectCores(), 
                   mc.set.seed = TRUE)
to_mask <- get_mask(AML$n, SS_res, 0.95)
```

Then we refit the model by masking the non-selected elements.

``` r
AML_Theta <- learn_MHN(AML, gamma = gamma, to_mask = to_mask)
```

``` r
save(AML_Theta, file = "AML_Theta.RData")
```

## 3.3 Plot the network

Now we can plot the learned MHN. The diagonal entries represent the
baseline rates of mutations to occur, independent of other events. The
off-diagonal entries represent the interactions between mutations.

``` r
plot_Theta(AML_Theta, AML$mutations)
```

![](Demo_files/figure-gfm/plot_AML_Theta_full-1.png)<!-- -->

    ## TableGrob (1 x 2) "arrange": 2 grobs
    ##   z     cells    name           grob
    ## 1 1 (1-1,1-1) arrange gtable[layout]
    ## 2 2 (1-1,2-2) arrange gtable[layout]

If we focus on the entries with non-zero off-diagonal entries, we get:

``` r
plot_Theta(AML_Theta, AML$mutations, full = FALSE)
```

![](Demo_files/figure-gfm/plot_AML_Theta-1.png)<!-- -->

Given the estimated MHN, we can plot the most probable evolutionary
trajectories of a given length using the `plot_pathways` function:

``` r
plot_pathways(AML_Theta, mutations = AML$mutations, top_M = 10, n_order = 4)
```

![](Demo_files/figure-gfm/plot_AML_pathways-1.png)<!-- -->

Given a particular tree and the estimated MHN, we can also find the next
most probable mutational events using the `plot_next_mutations`
function.

``` r
tree_df <- AML$tree_df[AML$tree_df$Tree_ID == match("AML-09_AML-09-001", AML$tree_labels),]
plot_next_mutations(AML$n, 
                    tree_df,
                    AML_Theta,
                    AML$mutations,
                    "AML-09_AML-09-001",
                    8)
```

    ## Top 8 most probable mutational events that will happen next:
    ## The next most probable node: Root->NPM1->KRAS->NRAS 
    ## Probability: 16.461 %
    ## The next most probable node: Root->NPM1->PTPN11 
    ## Probability: 10.353 %
    ## The next most probable node: Root->NPM1->KRAS->PTPN11 
    ## Probability: 10.353 %
    ## The next most probable node: Root->NPM1->KRAS->FLT3 
    ## Probability: 9.202 %
    ## The next most probable node: Root->NPM1->NRAS 
    ## Probability: 8.693 %
    ## The next most probable node: Root->NPM1->FLT3->WT1 
    ## Probability: 7.185 %
    ## The next most probable node: Root->NPM1->FLT3->KRAS 
    ## Probability: 5.511 %
    ## The next most probable node: Root->DNMT3A 
    ## Probability: 2.386 %

![](Demo_files/figure-gfm/next_mutations-1.png)<!-- -->
