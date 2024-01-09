<div align="left">
  <img src="https://github.com/cbg-ethz/TreeMHN/blob/main/figures/TreeMHN_Logo.png", width="400px">
</div>
<p></p>

This R package implements the TreeMHN model for the joint inference of exclusivity patterns and recurrent trajectories from tumor mutation trees. ([bioRxiv preprint](https://doi.org/10.1101/2021.11.04.467347))

## Quick start

TreeMHN takes as input a set of independent tumor mutation trees containing a total number of n mutations. The format is a dataframe with five columns:

- `Patient_ID`: IDs of patients, unique for each patient;

- `Tree_ID`: IDs of mutation trees, unique within each patient;

- `Node_ID`: IDs of each node in the tree, including the root node (with ID "1"), unique for each node;

- `Mutation_ID`: IDs of each mutational event. The root node has a mutation ID of "0", and other mutation IDs can be duplicated in the tree to allow for parallel mutations;

- `Parent_ID`: IDs of the parent node ID. The root node has itself as parent (ID "1").

The output is an n-by-n matrix representing the Mutual Hazard Network. The diagonal entries of this matrix indicate how often each mutation will occur and fixate independent of the other mutations. The off-diagonal entries encode the exclusivity and co-occurrence patterns of mutations. Conditioned on the estimated matrix, we can compute the probabilities of different evolutionary trajectories or evaluate the most likely next mutational events given a tumor tree.

Please see [Demo.md](https://github.com/cbg-ethz/TreeMHN/blob/main/Example/Demo.md) for more details.

## Installation

For Mac users, please compile the package with g++ instead of clang. To do this, you need to first install gcc using [Homebrew](https://formulae.brew.sh/formula/gcc):

```
brew install gcc
```

Then, create `~/.R/Makevars` with entry

```
CXX=$(brew --prefix)/bin/g++-[INSTALLED VERSION]
```

For all users, install the `devtools` package in R and run

```
devtools::install_github("cbg-ethz/TreeMHN")
```

The installation typically takes around one minute to finish.


<div align="left">
  <img src="https://github.com/cbg-ethz/TreeMHN/blob/main/figures/Affiliations.png", width="800px">
</div>
<p></p>
