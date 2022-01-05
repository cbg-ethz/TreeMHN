# TreeMHN

This package implements the TreeMHN model for the joint inference of repeated evolutionary trajectories and patterns of clonal exclusivity or co-occurrence from tumor mutation trees. ([bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2021.11.04.467347v1))

## Installation

In R, install the `devtools` package and run

```
devtools::install_github("cbg-ethz/TreeMHN")
```

For Mac users, please compile the package with g++ instead of clang by installing gcc using [Homebrew](https://formulae.brew.sh/formula/gcc)

```
brew install gcc
```

and creating `~/.R/Makevars` with entry

```
CXX=g++-11
```


## Example

Please see `Demo.md` for more details.
