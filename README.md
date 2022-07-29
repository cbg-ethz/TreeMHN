# TreeMHN

This package implements the TreeMHN model for the Joint inference of exclusivity patterns and recurrent trajectories from tumor mutation trees. ([bioRxiv preprint](https://doi.org/10.1101/2021.11.04.467347))

## Installation

For Mac users, please compile the package with g++ instead of clang by installing gcc using [Homebrew](https://formulae.brew.sh/formula/gcc)

```
brew install gcc
```

and creating `~/.R/Makevars` with entry

```
CXX=g++-11
```

In R, install the `devtools` package and run

```
devtools::install_github("cbg-ethz/TreeMHN")
```

## Example

Please see `Demo.md` for more details.
