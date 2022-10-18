<div align="left">
  <img src="https://github.com/cbg-ethz/TreeMHN/blob/main/figures/TreeMHN_Logo.png", width="400px">
</div>
<p></p>

This package implements the TreeMHN model for the Joint inference of exclusivity patterns and recurrent trajectories from tumor mutation trees. ([bioRxiv preprint](https://doi.org/10.1101/2021.11.04.467347))

## Installation

For Mac users, please compile the package with g++ instead of clang. To do this, you need to first install gcc using [Homebrew](https://formulae.brew.sh/formula/gcc):

```
brew install gcc
```

Then, create `~/.R/Makevars` with entry

```
CXX=$(brew --prefix)/bin/g++-[INSTALLED VERSION]
```

In R, install the `devtools` package and run

```
devtools::install_github("cbg-ethz/TreeMHN")
```

## Example

Please see `Demo.md` for more details.
