% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assessment.R
\name{trees_to_revolver_W}
\alias{trees_to_revolver_W}
\title{trees_to_revolver_W(n, tree_df)
This function computes the row-normalized W matrix of the REVOLVER algorithm}
\usage{
trees_to_revolver_W(n, tree_df)
}
\arguments{
\item{n}{Number of mutational events}

\item{tree_df}{A data frame with the following columns:
\itemize{
\item Tree_ID: IDs of mutation trees, unique for each patient
\item Node_ID: IDs of each node in the tree, including the root node (with ID "1"), unique for each node
\item Mutation_ID: IDs of each mutational event, the root node has a mutation ID of "0",
other mutation IDs can be duplicated in the tree to allow for parallel mutations
\item Parent_ID: IDs of the parent node ID. The root node has itself as parent (ID "1").
}}
}
\value{
Row-normalized W matrix of the REVOLVER algorithm.
}
\description{
trees_to_revolver_W(n, tree_df)
This function computes the row-normalized W matrix of the REVOLVER algorithm
}
\author{
Xiang Ge Luo
}
