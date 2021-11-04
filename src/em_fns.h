#ifndef EM_FUNCTIONS_H
#define EM_FUNCTIONS_H

#include <RcppArmadillo.h>
#include "utils.h"
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <string>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double get_exit_lambda(int n, IntegerVector genotype, arma::mat Theta, const List &node_labels);
arma::sp_mat build_tr_mat(int n, arma::mat Theta, const IntegerMatrix &genotypes, const List &node_labels);
double compute_obs_ll(arma::sp_mat tr_mat, double lambda_s);
int find_node_idx(std::vector<int> node, const List &node_labels);
arma::sp_mat dS_dlambda(std::vector<int> node, bool in_tree, const IntegerMatrix &genotypes, const List &node_labels);
arma::sp_mat build_AS(arma::sp_mat tr_mat, arma::sp_mat d_tr_mat);
long double compute_dp_dlambda(arma::sp_mat tr_mat, std::vector<int> node,
                               bool in_tree, double lambda_s,const IntegerMatrix &genotypes, const List &node_labels);
void tree_E_step_(const arma::mat Theta, NumericVector &time_diffs, double prob_p, int n,
                  const IntegerVector &nodes, const List &children,
                  const LogicalVector &in_tree, arma::sp_mat tr_mat, double lambda_s,
                  const IntegerMatrix &genotypes, const List &node_labels,
                  std::vector<int> pathway = {}, int current_pos = 0);
NumericVector tree_E_step(const arma::mat Theta, int n, double lambda_s,
                          const IntegerVector &nodes, const List &children, const LogicalVector &in_tree,
                          const IntegerMatrix &genotypes, const List &node_labels);
double tree_imp_samp_once(const NumericVector Theta, const IntegerVector &nodes,
                          const List &children, const LogicalVector &in_tree,
                          int n, NumericVector &time_diffs, double t_s,
                          IntegerVector pathway = {}, int current_pos = 0, double t_pa = 0);
NumericVector tree_importance_sampling(NumericVector Theta, const IntegerVector &nodes,
                                       const List &children, const LogicalVector &in_tree,
                                       int n, int M, double lambda_s);


#endif //EM_FUNCTIONS_H
