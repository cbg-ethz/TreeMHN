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

double get_exit_lambda(const int &n, const IntegerVector &genotype, const arma::mat &Theta, const List &node_labels);
arma::sp_mat build_tr_mat(const int &n, const arma::mat &Theta, const IntegerMatrix &genotypes, 
                          const List &node_labels, const double &lambda_s);
double compute_obs_ll(const arma::sp_mat &tr_mat, const double &lambda_s);
int find_node_idx(const std::vector<int> &node, const List &node_labels);
arma::sp_mat dS_dlambda(std::vector<int> node, const bool &in_tree, const IntegerMatrix &genotypes, const List &node_labels);
arma::sp_mat build_AS(const arma::sp_mat &tr_mat, const arma::sp_mat &d_tr_mat);
long double compute_dp_dlambda(const arma::sp_mat &tr_mat, const std::vector<int> &node, 
                               const bool &in_tree, const double &lambda_s,
                               const IntegerMatrix &genotypes, const List &node_labels);
void tree_E_step_(const arma::mat &Theta, NumericVector &time_diffs, 
                  const double &prob_p, const int &n, const IntegerVector &nodes, 
                  const List &children, const LogicalVector &in_tree, const arma::sp_mat &tr_mat, 
                  const double &lambda_s, const IntegerMatrix &genotypes, const List &node_labels,
                  std::vector<int> pathway = {}, int current_pos = 0);
NumericVector tree_E_step(const arma::mat &Theta, const int &n, const double &lambda_s,
                          const IntegerVector &nodes, const List &children, const LogicalVector &in_tree,
                          double &obs_ll, const IntegerMatrix &genotypes, const List &node_labels);
double tree_imp_samp_once(const arma::mat &Theta, const IntegerVector &nodes,
                          const List &children, const LogicalVector &in_tree,
                          int n, NumericVector &time_diffs, const double &t_s,
                          IntegerVector pathway = {}, int current_pos = 0, double t_pa = 0);
NumericVector tree_importance_sampling(const arma::mat &Theta, const IntegerVector &nodes,
                                       const List &children, const LogicalVector &in_tree,
                                       const int &n, const int &M, const double &lambda_s,
                                       double &obs_ll);
void update_timed_trees(const int &n, const int &N, List &timed_trees, 
                        const arma::mat &Theta, const double &lambda_s, const int &M,
                        NumericVector &log_prob_vec, const List &comp_geno_vec, 
                        const List &node_labels_vec, const int &nr_exact);


#endif //EM_FUNCTIONS_H
