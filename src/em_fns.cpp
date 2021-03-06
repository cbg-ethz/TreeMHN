#include <RcppArmadillo.h>
#include "utils.h"
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <string>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double get_exit_lambda(const int &n, const IntegerVector &genotype, const arma::mat &Theta, const List &node_labels) {

  if (sum(genotype) == 0) {
    return sum(exp(Theta.diag()));
  } else {

    double exit_lambda {0};

    // Exclude all existing nodes
    std::unordered_set<std::string> nodes_added {};
    for (int i {0}; i < genotype.length(); ++i) {
      if (genotype.at(i) == 1) {
        std::vector<int> node = node_labels.at(i);
        nodes_added.insert(node_to_string(node));
      }
    }

    // Add nodes
    for (int i {0}; i < genotype.length(); ++i) { //Loop through all nodes in the tree
      if (genotype.at(i) == 1) { //Check that the node is in the tree

        std::vector<int> node = node_labels.at(i); //Extract node label
        // Rcout << "Now at node : " << node_to_string(node) << "\n";
        std::vector<int> nodes_to_add = get_exit_set(n, node); //The set of nodes that could happen next are those that have not happened before
        // Add children
        for (int j : nodes_to_add) {
          node.push_back(j);
          std::string node_tag = node_to_string(node);
          auto temp = nodes_added.insert(node_tag);
          if (temp.second) { //if the node does not exist already, then add it
            // Rcout << "Added child : " << node_tag << "\n";
            exit_lambda += get_lambda(node, Theta);
          }
          node.pop_back();
        }

        // Add siblings
        node.pop_back();
        std::unordered_set<std::string>::const_iterator found = nodes_added.find(node_to_string(node));
        if (found == nodes_added.end()) {
          // if (!(nodes_added.contains(node_to_string(node)))) { //c++20
          for (int j : nodes_to_add) {
            node.push_back(j);
            std::string node_tag = node_to_string(node);
            auto temp = nodes_added.insert(node_tag);
            if (temp.second) {
              // Rcout << "Added sibling : " << node_tag << "\n";
              exit_lambda += get_lambda(node, Theta);
            }
            node.pop_back();
          }
        }
      }
    }
    return exit_lambda;
  }
}

// [[Rcpp::export]]
arma::sp_mat build_tr_mat(const int &n, const arma::mat &Theta, const IntegerMatrix &genotypes, 
                          const List &node_labels, const double &lambda_s) {

  int p = genotypes.ncol();
  int L = genotypes.nrow();
  arma::sp_mat tr_matrix(L,L);

  if (L == 1) {
    IntegerVector genotype = genotypes(0,_);
    tr_matrix(0,0) = get_exit_lambda(n, genotype, Theta, node_labels);

  } else {
    std::vector<long int> indices_in_lattice {};
    for (int i {0}; i < L; ++i) {
      long int row_sum = 0;
      for (int j = 0; j < p; ++j) {
        if (genotypes(i, j) == 1) {
          row_sum = row_sum + pow(2, p - j - 1);
        }
      }
      indices_in_lattice.push_back(row_sum);
    }

    for (int i {0}; i < L; ++i) {

      IntegerVector genotype = genotypes(i,_);
      tr_matrix(i,i) = get_exit_lambda(n, genotype, Theta, node_labels);

      for (int j = 0; j < p; ++j) {
        if (genotypes(i, j) == 0) {
          long int next_geno_idx = indices_in_lattice[i] + pow(2, p - j - 1);
          std::vector<long int>::iterator next_geno_idx_itr = std::find(indices_in_lattice.begin(), indices_in_lattice.end(), next_geno_idx);
          if (next_geno_idx_itr != indices_in_lattice.end()) {
            std::vector<int> node = node_labels.at(j);
            tr_matrix(next_geno_idx_itr-indices_in_lattice.begin(), i) = -get_lambda(node, Theta);
          }
        }
      }
    }

  }

  tr_matrix.diag() += lambda_s; // matrix V_T
  return tr_matrix;
}

// [[Rcpp::export]]
double compute_obs_ll(const arma::sp_mat &tr_mat, const double &lambda_s) {
  int p = tr_mat.n_rows;
  if (p == 1) {
    return log(lambda_s) - log(tr_mat(0,0));
  } else {
    arma::mat tr_mat_ = arma::mat(tr_mat);
    arma::vec b(p, fill::zeros);
    b(0) = 1;
    arma::vec temp = solve(trimatl(tr_mat_), b, arma::solve_opts::allow_ugly);
    return log(lambda_s) + log(temp(p-1) + 1e-10);
  }
}


int find_node_idx(const std::vector<int> &node, const List &node_labels) {

  for (int i {0}; i < node_labels.length(); ++i) {
    std::vector<int> temp = node_labels.at(i);
    if (node == temp) {
      return i;
    }
  }

  return 0;
}


arma::sp_mat dS_dlambda(std::vector<int> node, const bool &in_tree,
                        const IntegerMatrix &genotypes, const List &node_labels) {

  int p = genotypes.ncol();
  int L = genotypes.nrow();
  arma::sp_mat d_tr_mat(L,L);

  if (in_tree) {

    std::vector<long int> indices_in_lattice {};
    for (int i {0}; i < L; ++i) {
      long int row_sum = 0;
      for (int j = 0; j < p; ++j) {
        if (genotypes(i, j) == 1) {
          row_sum = row_sum + pow(2, p - j - 1);
        }
      }
      indices_in_lattice.push_back(row_sum);
    }

    int idx = find_node_idx(node, node_labels);
    if (node.size() == 1) {
      for (int i {0}; i < L; ++i) {
        IntegerVector genotype = genotypes(i,_);
        if (genotypes(i, idx) == 0) {
          d_tr_mat(i,i) = 1;

          long int next_geno_idx = indices_in_lattice[i] + pow(2, p - idx - 1);
          std::vector<long int>::iterator next_geno_idx_itr = std::find(indices_in_lattice.begin(), indices_in_lattice.end(), next_geno_idx);
          if (next_geno_idx_itr != indices_in_lattice.end()) {
            d_tr_mat(next_geno_idx_itr-indices_in_lattice.begin(), i) = -1;
          }
        }
      }
    } else {
      node.pop_back();
      int pa_idx = find_node_idx(node, node_labels);
      for (int i {0}; i < L; ++i) {
        IntegerVector genotype = genotypes(i,_);
        if ((genotypes(i, idx) == 0) && (genotypes(i, pa_idx) == 1)) {
          d_tr_mat(i,i) = 1;

          long int next_geno_idx = indices_in_lattice[i] + pow(2, p - idx - 1);
          std::vector<long int>::iterator next_geno_idx_itr = std::find(indices_in_lattice.begin(), indices_in_lattice.end(), next_geno_idx);
          if (next_geno_idx_itr != indices_in_lattice.end()) {
            d_tr_mat(next_geno_idx_itr-indices_in_lattice.begin(), i) = -1;
          }
        }
      }
    }

  } else {

    if (node.size() == 1) {
      d_tr_mat.diag() += 1;

    } else {
      node.pop_back();
      int pa_idx = find_node_idx(node, node_labels);
      for (int i {0}; i < L; ++i) {
        IntegerVector genotype = genotypes(i,_);
        if (genotypes(i, pa_idx) == 1) {
          d_tr_mat(i,i) = 1;
        }
      }
    }
  }

  return d_tr_mat;
}

arma::sp_mat build_AS(const arma::sp_mat &tr_mat, const arma::sp_mat &d_tr_mat) {

  int L = tr_mat.n_rows;
  arma::sp_mat zero_mat(L,L);
  arma::sp_mat top_mat = join_rows(tr_mat, zero_mat);
  arma::sp_mat bottom_mat = join_rows(d_tr_mat, tr_mat);
  arma::sp_mat AS = join_cols(top_mat, bottom_mat);
  return AS;

}

long double compute_dp_dlambda(const arma::sp_mat &tr_mat, const std::vector<int> &node, 
                               const bool &in_tree, const double &lambda_s,
                               const IntegerMatrix &genotypes, const List &node_labels) {

  arma::sp_mat d_tr_mat = dS_dlambda(node, in_tree, genotypes, node_labels);
  arma::sp_mat AS = build_AS(tr_mat, d_tr_mat);
  int L = AS.n_rows;
  arma::mat AS_ = arma::mat(AS);
  arma::vec b(L, fill::zeros);
  b(0) = 1;
  arma::vec temp = solve(trimatl(AS_), b, arma::solve_opts::allow_ugly);
  return temp(L-1) * lambda_s;
}

void tree_E_step_(const arma::mat &Theta, NumericVector &time_diffs, 
                  const double &prob_p, const int &n, const IntegerVector &nodes, 
                  const List &children, const LogicalVector &in_tree, const arma::sp_mat &tr_mat, 
                  const double &lambda_s, const IntegerMatrix &genotypes, const List &node_labels,
                  std::vector<int> pathway = {}, int current_pos = 0) {

  if (current_pos != 0) {
    int i = nodes.at(current_pos) - 1; // c++ indexing starts from 0

    std::vector<int> node = pathway;
    node.push_back(i + 1);
    long double dp_dl = compute_dp_dlambda(tr_mat, node, in_tree.at(current_pos), lambda_s, genotypes, node_labels);
    if (pathway.size() == 0) {

      time_diffs.at(current_pos) = exp(- Theta(i,i)) - dp_dl / prob_p;

    } else {

      double lambda_pi_i = get_lambda(node, Theta);
      time_diffs.at(current_pos) = 1 / lambda_pi_i - dp_dl / prob_p;

    }
  }

  IntegerVector ch {};
  if (children.at(current_pos) != R_NilValue) {

    ch = children.at(current_pos);

    if (current_pos != 0) {
      pathway.push_back(nodes.at(current_pos));
    }

    for (int pos : ch) {
      pos -= 1; // c++ indexing starts from 0
      tree_E_step_(Theta, time_diffs, prob_p, n,
                   nodes, children, in_tree, tr_mat,
                   lambda_s, genotypes, node_labels,
                   pathway, pos);
    }

  }
}

NumericVector tree_E_step(const arma::mat &Theta, const int &n, const double &lambda_s,
                          const IntegerVector &nodes, const List &children, const LogicalVector &in_tree,
                          double &obs_ll, const IntegerMatrix &genotypes, const List &node_labels) {

  arma::sp_mat tr_mat = build_tr_mat(n, Theta, genotypes, node_labels, lambda_s);
  NumericVector time_diffs(nodes.length());
  obs_ll = compute_obs_ll(tr_mat, lambda_s);
  tree_E_step_(Theta, time_diffs, exp(obs_ll), n,
               nodes, children, in_tree, tr_mat,
               lambda_s, genotypes, node_labels);
  return time_diffs;

}


double tree_imp_samp_once(const arma::mat &Theta, const IntegerVector &nodes,
                          const List &children, const LogicalVector &in_tree,
                          int n, NumericVector &time_diffs, const double &t_s,
                          IntegerVector pathway = {}, int current_pos = 0, double t_pa = 0) {

  double log_importance_weight = 0;

  if (current_pos != 0) {
    int i = nodes.at(current_pos) - 1; // c++ indexing starts from 0

    // compute parameter lambda
    double lambda_i;
    if (pathway.length() == 0) {
      lambda_i = exp(Theta(i,i));
    } else {
      double temp {0};
      for (int j : pathway) {
        temp += Theta(i, j - 1); // c++ indexing starts from 0
      }
      lambda_i = exp(Theta(i,i) + temp);
    }

    double texp_interval = t_s - t_pa;

    // sample time differences and compute importance weights
    if (in_tree.at(current_pos)) {
      time_diffs.at(current_pos) = rtexp(lambda_i, 0, texp_interval);
      log_importance_weight = log(1 - exp(- lambda_i * texp_interval));

    } else {
      time_diffs.at(current_pos) = rtexp(lambda_i, texp_interval, R_PosInf);
      log_importance_weight = - lambda_i * texp_interval;

    }
  }

  IntegerVector ch {};
  if (children.at(current_pos) != R_NilValue) {

    ch = children.at(current_pos);

    if (current_pos != 0) {
      pathway.push_back(nodes.at(current_pos));
    }

    for (int pos : ch) {
      pos -= 1; // c++ indexing starts from 0
      log_importance_weight += tree_imp_samp_once(Theta, nodes, children, in_tree,
                                                  n, time_diffs, t_s, pathway, pos,
                                                  t_pa + time_diffs.at(current_pos));
    }

    return log_importance_weight;

  } else {

    return log_importance_weight;

  }

}

NumericVector tree_importance_sampling(const arma::mat &Theta, const IntegerVector &nodes,
                                       const List &children, const LogicalVector &in_tree,
                                       const int &n, const int &M, const double &lambda_s,
                                       double &obs_ll) {

  NumericVector expected_time_diffs(nodes.length());
  double importance_weight_sum {0};

  for (int m {0}; m < M; ++m) {

    double t_s = rexp(1, lambda_s)[0];
    NumericVector time_diffs(nodes.length());
    double log_importance_weight = tree_imp_samp_once(Theta, nodes, children, in_tree, n, time_diffs, t_s);
    double importance_weight = exp(log_importance_weight);
    time_diffs = time_diffs * importance_weight;
    expected_time_diffs = expected_time_diffs + time_diffs;
    importance_weight_sum += importance_weight;

  }

  expected_time_diffs = expected_time_diffs / (importance_weight_sum);
  obs_ll = log(importance_weight_sum) - log(M);

  return expected_time_diffs;
}


// [[Rcpp::export]]
void update_timed_trees(const int &n, const int &N, List &timed_trees, 
                        const arma::mat &Theta, const double &lambda_s, const int &M,
                        NumericVector &log_prob_vec, const List &comp_geno_vec, 
                        const List &node_labels_vec, const int &nr_exact) {
  
  // Loop through all trees using exact inference
  for (int i {0}; i < nr_exact; ++i) {
    
    List tree = timed_trees.at(i);
    tree["time_diffs"] = tree_E_step(Theta, n, lambda_s, tree["nodes"], tree["children"], 
                           tree["in_tree"], log_prob_vec.at(i), comp_geno_vec.at(i), node_labels_vec.at(i));
    timed_trees.at(i) = tree;
    
  }
  
  // Loop through all trees using approximate inference
  for (int i {nr_exact}; i < N; ++i) {
    
    List tree = timed_trees.at(i);
    tree["time_diffs"] = tree_importance_sampling(Theta, tree["nodes"], tree["children"], 
                                        tree["in_tree"], n, M, lambda_s, log_prob_vec.at(i));
    timed_trees.at(i) = tree;
    
  }
  
}


