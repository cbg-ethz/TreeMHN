#include <RcppArmadillo.h>
#include "utils.h"
#include "em_fns.h"
#include <iostream>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double l1_penalty(arma::mat Theta, double gamma) {
  Theta.diag(0).fill(0);
  return gamma * accu(sqrt(Theta % Theta));
}


NumericMatrix with_l1_penalty_grad(arma::mat Theta, double gamma, arma::mat Theta_grad, int n) {
  Theta.diag(0).fill(0);
  Theta_grad = Theta_grad - gamma * arma::sign(Theta);
  return wrap(Theta_grad);
}


double full_tree_log_score(arma::mat Theta, const IntegerVector &nodes, const List &children,
                           const NumericVector &time_diffs, IntegerVector pathway = {}, int current_pos = 0) {

  double log_score {0};

  if (current_pos != 0) {

    int i = nodes.at(current_pos) - 1; // c++ indexing starts from 0

    if (pathway.length() == 0) {
      log_score += Theta(i,i) - exp(Theta(i,i)) * time_diffs.at(current_pos);
    } else {

      double temp {0};
      for (int j : pathway) {
        temp += Theta(i, j - 1); // c++ indexing starts from 0
      }
      log_score += Theta(i,i) + temp - exp(Theta(i,i) + temp) * time_diffs.at(current_pos);
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
      log_score += full_tree_log_score(Theta, nodes, children, time_diffs, pathway, pos);
    }

    return log_score;

  } else {

    return log_score;

  }

}

arma::mat full_tree_grad(arma::mat Theta, int n, const IntegerVector &nodes, const List &children,
                         const NumericVector &time_diffs, IntegerVector pathway = {}, int current_pos = 0) {

  arma::mat Theta_grad(n,n,fill::zeros);

  if (current_pos != 0) {

    int i = nodes.at(current_pos) - 1; // c++ indexing starts from 0

    if (pathway.length() == 0) {

      Theta_grad(i,i) += 1 - exp(Theta(i,i)) * time_diffs.at(current_pos);

    } else {

      double temp {0};

      for (int j : pathway) {
        temp += Theta(i, j - 1); // c++ indexing starts from 0
      }

      Theta_grad(i,i) += 1 - exp(Theta(i,i) + temp) * time_diffs.at(current_pos);

      for (int j : pathway) {
        Theta_grad(i, j - 1) += 1 - exp(Theta(i,i) + temp) * time_diffs.at(current_pos);
      }

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
      Theta_grad += full_tree_grad(Theta, n, nodes, children, time_diffs, pathway, pos);
    }

    return Theta_grad;

  } else {

    return Theta_grad;

  }
}

// [[Rcpp::export]]
double full_MHN_objective_(NumericVector Theta, const List &trees, double gamma,
                           int n, int N, double lambda_s, IntegerVector to_mask, NumericVector weights) {

  double log_score {0};

  // Convert Theta vector to matrix
  Theta.attr("dim") = Dimension(n, n);
  arma::mat Theta_ = as<arma::mat>(Theta);

  // Mask elements if needed
  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_(i - 1) = 0;
    }
  }

  // Loop through all trees
  for (int i {0}; i < trees.length(); ++i) {
    List tree = trees.at(i);
    IntegerVector nodes = tree["nodes"];
    List children = tree["children"];
    NumericVector time_diffs =  tree["time_diffs"];
    double weight = weights.at(i); //New
    log_score += weight * full_tree_log_score(Theta_, nodes, children, time_diffs); //New
  }

  if (log_score == -R_PosInf) {
    log_score = -1e10;
  }

  if (log_score == R_PosInf) {
    log_score = -1e2;
  }

  log_score = log_score - l1_penalty(Theta_, gamma);
  return log_score;
}


// [[Rcpp::export]]
NumericMatrix full_MHN_grad_(NumericVector Theta, const List &trees, double gamma,
                             int n, int N, double lambda_s, IntegerVector to_mask, NumericVector weights) {

  // Convert Theta vector to matrix
  Theta.attr("dim") = Dimension(n, n);

  arma::mat Theta_grad(n,n,fill::zeros);
  arma::mat Theta_ = as<arma::mat>(Theta);

  // Mask elements if needed
  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_(i - 1) = 0;
    }
  }

  // Loop through all trees
  for (int i {0}; i < trees.length(); ++i) {
    List tree = trees.at(i);
    IntegerVector nodes = tree["nodes"];
    List children = tree["children"];
    NumericVector time_diffs =  tree["time_diffs"];
    double weight = weights.at(i); //New
    Theta_grad += weight * full_tree_grad(Theta_, n, nodes, children, time_diffs); //New
  }

  NumericMatrix Theta_grad_ = with_l1_penalty_grad(Theta_, gamma, Theta_grad, n);

  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_grad(i - 1) = 0;
    }
  }

  return Theta_grad_;
}

// [[Rcpp::export]]
double obs_MHN_objective_(NumericVector Theta, int n, int N, double lambda_s,
                          const List &trees, double gamma, List &obj_grad_help, 
                          IntegerVector to_mask, NumericVector weights) {

  // Convert Theta vector to matrix
  Theta.attr("dim") = Dimension(n, n);
  arma::mat Theta_ = as<arma::mat>(Theta);

  // Mask elements if needed
  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_(i - 1) = 0;
    }
  }

  List tr_mat_vec = obj_grad_help["tr_mat_vec"];
  List comp_geno_vec = obj_grad_help["comp_geno_vec"];
  List node_labels_vec = obj_grad_help["node_labels_vec"];
  NumericVector log_prob_vec = obj_grad_help["log_prob_vec"];

  double log_score {0};
  for (int i {0}; i < N; ++i) {
    IntegerMatrix genotypes = comp_geno_vec.at(i);
    List node_labels = node_labels_vec.at(i);
    arma::sp_mat tr_mat = build_tr_mat(n, Theta_, genotypes, node_labels);
    tr_mat_vec.at(i) = tr_mat;
    double p = compute_obs_ll(tr_mat, lambda_s);
    double weight = weights.at(i); //New
    log_score += weight * p; //New
    log_prob_vec.at(i) = p;
  }
  obj_grad_help["tr_mat_vec"] = tr_mat_vec;
  obj_grad_help["log_prob_vec"] = log_prob_vec;

  log_score = log_score - l1_penalty(Theta_, gamma);
  return log_score;

}

arma::mat obs_tree_grad(const arma::mat Theta, int n, const IntegerVector &nodes, const List &children,
                        const LogicalVector &in_tree, arma::sp_mat tr_mat, double lambda_s,
                        const IntegerMatrix &genotypes, const List &node_labels,
                        std::vector<int> pathway = {}, int current_pos = 0) {

  arma::mat Theta_grad(n,n,fill::zeros);

  if (current_pos != 0) {
    int i = nodes.at(current_pos) - 1; // c++ indexing starts from 0

    std::vector<int> node = pathway;
    node.push_back(i + 1);
    long double dp_dl = compute_dp_dlambda(tr_mat, node, in_tree.at(current_pos), lambda_s, genotypes, node_labels);
    if (pathway.size() == 0) {
      Theta_grad(i,i) += dp_dl * exp(Theta(i, i));
    } else {
      double lambda_pi_i = get_lambda(node, Theta);
      Theta_grad(i,i) += dp_dl * lambda_pi_i;
      for (int j : pathway) {
        Theta_grad(i, j - 1) += dp_dl * lambda_pi_i;
      }
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
      Theta_grad += obs_tree_grad(Theta, n, nodes, children, in_tree,
                                  tr_mat, lambda_s, genotypes, node_labels,
                                  pathway, pos);
    }

    return Theta_grad;

  } else {

    return Theta_grad;

  }
}

// [[Rcpp::export]]
NumericMatrix obs_MHN_grad_(NumericVector Theta, int n, int N, double lambda_s,
                            const List &trees, double gamma, const List &obj_grad_help, 
                            IntegerVector to_mask, NumericVector weights) {

  // Convert Theta vector to matrix
  Theta.attr("dim") = Dimension(n, n);
  arma::mat Theta_ = as<arma::mat>(Theta);

  // Mask elements if needed
  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_(i - 1) = 0;
    }
  }

  arma::mat Theta_grad(n,n,fill::zeros);
  List tr_mat_vec = obj_grad_help["tr_mat_vec"];
  List comp_geno_vec = obj_grad_help["comp_geno_vec"];
  List node_labels_vec = obj_grad_help["node_labels_vec"];
  NumericVector log_prob_vec = obj_grad_help["log_prob_vec"];

  // Loop through all trees
  for (int i {0}; i < N; ++i) {
    List tree = trees.at(i);
    IntegerVector nodes = tree["nodes"];
    List children = tree["children"];
    LogicalVector in_tree =  tree["in_tree"];
    arma::sp_mat tr_mat = tr_mat_vec.at(i);
    IntegerMatrix genotypes = comp_geno_vec.at(i);
    List node_labels = node_labels_vec.at(i);
    double weight = weights.at(i); //New

    arma::mat temp = obs_tree_grad(Theta_, n, nodes, children, in_tree, tr_mat,
                                   lambda_s, genotypes, node_labels);

    double p = log_prob_vec.at(i);
    Theta_grad += weight * arma::sign(temp) % exp(log(arma::abs(temp)) - p); //New
  }

  NumericMatrix Theta_grad_ = with_l1_penalty_grad(Theta_, gamma, Theta_grad, n);

  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_grad(i - 1) = 0;
    }
  }

  return Theta_grad_;
}
