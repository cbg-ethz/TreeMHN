#include <RcppArmadillo.h>
#include "utils.h"
#include "em_fns.h"
#include <iostream>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double l1_penalty(arma::mat Theta, const double &gamma) {
  Theta.diag(0).fill(0);
  return gamma * accu(sqrt(Theta % Theta));
}


arma::mat with_l1_penalty_grad(arma::mat Theta, const double &gamma, 
                               const arma::mat &Theta_grad, const int &n) {
  Theta.diag(0).fill(0);
  return Theta_grad - gamma * arma::sign(Theta);
}


double full_tree_log_score(const arma::mat &Theta, const IntegerVector &nodes, const List &children,
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

arma::mat full_tree_grad(const arma::mat &Theta, const int &n, const IntegerVector &nodes, const List &children,
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

arma::mat obs_tree_grad(const arma::mat &Theta, const int &n, const IntegerVector &nodes, const List &children,
                        const LogicalVector &in_tree, const arma::sp_mat &tr_mat, const double &lambda_s,
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


double prob_empty_tree(const arma::vec &lambdas, const double &lambda_s) {
  
  return lambda_s / (lambda_s + sum(lambdas));
  
}


double prob_one_tree(const double &n, const arma::mat &exp_Theta, const double &lambda_s) {
  
  double p = 0;
  double sum_lambdas = sum(exp_Theta.diag(0));
  for (int i {0}; i < n; ++i) {
    double temp_denom = 0;
    for (int j {0}; j < n; ++j) {
      if (j != i) {
        temp_denom += exp_Theta(j, j) * (1 + exp_Theta(j, i));
      }
    }
    p += exp_Theta(i, i) / (lambda_s + sum_lambdas) * lambda_s / (lambda_s + temp_denom);
  }
  return p;
  
}


// [[Rcpp::export]]
double full_MHN_objective(const arma::vec &Theta, const List &trees, const double &gamma,
                          const int &n, const int &N, const double &lambda_s, 
                          const IntegerVector &to_mask, const NumericVector &weights,
                          const int &N_patients, double smallest_tree_size = 1) {
  
  arma::mat Theta_ = arma::reshape(Theta, n, n);
  
  // Mask elements if needed
  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_(i - 1) = 0;
    }
  }
  
  double log_score = - l1_penalty(Theta_, gamma);;
  
  // Loop through all trees
  for (int i {0}; i < trees.length(); ++i) {
    List tree = trees.at(i);
    log_score += weights.at(i) * full_tree_log_score(Theta_, tree["nodes"], tree["children"], tree["time_diffs"]);
  }
  
  arma::mat exp_Theta = exp(Theta_);
  arma::vec lambdas = exp_Theta.diag(0);
  
  if (smallest_tree_size == 1) {
    log_score -= N_patients * log(1 - prob_empty_tree(lambdas, lambda_s));
  }
  
  if (smallest_tree_size == 2) {
    log_score -= N_patients * log(1 - prob_empty_tree(lambdas, lambda_s)) - prob_one_tree(n, exp_Theta, lambda_s);
  }
  
  return log_score;
  
}

// [[Rcpp::export]]
arma::mat full_MHN_grad(const arma::vec &Theta, const List &trees, const double &gamma,
                        const int &n, const int &N, const double &lambda_s, 
                        const IntegerVector &to_mask, const NumericVector &weights,
                        const int &N_patients, double smallest_tree_size = 1) {
  
  arma::mat Theta_ = arma::reshape(Theta, n, n);
  
  // Mask elements if needed
  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_(i - 1) = 0;
    }
  }
  
  arma::mat Theta_grad(n,n,fill::zeros);
  // Loop through all trees
  for (int i {0}; i < trees.length(); ++i) {
    List tree = trees.at(i);
    Theta_grad += weights.at(i) * full_tree_grad(Theta_, n, tree["nodes"], tree["children"], tree["time_diffs"]);
  }
  Theta_grad = with_l1_penalty_grad(Theta_, gamma, Theta_grad, n); 
  
  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_grad(i - 1) = 0;
    }
  }
  
  arma::mat exp_Theta = exp(Theta_);
  arma::vec lambdas = exp_Theta.diag(0);
  double p_empty = prob_empty_tree(lambdas, lambda_s);
  
  if (smallest_tree_size >= 1) {
    Theta_grad.diag(0) -= N_patients * pow(p_empty, 2) / lambda_s / (1 - p_empty + 1e-10) * lambdas;
  }
  
  // if (smallest_tree_size == 2) {
  //   double p_one = prob_one_tree(n, exp_Theta, lambda_s);
  //   Theta_grad.diag(0) -= N_patients * (diag(p_empty^2 / lambda_s * lambdas) - dp_one(n, exp_Theta, lambda_s)) / (1 - p_empty - p_one + 1e-10)
  // }
  
  return Theta_grad;
  
}


// [[Rcpp::export]]
double obs_MHN_objective(const arma::vec &Theta, const int &n, const int &N, 
                         const double &lambda_s, const List &trees, const double &gamma, 
                         List &tr_mat_vec, NumericVector &log_prob_vec,
                         const List &comp_geno_vec, const List &node_labels_vec, 
                         const IntegerVector &to_mask, const NumericVector &weights,
                         const int &N_patients, double smallest_tree_size = 1) {
  
  arma::mat Theta_ = arma::reshape(Theta, n, n);
  
  // Mask elements if needed
  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_(i - 1) = 0;
    }
  }
  
  double log_score = - l1_penalty(Theta_, gamma);
  for (int i {0}; i < N; ++i) {

    arma::sp_mat tr_mat = build_tr_mat(n, Theta_, comp_geno_vec.at(i), node_labels_vec.at(i), lambda_s);
    tr_mat_vec.at(i) = tr_mat;
    double p = compute_obs_ll(tr_mat, lambda_s);
    log_score += weights.at(i) * p;
    log_prob_vec.at(i) = p;
  }
  
  arma::mat exp_Theta = exp(Theta_);
  arma::vec lambdas = exp_Theta.diag(0);
  
  if (smallest_tree_size == 1) {
    log_score -= sum(weights) * log(1 - prob_empty_tree(lambdas, lambda_s));
  }
  
  if (smallest_tree_size == 2) {
    log_score -= sum(weights) * log(1 - prob_empty_tree(lambdas, lambda_s) - prob_one_tree(n, exp_Theta, lambda_s));
  }
  
  return log_score;
  
}


// [[Rcpp::export]]
arma::mat obs_MHN_grad(const arma::vec &Theta, const int &n, const int &N, 
                       const double &lambda_s, const List &trees, const double &gamma, 
                       const List &tr_mat_vec, const NumericVector &log_prob_vec,
                       const List &comp_geno_vec, const List &node_labels_vec, 
                       const IntegerVector &to_mask, const NumericVector &weights,
                       const int &N_patients, double smallest_tree_size = 1) {
  
  arma::mat Theta_ = arma::reshape(Theta, n, n);
  
  // Mask elements if needed
  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_(i - 1) = 0;
    }
  }
  
  arma::mat Theta_grad(n,n,fill::zeros);
  // Loop through all trees
  for (int i {0}; i < N; ++i) {
    
    List tree = trees.at(i);
    arma::mat temp = obs_tree_grad(Theta_, n, tree["nodes"], tree["children"], 
                                   tree["in_tree"], tr_mat_vec.at(i),
                                   lambda_s, comp_geno_vec.at(i), node_labels_vec.at(i));
    
    double p = log_prob_vec.at(i);
    Theta_grad += weights.at(i) * arma::sign(temp) % exp(log(arma::abs(temp)) - p);
  }
  Theta_grad = with_l1_penalty_grad(Theta_, gamma, Theta_grad, n);
  
  if (to_mask.length() != 0) {
    for (int i : to_mask) {
      Theta_grad(i - 1) = 0;
    }
  }
  
  arma::mat exp_Theta = exp(Theta_);
  arma::vec lambdas = exp_Theta.diag(0);
  double p_empty = prob_empty_tree(lambdas, lambda_s);
  
  if (smallest_tree_size >= 1) {
    Theta_grad.diag(0) -= sum(weights) * pow(p_empty, 2) / lambda_s / (1 - p_empty + 1e-10) * lambdas;
  }
  
  // if (smallest_tree_size == 2) {
  //   double p_one = prob_one_tree(n, exp_Theta, lambda_s);
  //   Theta_grad.diag(0) -= N_patients * (diag(p_empty^2 / lambda_s * lambdas) - dp_one(n, exp_Theta, lambda_s)) / (1 - p_empty - p_one + 1e-10)
  // }
  
  return Theta_grad;
}

