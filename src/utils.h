#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double rtexp(const double &lambda, const double &lbound, const double &ubound);
std::vector<int> my_setdiff(std::vector<int> v1, std::vector<int> v2);
double get_lambda(std::vector<int> node, const arma::mat &Theta);
std::string node_to_string(const std::vector<int> &node);
std::vector<int> get_exit_set(const int &n, const std::vector<int> &node);

#endif //UTILS_H
