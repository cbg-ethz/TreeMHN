#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double rtexp(double lambda, double lbound, double ubound);
std::vector<int> my_setdiff(std::vector<int> v1, std::vector<int> v2);
double get_lambda(std::vector<int> node, arma::mat Theta);
std::string node_to_string(std::vector<int> node);
std::vector<int> get_exit_set(int n, std::vector<int> node);

#endif //UTILS_H
