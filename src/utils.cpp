#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double rtexp(double lambda, double lbound, double ubound) {

  if (lbound > ubound) {
    // throw std::runtime_error("Something went wrong in rtexp!");
    return 1e-10;
  }

  double u = runif(1)[0];
  double texp  = -1 / lambda * (log(exp(-lambda * lbound) -
                                (exp(-lambda * lbound) - exp(-lambda * ubound)) * u));

  double eps = abs((ubound - lbound) / 2);
  eps = (eps == R_PosInf) ? lbound * 2 : eps;
  texp = (texp == 0) ? eps : texp;
  texp = (texp == R_PosInf) ? eps : texp;

  return texp;
}


std::vector<int> my_setdiff(std::vector<int> v1, std::vector<int> v2) {
  std::vector<int> new_vec {};
  std::sort(v2.begin(),v2.end());
  std::set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(),
                      std::inserter(new_vec, new_vec.begin()));
  return new_vec;
}

// [[Rcpp::export]]
double get_lambda(std::vector<int> node, arma::mat Theta) {
  double lambda {0};
  if (node.size() == 1) {
    int i = node.at(0) - 1;
    lambda = exp(Theta(i,i));
  } else {
    int i = node.back() - 1;
    node.pop_back();
    double temp {0};
    for (int j : node) {
      temp += Theta(i,j - 1);
    }
    lambda = exp(Theta(i,i) + temp);
  }
  return lambda;
}

std::string node_to_string(std::vector<int> node) {
  std::string s {""};
  for (int i: node) {
    s += std::to_string(i);
    s += "_";
  }
  return s;
}


std::vector<int> get_exit_set(int n, std::vector<int> node) {

  std::vector<int> temp(n);
  std::iota(temp.begin(), temp.end(), 1);
  return my_setdiff(temp, node);

}
