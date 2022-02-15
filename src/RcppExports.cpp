// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// build_tr_mat
arma::sp_mat build_tr_mat(const int& n, const arma::mat& Theta, const IntegerMatrix& genotypes, const List& node_labels, const double& lambda_s);
RcppExport SEXP _TreeMHN_build_tr_mat(SEXP nSEXP, SEXP ThetaSEXP, SEXP genotypesSEXP, SEXP node_labelsSEXP, SEXP lambda_sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type genotypes(genotypesSEXP);
    Rcpp::traits::input_parameter< const List& >::type node_labels(node_labelsSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda_s(lambda_sSEXP);
    rcpp_result_gen = Rcpp::wrap(build_tr_mat(n, Theta, genotypes, node_labels, lambda_s));
    return rcpp_result_gen;
END_RCPP
}
// compute_obs_ll
double compute_obs_ll(const arma::sp_mat& tr_mat, const double& lambda_s);
RcppExport SEXP _TreeMHN_compute_obs_ll(SEXP tr_matSEXP, SEXP lambda_sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type tr_mat(tr_matSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda_s(lambda_sSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_obs_ll(tr_mat, lambda_s));
    return rcpp_result_gen;
END_RCPP
}
// update_timed_trees
void update_timed_trees(const int& n, const int& N, List& timed_trees, const arma::mat& Theta, const double& lambda_s, const int& M, const List& comp_geno_vec, const List& node_labels_vec, const int& nr_exact);
RcppExport SEXP _TreeMHN_update_timed_trees(SEXP nSEXP, SEXP NSEXP, SEXP timed_treesSEXP, SEXP ThetaSEXP, SEXP lambda_sSEXP, SEXP MSEXP, SEXP comp_geno_vecSEXP, SEXP node_labels_vecSEXP, SEXP nr_exactSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< List& >::type timed_trees(timed_treesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda_s(lambda_sSEXP);
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const List& >::type comp_geno_vec(comp_geno_vecSEXP);
    Rcpp::traits::input_parameter< const List& >::type node_labels_vec(node_labels_vecSEXP);
    Rcpp::traits::input_parameter< const int& >::type nr_exact(nr_exactSEXP);
    update_timed_trees(n, N, timed_trees, Theta, lambda_s, M, comp_geno_vec, node_labels_vec, nr_exact);
    return R_NilValue;
END_RCPP
}
// full_MHN_objective
double full_MHN_objective(const arma::vec& Theta, const List& trees, const double& gamma, const int& n, const int& N, const double& lambda_s, const IntegerVector& to_mask, const NumericVector& weights, const int& N_patients, double smallest_tree_size);
RcppExport SEXP _TreeMHN_full_MHN_objective(SEXP ThetaSEXP, SEXP treesSEXP, SEXP gammaSEXP, SEXP nSEXP, SEXP NSEXP, SEXP lambda_sSEXP, SEXP to_maskSEXP, SEXP weightsSEXP, SEXP N_patientsSEXP, SEXP smallest_tree_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type trees(treesSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda_s(lambda_sSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type to_mask(to_maskSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const int& >::type N_patients(N_patientsSEXP);
    Rcpp::traits::input_parameter< double >::type smallest_tree_size(smallest_tree_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(full_MHN_objective(Theta, trees, gamma, n, N, lambda_s, to_mask, weights, N_patients, smallest_tree_size));
    return rcpp_result_gen;
END_RCPP
}
// full_MHN_grad
arma::mat full_MHN_grad(const arma::vec& Theta, const List& trees, const double& gamma, const int& n, const int& N, const double& lambda_s, const IntegerVector& to_mask, const NumericVector& weights, const int& N_patients, double smallest_tree_size);
RcppExport SEXP _TreeMHN_full_MHN_grad(SEXP ThetaSEXP, SEXP treesSEXP, SEXP gammaSEXP, SEXP nSEXP, SEXP NSEXP, SEXP lambda_sSEXP, SEXP to_maskSEXP, SEXP weightsSEXP, SEXP N_patientsSEXP, SEXP smallest_tree_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const List& >::type trees(treesSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda_s(lambda_sSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type to_mask(to_maskSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const int& >::type N_patients(N_patientsSEXP);
    Rcpp::traits::input_parameter< double >::type smallest_tree_size(smallest_tree_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(full_MHN_grad(Theta, trees, gamma, n, N, lambda_s, to_mask, weights, N_patients, smallest_tree_size));
    return rcpp_result_gen;
END_RCPP
}
// obs_MHN_objective
double obs_MHN_objective(const arma::vec& Theta, const int& n, const int& N, const double& lambda_s, const List& trees, const double& gamma, List& tr_mat_vec, NumericVector& log_prob_vec, const List& comp_geno_vec, const List& node_labels_vec, const IntegerVector& to_mask, const NumericVector& weights, const int& N_patients, double smallest_tree_size);
RcppExport SEXP _TreeMHN_obs_MHN_objective(SEXP ThetaSEXP, SEXP nSEXP, SEXP NSEXP, SEXP lambda_sSEXP, SEXP treesSEXP, SEXP gammaSEXP, SEXP tr_mat_vecSEXP, SEXP log_prob_vecSEXP, SEXP comp_geno_vecSEXP, SEXP node_labels_vecSEXP, SEXP to_maskSEXP, SEXP weightsSEXP, SEXP N_patientsSEXP, SEXP smallest_tree_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda_s(lambda_sSEXP);
    Rcpp::traits::input_parameter< const List& >::type trees(treesSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< List& >::type tr_mat_vec(tr_mat_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type log_prob_vec(log_prob_vecSEXP);
    Rcpp::traits::input_parameter< const List& >::type comp_geno_vec(comp_geno_vecSEXP);
    Rcpp::traits::input_parameter< const List& >::type node_labels_vec(node_labels_vecSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type to_mask(to_maskSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const int& >::type N_patients(N_patientsSEXP);
    Rcpp::traits::input_parameter< double >::type smallest_tree_size(smallest_tree_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(obs_MHN_objective(Theta, n, N, lambda_s, trees, gamma, tr_mat_vec, log_prob_vec, comp_geno_vec, node_labels_vec, to_mask, weights, N_patients, smallest_tree_size));
    return rcpp_result_gen;
END_RCPP
}
// obs_MHN_grad
arma::mat obs_MHN_grad(const arma::vec& Theta, const int& n, const int& N, const double& lambda_s, const List& trees, const double& gamma, const List& tr_mat_vec, const NumericVector& log_prob_vec, const List& comp_geno_vec, const List& node_labels_vec, const IntegerVector& to_mask, const NumericVector& weights, const int& N_patients, double smallest_tree_size);
RcppExport SEXP _TreeMHN_obs_MHN_grad(SEXP ThetaSEXP, SEXP nSEXP, SEXP NSEXP, SEXP lambda_sSEXP, SEXP treesSEXP, SEXP gammaSEXP, SEXP tr_mat_vecSEXP, SEXP log_prob_vecSEXP, SEXP comp_geno_vecSEXP, SEXP node_labels_vecSEXP, SEXP to_maskSEXP, SEXP weightsSEXP, SEXP N_patientsSEXP, SEXP smallest_tree_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda_s(lambda_sSEXP);
    Rcpp::traits::input_parameter< const List& >::type trees(treesSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const List& >::type tr_mat_vec(tr_mat_vecSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type log_prob_vec(log_prob_vecSEXP);
    Rcpp::traits::input_parameter< const List& >::type comp_geno_vec(comp_geno_vecSEXP);
    Rcpp::traits::input_parameter< const List& >::type node_labels_vec(node_labels_vecSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type to_mask(to_maskSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const int& >::type N_patients(N_patientsSEXP);
    Rcpp::traits::input_parameter< double >::type smallest_tree_size(smallest_tree_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(obs_MHN_grad(Theta, n, N, lambda_s, trees, gamma, tr_mat_vec, log_prob_vec, comp_geno_vec, node_labels_vec, to_mask, weights, N_patients, smallest_tree_size));
    return rcpp_result_gen;
END_RCPP
}
// get_augmented_trees
List get_augmented_trees(int n, const List& trees);
RcppExport SEXP _TreeMHN_get_augmented_trees(SEXP nSEXP, SEXP treesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const List& >::type trees(treesSEXP);
    rcpp_result_gen = Rcpp::wrap(get_augmented_trees(n, trees));
    return rcpp_result_gen;
END_RCPP
}
// parse_trees
List parse_trees(std::string path);
RcppExport SEXP _TreeMHN_parse_trees(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(parse_trees(path));
    return rcpp_result_gen;
END_RCPP
}
// get_lambda
double get_lambda(std::vector<int> node, const arma::mat& Theta);
RcppExport SEXP _TreeMHN_get_lambda(SEXP nodeSEXP, SEXP ThetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type node(nodeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Theta(ThetaSEXP);
    rcpp_result_gen = Rcpp::wrap(get_lambda(node, Theta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TreeMHN_build_tr_mat", (DL_FUNC) &_TreeMHN_build_tr_mat, 5},
    {"_TreeMHN_compute_obs_ll", (DL_FUNC) &_TreeMHN_compute_obs_ll, 2},
    {"_TreeMHN_update_timed_trees", (DL_FUNC) &_TreeMHN_update_timed_trees, 9},
    {"_TreeMHN_full_MHN_objective", (DL_FUNC) &_TreeMHN_full_MHN_objective, 10},
    {"_TreeMHN_full_MHN_grad", (DL_FUNC) &_TreeMHN_full_MHN_grad, 10},
    {"_TreeMHN_obs_MHN_objective", (DL_FUNC) &_TreeMHN_obs_MHN_objective, 14},
    {"_TreeMHN_obs_MHN_grad", (DL_FUNC) &_TreeMHN_obs_MHN_grad, 14},
    {"_TreeMHN_get_augmented_trees", (DL_FUNC) &_TreeMHN_get_augmented_trees, 2},
    {"_TreeMHN_parse_trees", (DL_FUNC) &_TreeMHN_parse_trees, 1},
    {"_TreeMHN_get_lambda", (DL_FUNC) &_TreeMHN_get_lambda, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_TreeMHN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
