#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export()]]
arma::mat colNormalize_dense_cpp(arma::mat& x, const arma::uword L) {
    arma::mat result(x);
    for (int j = 0; j < x.n_cols; ++j) {
        double norm = arma::norm(x.col(j), L);
        for (int i = 0; i < x.n_rows; ++i) {
            result(i, j) /= norm;
        }
    }
    return result;
}

// x: n features by n selected factors
// group: n-selected-factor integer vector, pre-transformed to 0-base,
// from upstream kmeans clustering.
// [[Rcpp::export()]]
arma::mat colAggregateMedian_dense_cpp(const arma::mat& x, const arma::uvec& group, const arma::uword n) {
    arma::mat result(x.n_rows, n);
    for (int i = 0; i < n; ++i) {
        arma::uvec idx = arma::find(group == i);
        arma::mat sub_x = x.cols(idx);
        arma::vec median = arma::median(sub_x, 1);
        result.col(i) = median;
    }
    return result;
}


