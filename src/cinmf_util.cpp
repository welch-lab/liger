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

// all_data: n features by nFactor * nRep factors
// knn: nFactor * nRep factors by k + 1 nearest neighbors. The first column is
// the index of the factor itself
// threshold: threshold of maximum mean euclidean distance from a factor to its
// nearest neighbors
// Return: The one-based index of selected factors
// [[Rcpp::export()]]
Rcpp::LogicalVector select_factor_cpp(
        arma::mat all_data,
        arma::mat knn,
        double threshold) {
    int k = knn.n_cols - 1;
    arma::vec obs(all_data.n_rows);
    arma::vec neighbor(all_data.n_rows);
    double dist, mean_dist;
    Rcpp::LogicalVector selected_factors(all_data.n_cols);
    for (int j = 0; j < knn.n_rows; ++j) {
        obs = all_data.col(j);
        mean_dist = 0;
        for (int i = 1; i < knn.n_cols; ++i) {
            neighbor = all_data.col(knn(j, i));
            dist = arma::norm(obs - neighbor, 2);
            mean_dist += dist;
        }
        mean_dist /= k;
        if (mean_dist < threshold) {
            selected_factors(j) = true;
        } else {
            selected_factors(j) = false;
        }
    }
    return selected_factors;
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


