#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// H - matrix of size k x n
// W - matrix of size m x k
// V - matrix of size m x k
// E - sparse matrix of size m x n
// [[Rcpp::export]]
double objErr_i(
        const arma::mat& H,
        const arma::mat& W,
        const arma::mat& V,
        const arma::sp_mat& E,
        const double& lambda) {
    arma::mat L = W + V; // m x k
    double sqnormE = arma::norm<arma::sp_mat>(E, "fro");
    sqnormE *= sqnormE;
    arma::mat LtL = L.t() * L; // k x k
    arma::mat HtH = H * H.t(); // k x k
    double TrLtLHtH = arma::trace(LtL * HtH);
    arma::mat EtL = E.t() * L; // n x k
    double TrHtEtL = arma::trace(H * EtL);
    arma::mat VtV = V.t() * V; // k x k
    double TrVtVHtH = arma::trace(VtV * HtH);
    double err = sqnormE - 2 * TrHtEtL + TrLtLHtH + lambda * TrVtVHtH;
    return err;
}
