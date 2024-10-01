#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// Credit to https://github.com/immunogenomics/harmony/blob/master/src/harmony.cpp
// Z_orig - D (PC) x N (cell) matrix, input embedding
// R      - K (cluster) x N matrix, cluster assignment probability matrix
// lambda - B element vector, ridge regression penalty
// Phi    - B x N sparse matrix, generated with Matrix::fac2sparse(lig$dataset)
// B      - number of batches
// N      - number of cells
// [[Rcpp::export()]]
arma::mat moe_correct_ridge_cpp(
        const arma::mat& Z_orig,
        const arma::mat& R,
        const arma::colvec& lambda,
        const arma::sp_mat& Phi,
        const unsigned int B,
        const unsigned int N
) {
    arma::sp_mat intcpt = zeros<arma::sp_mat>(1, N);
    intcpt = intcpt + 1;

    arma::sp_mat Phi_moe = join_cols(intcpt, Phi);
    arma::sp_mat Phi_moe_t = Phi_moe.t();

    arma::sp_mat _Rk(N, N);
    arma::sp_mat lambda_mat(B + 1, B + 1);
    arma::colvec lambda_moe = arma::zeros<arma::colvec>(B + 1);
    lambda_moe.subvec(1, B) = lambda;
    lambda_mat.diag() = lambda_moe;
    arma::mat Z_corr = Z_orig;
    unsigned int K = R.n_rows;
    arma::mat W(B + 1, K);

    std::vector<arma::uvec>index;
    // Create index
    std::vector<unsigned>counters;
    arma::vec sizes(sum(Phi, 1));
    // std::cout << sizes << std::endl;
    for (unsigned i = 0; i < sizes.n_elem; i++) {
        arma::uvec a(int(sizes(i)));
        index.push_back(a);
        counters.push_back(0);
    }

    arma::sp_mat::const_iterator it =     Phi.begin();
    arma::sp_mat::const_iterator it_end = Phi.end();
    for(; it != it_end; ++it)
    {
        unsigned int row_idx = it.row();
        unsigned int col_idx = it.col();
        index[row_idx](counters[row_idx]++) = col_idx;
    }

    // Progress p(K, verbose);
    for (unsigned k = 0; k < K; k++) {
        // p.increment();
        // if (Progress::check_abort())
        // return;
        // if (lambda_estimation) {
        // lambda_mat.diag() = find_lambda_cpp(alpha, E.row(k).t());
        // }
        _Rk.diag() = R.row(k);
        arma::sp_mat Phi_Rk = Phi_moe * _Rk;
        arma::mat inv_cov(arma::inv(arma::mat(Phi_Rk * Phi_moe_t + lambda_mat)));

        // Calculate R-scaled PCs once
        arma::mat Z_tmp = Z_orig.each_row() % R.row(k);
        // Generate the betas contribution of the intercept using the data
        // This erases whatever was written before in W
        W = inv_cov.unsafe_col(0) * sum(Z_tmp, 1).t();
        // Calculate betas by calculating each batch contribution
        for(unsigned b=0; b < B; b++) {
            // inv_conv is B+1xB+1 whereas index is B long
            W += inv_cov.unsafe_col(b + 1) * sum(Z_tmp.cols(index[b]), 1).t();
        }
        W.row(0).zeros(); // do not remove the intercept
        Z_corr -= W.t() * Phi_Rk;
    }
    return Z_corr;
}
