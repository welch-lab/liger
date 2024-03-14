#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

// Algorithm adopted from https://github.com/immunogenomics/presto
// With corrections and removed supports for dense matrices which do not seem
// necessary for LIGER

// X - feature x cell
// output - cell x featureRank
std::list<float> cpp_in_place_rank_mean(arma::vec& v_temp, int idx_begin,
                                        int idx_end) {
    std::list<float> ties;

    if (idx_begin > idx_end) return ties;
    std::vector<pair<float, size_t> > v_sort(idx_end - idx_begin + 1);
    for (size_t i = idx_begin; i <= idx_end; i++) {
        v_sort[i - idx_begin] = make_pair(v_temp[i], i - idx_begin);
    }


    sort(v_sort.begin(), v_sort.end());

    float rank_sum = 0, n = 1;
    size_t i;
    for (i = 1U; i < v_sort.size(); i++) {
        if (v_sort[i].first != v_sort[i - 1].first) {
            // if current val != prev val
            // set prev val to something
            for (unsigned j = 0; j < n; j++) {
                v_temp[v_sort[i - 1 - j].second + idx_begin] =
                    (rank_sum / n) + 1;
            }
            // restart count ranks
            rank_sum = i;
            if (n > 1) ties.push_back(n);
            n = 1;
        } else {
            // if curr val is a tie,
            // don't set anything yet, start computing mean rank
            rank_sum += i;
            n++;
        }
    }
    if (n > 1) ties.push_back(n);
    // set the last element(s)
    for (unsigned j = 0; j < n; j++)
        v_temp[v_sort[i - 1 - j].second + idx_begin] = (rank_sum / n) + 1;

    return ties;
}

// [[Rcpp::export]]
std::vector<std::list<float> > cpp_rank_matrix_dgc(
        arma::vec& x, const arma::vec& p, int nrow, int ncol) {
    vector<list<float> > ties(ncol);
    int n_zero;
    for (int i = 0; i < ncol; i++) {
        n_zero = nrow - (p[i+1] - p[i]);
        if (p[i+1] == p[i]) {
            ties[i].push_back(n_zero);
            continue;
        }
        ties[i] = cpp_in_place_rank_mean(x, p[i], p[i + 1] - 1);
        ties[i].push_back(n_zero);
        x.rows(p[i], p[i + 1] - 1) += n_zero;
    }
    return ties;
}

// [[Rcpp::export]]
arma::mat rowAggregateSum_sparse(arma::sp_mat& X,
                                 const arma::uvec& groups,
                                 unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
        res(groups[it.row()], it.col()) += *it;
    }
    return res;
}

// [[Rcpp::export]]
arma::mat colAggregateSum_sparse(arma::sp_mat& X,
                                 const arma::uvec& groups,
                                 unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
        res(groups[it.col()], it.row()) += *it;
    }
    return res;
}

// Non-zero counting aggregate %%%%
// NNZ - Number of Non-Zero

// [[Rcpp::export]]
arma::mat colNNZAggr_sparse(arma::sp_mat& X,
                            const arma::uvec& groups,
                            unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
        if (*it > 0) res(groups[it.col()], it.row())++;
    }
    return res;
}

// [[Rcpp::export]]
arma::mat rowNNZAggr_sparse(arma::sp_mat& X,
                            const arma::uvec& groups,
                            unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (arma::sp_mat::iterator it = X.begin(); it != X.end(); ++it) {
        if (*it > 0) res(groups[it.row()], it.col())++;
    }
    return res;
}

