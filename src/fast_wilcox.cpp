#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Codes from Presto (https://github.com/immunogenomics/presto)

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;


// [[Rcpp::export]]
arma::mat cpp_sumGroups_dgc(const arma::vec& x, const arma::uvec& p, 
                            const arma::vec& i, unsigned ncol, 
                            const arma::uvec& groups, unsigned ngroups) {
    // Here, columns are genes
    arma::mat res = arma::zeros<arma::mat>(ngroups, ncol);
    for (unsigned c = 0; c < ncol; c++) {
        for (unsigned j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            // group_map gives the group num of that row
            res(groups[i[j]], c) += x[j];
        }
    }    
    return res;
}

// [[Rcpp::export]]
arma::mat cpp_sumGroups_dgc_T(const arma::vec& x, const arma::vec& p, 
                              const arma::vec& i, int ncol, int nrow, 
                              const arma::uvec& groups, int ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, nrow);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            // group_map gives the group num of that row
            res(groups[c], i[j]) += x[j];
        }
    }    
    return res;
}

// [[Rcpp::export]]
arma::mat cpp_sumGroups_dense(const arma::mat& X, const arma::uvec& groups,
                              unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (unsigned r = 0; r < X.n_rows; r++) {
        res.row(groups[r]) += sum(X.row(r), 0);
    }    
    return res;
}

// [[Rcpp::export]]
arma::mat cpp_sumGroups_dense_T(const arma::mat& X, const arma::uvec& groups,
                                unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (unsigned c = 0; c < X.n_cols; c++) {
        res.row(groups[c]) += sum(X.col(c), 1).t();
    }    
    return res;
}


// [[Rcpp::export]]
arma::mat cpp_nnzeroGroups_dense(const arma::mat& X, const arma::uvec& groups,
                                 unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_cols);
    for (unsigned c = 0; c < X.n_cols; c++) {
        for (unsigned r = 0; r < X.n_rows; r++) {
            if (X(r, c) != 0)
                res(groups[r], c)++;
        }
    }    
    return res;
}

// [[Rcpp::export]]
arma::mat cpp_nnzeroGroups_dense_T(const arma::mat& X,
                                   const arma::uvec& groups, 
                                   unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, X.n_rows);
    for (unsigned c = 0; c < X.n_cols; c++) {
        for (unsigned r = 0; r < X.n_rows; r++) {
            if (X(r, c) != 0)
                res(groups[c], r)++;
//             res.row(groups[c]) += sum(X.col(c), 1).t();
        }
    }    
    return res;
}



// [[Rcpp::export]]
arma::mat cpp_nnzeroGroups_dgc(const arma::uvec& p, const arma::vec& i,
                               unsigned ncol, const arma::uvec& groups,
                               unsigned ngroups) {
    arma::mat res = arma::zeros<arma::mat>(ngroups, ncol);
    for (unsigned c = 0; c < ncol; c++) {
        for (unsigned j = p[c]; j < p[c + 1]; j++) {
            res(groups[i[j]], c)++;
        }
    }    
    return res;
}


// [[Rcpp::export]]
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
        if (p[i+1] == p[i]) continue;
        n_zero = nrow - (p[i+1] - p[i]);
        ties[i] = cpp_in_place_rank_mean(x, p[i], p[i + 1] - 1);
        ties[i].push_back(n_zero);
        x.rows(p[i], p[i + 1] - 1) += n_zero;
    }
    return ties;
}


// [[Rcpp::export]]
Rcpp::List cpp_rank_matrix_dense(arma::mat& X) {
    // sizes of tied groups
    arma::inplace_trans(X);
    vector<list<float> > ties(X.n_cols);
    
    std::vector<pair<float, size_t> > v_sort(X.n_rows);
    for (unsigned c = 0; c < X.n_cols; c++) {
        for (size_t i = 0; i < X.n_rows; i++) {
            v_sort[i] = make_pair(X.col(c)[i], i);
        }
        sort(v_sort.begin(), v_sort.end());

        float rank_sum = 0, n = 1;
        size_t i;
        for (i = 1U; i < v_sort.size(); i++) {
            if (v_sort[i].first != v_sort[i - 1].first) {
                // if current val != prev val
                // set prev val to something
                for (unsigned j = 0; j < n; j++) {
                    X.col(c)[v_sort[i - 1 - j].second] = (rank_sum / n) + 1;  
                }            
                // restart count ranks
                rank_sum = i;
                if (n > 1) ties[c].push_back(n);
                n = 1;
            } else {
                // if curr val is a tie, 
                // don't set anything yet, start computing mean rank
                rank_sum += i;
                n++;
            }
        }
        // set the last element(s)
        for (unsigned j = 0; j < n; j++)
            X.col(c)[v_sort[i - 1 - j].second] = (rank_sum / n) + 1;
    }
    return Rcpp::List::create(Named("X_ranked") = X, Named("ties") = ties);
}




// [[Rcpp::export]]
arma::mat cpp_nnzeroGroups_dgc_T(const arma::vec& p, const arma::vec& i, 
                                 int ncol, int nrow, const arma::uvec& groups,
                                 int ngroups) {
    // Here, columns are samples
    arma::mat res = arma::zeros<arma::mat>(ngroups, nrow);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            // group_map gives the group num of that row
            res(groups[c], i[j])++;
        }
    }    
    return res;
}
