#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// ========================= Used for scaleNotCenter ===========================

// [[Rcpp::export]]
arma::sp_mat scaleNotCenter_byRow_rcpp(arma::sp_mat x) {
  arma::uword nrow = x.n_rows, ncol = x.n_cols;
  arma::vec sum_of_squares(nrow);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    sum_of_squares[i.row()] += (*i)*(*i);
  }
  sum_of_squares /= (ncol - 1);
  sum_of_squares = arma::sqrt(sum_of_squares);
  for (arma::sp_mat::iterator i = x.begin(); i != x.end(); ++i) {
    if (sum_of_squares[i.row()] == 0) {
      *i = 0;
    } else {
      *i /= sum_of_squares[i.row()];
    }
  }
  return x;
}

// Do scaleNotCenter on each dataset annotated by `ann` (`n` levels split by
// column), without literally splitting them and merge back.
// x - matrix containing multiple dataset
// ann - interger vector annotating the dataset belongin of each column (cell)
// Must be zero-based
// n - number of datasets
// [[Rcpp::export]]
arma::sp_mat scaleNotCenter_byRow_perDataset_rcpp(arma::sp_mat x,
                                                  const arma::uvec& ann,
                                                  const arma::uword& n) {
  // Count number of cells per group
  arma::uvec n_per_group = arma::zeros<arma::uvec>(n);
  for (unsigned int i=0; i<ann.size(); ++i) n_per_group[ann[i]]++;
  // Calculate sum of squares of each gene in each group
  arma::mat sum_of_squares(x.n_rows, n);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    sum_of_squares(i.row(), ann[i.col()]) += (*i)*(*i);
  }
  for (unsigned int i=0; i<n; ++i) {
    sum_of_squares.col(i) /= n_per_group[ann[i]] - 1;
  }
  sum_of_squares = arma::sqrt(sum_of_squares);
  for (arma::sp_mat::iterator i = x.begin(); i != x.end(); ++i) {
    if (sum_of_squares(i.row(), ann[i.col()]) == 0) {
      *i = 0;
    } else {
      *i /= sum_of_squares(i.row(), ann[i.col()]);
    }
  }
  return x;
}

// [[Rcpp::export]]
NumericVector rowVars_sparse_rcpp(const arma::sp_mat& x,
                                  const NumericVector& means) {
  int nrow = x.n_rows, ncol = x.n_cols;

  NumericVector vars(nrow);
  NumericVector nonzero_vals(nrow);

  for(arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    vars(i.row()) += (*i-means(i.row()))*(*i-means(i.row()));
    nonzero_vals(i.row())++;
  }
  // Add back square mean error for zero elements
  // const_iterator only loops over nonzero elements
  for (int i = 0; i < nrow; ++i) {
    vars(i) += (ncol - nonzero_vals(i))*means(i)*means(i);
    vars(i) /= ncol - 1;
  }
  return vars;
}

// [[Rcpp::export]]
arma::sp_mat rowDivide_rcpp(arma::sp_mat x, const arma::vec& v) {
  if (x.n_rows != v.size()) {
    Rcpp::stop("nrow(x) != length(v)");
  }
  for (arma::sp_mat::iterator i = x.begin(); i != x.end(); ++i) {
    if (v[i.row()] == 0) *i = 0;
    else *i /= v[i.row()];
  }
  return x;
}

// ========================= Used in selectGenes================================

// [[Rcpp::export]]
NumericVector sumSquaredDeviations(arma::sp_mat x, NumericVector means) {
  int nrow = x.n_rows, ncol = x.n_cols;

  NumericVector ssd(nrow);
  NumericVector nonzero_vals(nrow);

  for(arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i)
  {
    ssd(i.row()) += (*i-means(i.row()))*(*i-means(i.row()));
    nonzero_vals(i.row()) += 1;
  }
  // Add back square mean error for zero elements
  // const_iterator only loops over nonzero elements
  for (int i = 0; i < nrow; ++i)
  {
    ssd(i) += (ncol - nonzero_vals(i))*(means(i)*means(i));
  }
  return ssd;
}

//  [[Rcpp::export]]
NumericMatrix denseZScore(NumericMatrix & x, NumericVector m){
  int nrow = x.nrow(), ncol = x.ncol();
  NumericVector v(nrow);
  NumericVector nz(nrow);
  NumericMatrix Z = clone(x);
  NumericVector r(ncol);
  for(int i = 0; i < nrow; i++){
    r = Z(i, _);
    for(NumericVector::iterator j = r.begin(); j != r.end(); ++j){
      v(i) += (*j - m(i)) * (*j - m(i));
    }
    v(i) /= ncol - 1;
    v(i) = sqrt(v(i));
    for (NumericVector::iterator j = r.begin(); j != r.end(); ++j) {
      *j -= m(i);
      *j /= v(i);
    }
  }
  return Z;
}

//  [[Rcpp::export]]
NumericVector rowVarsDense(arma::mat x, arma::vec m) {
  int nrow = x.n_rows;
  int ncol = x.n_cols;
  NumericVector out(nrow);
  for (int i = 0; i < nrow; ++i) {
    arma::rowvec row_i = x.row(i);
    arma::vec diff = row_i.t() - m[i];
    out(i) = arma::accu(arma::square(diff)) / (ncol - 1);
  }
  return out;
}

/* standardize matrix rows using given mean and standard deviation,
 clip values larger than vmax to vmax,
 then return variance for each row */
// [[Rcpp::export(rng = false)]]
NumericVector SparseRowVarStd(arma::sp_mat x,
                              NumericVector mu,
                              NumericVector sd,
                              double vmax){
  x = x.t();
  NumericVector allVars(x.n_cols);
  NumericVector colSum(x.n_cols);
  NumericVector nNonZero(x.n_cols);
  for (arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i) {
    int k = i.col();
    if (sd[k] == 0) continue;
    nNonZero[k] += 1;
    colSum[k] += pow(std::min(vmax, (*i - mu[k]) / sd[k]), 2);
  }
  for (int k=0; k<x.n_cols; ++k) {
    allVars[k] = (colSum[k] + pow((0 - mu[k]) / sd[k], 2) * (x.n_rows - nNonZero[k])) / (x.n_rows - 1);
  }
  return(allVars);
}


// [[Rcpp::export()]]
arma::sp_mat colAggregateSums_sparse(const arma::sp_mat& x,
                                     const arma::uvec& group,
                                     const arma::uword& ngroups) {
    arma::sp_mat out(x.n_rows, ngroups);
    for (arma::sp_mat::const_iterator it = x.begin(); it != x.end(); ++it) {
        out(it.row(), group(it.col())) += *it;
    }
    return out;
}

// // [[Rcpp::export()]]
// arma::sp_mat colAggregateMeans_sparse(const arma::sp_mat& x,
//                                       const arma::uvec& group,
//                                       const arma::uword& ngroups,
//                                       const arma::uvec& groupSizes) {
//     arma::sp_mat out(x.n_rows, ngroups);
//     for (arma::sp_mat::const_iterator it = x.begin(); it != x.end(); ++it) {
//         out(it.row(), group(it.col())) += *it / groupSizes[group(it.col())];
//     }
//     return out;
// }

// [[Rcpp::export()]]
Rcpp::NumericVector sample_cpp(const int x, const int size) {
    arma::uvec rand = arma::randperm<arma::uvec>(x);
    arma::uvec head = rand.head(size) + 1; // 0-base to R's 1-base
    return Rcpp::NumericVector(head.begin(), head.end());
}
