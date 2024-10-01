#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat normalize_byCol_dense_rcpp(arma::mat x) {
  arma::vec sums = arma::sum(x, 0).t();
  for (arma::uword i = 0; i < x.n_cols; ++i) {
    if (sums[i] == 0) {
      x.col(i).fill(0);
    } else {
      x.col(i) /= sums[i];
    }
  }
  return x;
}

// Different from the one above, this function supports L-n normalization
// while the one above only does sum normalization, not even necessarily L-1
// [[Rcpp::export()]]
arma::mat colNormalize_dense_cpp(arma::mat& x, const arma::uword L) {
  arma::mat result(x);
  for (int j = 0; j < x.n_cols; ++j) {
    double norm = arma::norm(x.col(j), L);
    if (norm == 0) {
      continue;
    }
    for (int i = 0; i < x.n_rows; ++i) {
      result(i, j) /= norm;
    }
  }
  return result;
}

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

// [[Rcpp::export]]
arma::mat safe_scale(arma::mat x, bool center, bool scale) {
  arma::vec means = arma::mean(x, 0).t();
  if (center) {
    for (arma::uword i = 0; i < x.n_cols; ++i) {
      x.col(i) -= means[i];
    }
  }

  if (scale) {
    arma::vec sds = arma::zeros<arma::vec>(x.n_cols);
    for (arma::uword i = 0; i < x.n_cols; ++i) {
      for (arma::uword j = 0; j < x.n_rows; ++j) {
        sds[i] += x(j, i) * x(j, i);
      }
      sds[i] = std::sqrt(sds[i] / (x.n_rows - 1));

      if (sds[i] == 0) {
        x.col(i).fill(0);
      } else {
        x.col(i) /= sds[i];
      }
    }
  }
  return x;
}


// [[Rcpp::export]]
arma::mat scaleNotCenter_byCol_dense_rcpp(
  arma::mat x
) {
  arma::uword nrow = x.n_rows, ncol = x.n_cols;
  arma::vec sum_of_squares(ncol);
  for (arma::uword i = 0; i < nrow; ++i) {
    for (arma::uword j = 0; j < ncol; ++j) {
      sum_of_squares[j] += x(i, j) * x(i, j);
    }
  }
  sum_of_squares /= (nrow - 1);
  sum_of_squares = arma::sqrt(sum_of_squares);
  for (arma::uword i = 0; i < nrow; ++i) {
    for (arma::uword j = 0; j < ncol; ++j) {
      if (sum_of_squares[j] == 0) {
        x(i, j) = 0;
      } else {
        x(i, j) /= sum_of_squares[j];
      }
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

// [[Rcpp::export()]]
Rcpp::NumericVector sample_cpp(const int x, const int size) {
    arma::uvec rand = arma::randperm<arma::uvec>(x);
    arma::uvec head = rand.head(size) + 1; // 0-base to R's 1-base
    return Rcpp::NumericVector(head.begin(), head.end());
}


// psdBulk - the pseudo-bulk matrix to be updated
// sparseRaw - the raw expression data to be collapsed
// featureIdx - integer vector with NAs, length equals to nrow(sparseRaw),
//   indicating which row of `psdBulk` the current `it.row()` should be added
//   to. Zero-based before invoking
// repIdx - Integer vector with NAs, length equals to nrow(sparseRaw),
//   indicating which col of `psdBulk` the current `it.col()` should be added
//   to/ Zero-based before invoking
// [[Rcpp::export()]]
void updatePseudoBulkRcpp(
  Rcpp::NumericMatrix& psdBulk,
  const arma::sp_mat& sparseRaw,
  const Rcpp::IntegerVector& featureIdx,
  const Rcpp::IntegerVector& repIdx
) {
  for (arma::sp_mat::const_iterator it = sparseRaw.begin(); it != sparseRaw.end(); ++it) {
    if (featureIdx[it.row()] != NA_INTEGER && repIdx[it.col()] != NA_INTEGER) {
      psdBulk(featureIdx[it.row()], repIdx[it.col()]) += *it;
    }
  }
}

// out - dense matrix of n feature rows and 2 cols
// sparseRaw - the raw expression data to be counted
// featureIdx - integer vector with NAs, length equals to nrow(sparseRaw),
//   indicating which row of `out` the current `it.row()` should be incremented.
// groupVar - integer vector with NAs, length equals to ncol(sparseRaw),
//   indicating which col of `out` the current `it.col()` should be incremented.
// [[Rcpp::export()]]
void updateNCellExprRcpp(
  Rcpp::NumericMatrix& out,
  const arma::sp_mat& sparseRaw,
  const Rcpp::IntegerVector& featureIdx,
  const Rcpp::IntegerVector& groupVar
) {
  for (arma::sp_mat::const_iterator it = sparseRaw.begin(); it != sparseRaw.end(); ++it) {
    if (featureIdx[it.row()] != NA_INTEGER && groupVar[it.col()] != NA_INTEGER) {
      out(featureIdx[it.row()], groupVar[it.col()])++;
    }
  }
}
