#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::sp_mat scaleNotCenterFast(arma::sp_mat x) {
  int nrow = x.n_rows, ncol = x.n_cols;
  
  NumericVector sum_of_squares(ncol);
  for(arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i)
  {
    sum_of_squares(i.col()) += (*i)*(*i);
  }
  for (int i = 0; i < ncol; ++i)
  {
    sum_of_squares(i) = sqrt(sum_of_squares(i)/(nrow-1));
  }
  for(arma::sp_mat::iterator i = x.begin(); i != x.end(); ++i)
  {
    *i /= sum_of_squares(i.col());
  }
  return x;
}

// [[Rcpp::export]]
NumericVector rowMeansFast(arma::sp_mat x) {
  int nrow = x.n_rows, ncol = x.n_cols;
  
  NumericVector means(nrow);
  for (int i = 0; i < nrow; ++i)
  {
    means(i) = 0;
  }
  for(arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i)
  {
    means(i.row()) += *i;
  }
  for (int i = 0; i < nrow; ++i)
  {
    means(i) /= ncol;
  }
  return means;
}

// [[Rcpp::export]]
NumericVector rowVarsFast(arma::sp_mat x, NumericVector means) {
  int nrow = x.n_rows, ncol = x.n_cols;
  
  NumericVector vars(nrow);
  NumericVector nonzero_vals(nrow);
  
  for(arma::sp_mat::const_iterator i = x.begin(); i != x.end(); ++i)
  {
    vars(i.row()) += (*i-means(i.row()))*(*i-means(i.row())); 
    nonzero_vals(i.row()) += 1;
  }
  // Add back square mean error for zero elements
  // const_iterator only loops over nonzero elements 
  for (int i = 0; i < nrow; ++i)
  {
    vars(i) += (ncol - nonzero_vals(i))*(means(i)*means(i));
    vars(i) /= (ncol-1);
  }
  return vars;
}

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

