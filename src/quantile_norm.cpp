#include <RcppArmadillo.h>
#include <map>
#include <iostream>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
IntegerVector cluster_vote_rcpp(const arma::mat& nn_ranked, IntegerVector clusts) {
  unsigned int k = nn_ranked.n_cols;
  for (unsigned int i = 0; i < nn_ranked.n_rows; ++i){
    std::map<int,int> clust_counts;
    for(int j=0; j<k; ++j) {
      std::pair<std::map<int,int>::iterator,bool> flag = clust_counts.insert(std::pair<int,int>(clusts(nn_ranked(i,j)-1),1));
      if (!flag.second) //failure to insert, element already exists
      {
        ++clust_counts[clusts(nn_ranked(i,j)-1)];
      }
    }
    int max_clust = -1;
    int max_count = 0;
    for (std::map<int,int>::iterator it=clust_counts.begin(); it!=clust_counts.end(); ++it)
    {
      if(it->second > max_count)
      {
        max_clust = it->first;
        max_count = it->second;
      }
    }
    clusts(i) = max_clust;
  }
  return clusts;
}

// [[Rcpp::export]]
IntegerVector max_factor_rcpp(arma::mat H, const arma::uvec& dims_use, bool center=false) {
  // Column wise scaling
  for (unsigned int i = 0; i < H.n_cols; ++i) {
    double colMean = arma::mean(H.col(i));
    double colSdev = 1;
    if (center) {
      H.col(i) -= colMean;
    }
    colSdev = arma::accu(H.col(i)%H.col(i)) / (H.n_rows - 1);
    colSdev = sqrt(colSdev);
    H.col(i) /= colSdev;
  }
  // Find the max index of each row
  IntegerVector clusts(H.n_rows);
  for (unsigned int i = 0; i < H.n_rows; ++i) {
    int max_clust = -1;
    double max_val = 0;
    for (unsigned int k = 0; k < dims_use.size(); ++k) {
      // Assuming argument `dims_use` passed from R is 1-based instead of 0-based
      int j = dims_use[k] - 1;
      if (H(i, j) > max_val) {
        max_clust = j + 1;
        max_val = H(i, j);
      }
    }
    clusts[i] = max_clust;
  }
  return clusts;
}
