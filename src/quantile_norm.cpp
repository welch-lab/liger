#include <RcppEigen.h>
#include <map>
#include <iostream>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
IntegerVector cluster_vote(const Eigen::MatrixXd & nn_ranked,IntegerVector clusts) {
  int k = nn_ranked.cols();
  for(int i=0; i<nn_ranked.rows(); ++i){
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

//Note: code modified from Seurat package https://github.com/satijalab/seurat/blob/master/src/data_manipulation.cpp
// [[Rcpp::export]]
Eigen::MatrixXd scale_columns_fast(Eigen::MatrixXd mat, bool scale = true, bool center = true){
  for(int i=0; i < mat.cols(); ++i){
    Eigen::ArrayXd c = mat.col(i).array();
    double colMean = c.mean();
    double colSdev = 1;
    if(scale == true){
      if(center == true){
        colSdev = sqrt((c - colMean).square().sum() / (mat.rows() - 1));
      }
      else{
        colSdev = sqrt(c.square().sum() / (mat.rows() - 1));
      }
    }
    if(center == false){
      colMean = 0;
    }
    mat.col(i) = (c - colMean) / colSdev;
  }
  return mat;
}

// [[Rcpp::export]]
IntegerVector max_factor(Eigen::MatrixXd H, IntegerVector dims_use, bool center_cols=false) {
  H = scale_columns_fast(H,true,center_cols);
  IntegerVector clusts(H.rows());
  for(int i=0; i < H.rows(); ++i){
    int max_clust = -1;
    double max_val = 0;
    for(int k=0; k < dims_use.length(); ++k) {
      int j = dims_use(k)-1;
      if(H(i,j) > max_val)
      {
        max_clust = j+1;
        max_val = H(i,j);
      }
      
      clusts(i) = max_clust;  
    }
  }
  return clusts;
}
