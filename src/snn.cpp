#include <RcppArmadillo.h>
#include <progress.hpp>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export()]]
arma::sp_mat ComputeSNN(arma::umat& nn_idx, double prune) {
    int n = nn_idx.n_rows;
    int k = nn_idx.n_cols;
    arma::umat loc(2, nn_idx.n_elem);
    loc.row(0) = arma::vectorise<arma::umat>(arma::repmat(arma::linspace<arma::umat>(0, n - 1, n), 1, k)).t();
    loc.row(1) = arma::vectorise<arma::umat>(nn_idx).t() - 1;
    arma::vec val = arma::ones<arma::vec>(n*k);
    arma::sp_mat snn(loc, val, n, n);
    snn *= snn.t();
    for (arma::sp_mat::iterator it = snn.begin(); it != snn.end(); ++it) {
        *it /= k + (k - *it);
    }
    snn.for_each([prune](arma::sp_mat::elem_type& val) {
        if (val < prune) val = 0;
    });
    return snn;
}

//[[Rcpp::export]]
void WriteEdgeFile(arma::sp_mat snn, String filename, bool display_progress){
    if (display_progress == true) {
        Rcpp::Rcerr << "Writing SNN as edge file" << std::endl;
    }
    // Write out lower triangle
    std::ofstream output;
    output.open(filename);
    Progress p(snn.n_elem, display_progress);
    for (arma::sp_mat::const_iterator it = snn.begin(); it != snn.end(); ++it){
        p.increment();
        if(it.col() >= it.row()){
            continue;
        }
        output << std::setprecision(15) << it.col() << "\t" << it.row() << "\t" << *it << "\n";
    }
    output.close();
}

//[[Rcpp::export]]
arma::sp_mat DirectSNNToFile(arma::umat& nn_ranked,
                                            double prune, bool display_progress,
                                            String filename) {
    arma::sp_mat SNN = ComputeSNN(nn_ranked, prune);
    WriteEdgeFile(SNN, filename, display_progress);
    return SNN;
}
