#pragma once

#include <Rcpp.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rconfig.h>
#include <R_ext/Rdynload.h>
#include <vector>

#ifdef HAVE_VISIBILITY_ATTRIBUTE
  # define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
  # define attribute_hidden
#endif

#ifdef __cplusplus
extern "C" {
#endif


// [[Rcpp::depends(leidenAlg)]]

//' Refer to the R function find_partition()
//' For notes of the graph object, refer to https://igraph.org/c/doc/igraph-Basic.html
//'
//' @param edgelist The graph edge list
//' @param edgelist_length integer The length of the graph edge list
//' @param num_vertices integer The number of vertices in the graph
//' @param direction boolean Whether the graph is directed or undirected
//' @param edge_weights Vector of edge weights. In weighted graphs, a real number is assigned to each (directed or undirected) edge. For an unweighted graph, this is set to 1. Refer to igraph, weighted graphs.
//' @param resolution Integer resoluiton parameter controlling communities detected (default=1.0) Higher resolutions lead to more communities, while lower resolutions lead to fewer communities.
//' @param niter Number of iterations that the algorithm should be run for (default=2)
//' @return A vector of membership values
// [[Rcpp::export]]
std::vector<size_t> find_partition_rcpp(std::vector<int> &edgelist, int edgelist_length, int num_vertices, bool direction, std::vector<double> &edge_weights, double resolution, int niter)
{
    static SEXP (*fun)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP) = (SEXP(*)(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP))R_GetCCallable("leidenAlg", "find_partition_rcpp");
    return fun(edgelist, edgelist_length, num_vertices, direction, edge_weights, resolution, niter);
}

#ifdef __cplusplus
}
#endif
