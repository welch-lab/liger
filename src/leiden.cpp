#include <RcppArmadillo.h>
#include <R.h>
#include <Rdefines.h>
#include <vector>
#include <progress.hpp>
#include "igraph.h"
#include "igraph_constructors.h"
#include "GraphHelper.h"
#include "Optimiser.h"
#include "RBERVertexPartition.h"
#include "RBConfigurationVertexPartition.h"
#include "igraph_constructors.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppProgress)]]

// Debug Mode implies checking assertions.
#if defined(_GLIBCXX_ASSERTIONS)
# define _GLIBCXX_ASSERTIONS 0
#endif

// a wrapper for the Leidgen algorithm implementation (https://github.com/vtraag/leidenalg)

// https://stackoverflow.com/questions/10250438/using-stdvector-with-igraph
//
void Stl_To_Igraph_vector_t(std::vector<int>& vectR, igraph_vector_t* v) {
    size_t n = vectR.size();

    /* Make sure that there is enough space for the items in v */
    igraph_vector_resize(v, n);

    /* Copy all the items */
    for (size_t i = 0; i < n; i++) {
        VECTOR(*v)[i] = vectR[i];
    }
}


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
std::vector<size_t> leidenalg_find_partition_rcpp(std::vector<int>& edgelist, int edgelist_length, int num_vertices, bool direction, std::vector<double>& edge_weights, double resolution=1.0, int nStart = 10, int nIter = 2, int seed = 123, bool verbose = true) {
  if (verbose) {
    std::cerr << "Running Leiden Algorithm" << std::endl <<
      "Resolution    = " << resolution << std::endl <<
      "Max iteration = " << nIter << std::endl <<
      "Random starts = " << nStart << std::endl;
  }
  igraph_t g;
  igraph_vector_t edges;

  // initialize igraph_vector_t
  igraph_vector_init(&edges, edgelist_length);

  // questionable attempt to convert 'std::vector<int>' to 'igraph_vector_t' (not 'igraph_vector_int_t' as documented by 'igraph_create()')
  Stl_To_Igraph_vector_t(edgelist, &edges);

  igraph_create(&g, &edges, num_vertices, direction);

  Graph og(&g, edge_weights);

  std::vector<size_t> bestCluster;
  double maxQuality = -std::numeric_limits<double>::infinity();

  Progress pb(nStart, verbose);
  for (int i = 0; i < nStart; i++)
  {
    seed++;
    Optimiser o(seed);
    RBConfigurationVertexPartition p(&og,resolution);
    double val = 1;
    int iter = 0;
    while (val > 0 && (iter < nIter || nIter < 0)) {
      val = o.optimise_partition(&p);
      iter++;
    }
    double q = p.quality(resolution);
    if (q > maxQuality) {
      maxQuality = q;
      bestCluster = p.membership();
    }
    pb.increment();
  }


  // destroy the igraph_t 'g' and the igraph_vector_t 'edges'
  // https://igraph.org/c/html/latest/igraph-Data-structures.html#igraph_vector_destroy
  // https://igraph.org/c/html/latest/igraph-Generators.html#igraph_create
  igraph_destroy(&g);
  igraph_vector_destroy(&edges);

  return(bestCluster);
  //return(igraph_ecount(&g));
}

