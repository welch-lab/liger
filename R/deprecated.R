#' @importFrom Matrix colSums rowSums rowMeans t sparseMatrix
NULL

# These are deprecated functions likely to be removed in future versions.
# Documentation for these functions is incomplete.

#' Quantile align (normalize) factor loadings
#' 
#' This is a deprecated function. Calling 'quantile_norm' instead.
#'
#' This process builds a shared factor neighborhood graph to jointly cluster cells, then quantile
#' normalizes corresponding clusters.
#'
#' The first step, building the shared factor neighborhood graph, is performed in SNF(), and
#' produces a graph representation where edge weights between cells (across all datasets)
#' correspond to their similarity in the shared factor neighborhood space. An important parameter
#' here is knn_k, the number of neighbors used to build the shared factor space (see SNF()). Afterwards,
#' modularity-based community detection is performed on this graph (Louvain clustering) in order
#' to identify shared clusters across datasets. The method was first developed by Waltman and van Eck
#' (2013) and source code is available at http://www.ludowaltman.nl/slm/. The most important parameter
#' here is resolution, which corresponds to the number of communities detected.
#'
#' Next we perform quantile alignment for each dataset, factor, and cluster (by
#' stretching/compressing datasets' quantiles to better match those of the reference dataset). These
#' aligned factor loadings are combined into a single matrix and returned as H.norm.
#'
#' @param object \code{liger} object. Should run optimizeALS before calling.
#' @param knn_k Number of nearest neighbors for within-dataset knn graph (default 20).
#' @param k2 Horizon parameter for shared nearest factor graph. Distances to all but the k2 nearest
#'   neighbors are set to 0 (cuts down on memory usage for very large graphs). (default 500)
#' @param prune.thresh Minimum allowed edge weight. Any edges below this are removed (given weight
#'  0) (default 0.2)
#' @param ref_dataset Name of dataset to use as a "reference" for normalization. By default,
#'   the dataset with the largest number of cells is used.
#' @param min_cells Minimum number of cells to consider a cluster shared across datasets (default 2)
#' @param quantiles Number of quantiles to use for quantile normalization (default 50).
#' @param nstart Number of times to perform Louvain community detection with different random
#'   starts (default 10).
#' @param resolution Controls the number of communities detected. Higher resolution -> more
#'   communities. (default 1)
#' @param dims.use Indices of factors to use for shared nearest factor determination (default
#'   1:ncol(H[[1]])).
#' @param dist.use Distance metric to use in calculating nearest neighbors (default "CR").
#' @param center Centers the data when scaling factors (useful for less sparse modalities like
#'   methylation data). (default FALSE)
#' @param small.clust.thresh Extracts small clusters loading highly on single factor with fewer
#'   cells than this before regular alignment (default 0 -- no small cluster extraction).
#' @param id.number Number to use for identifying edge file (when running in parallel)
#'   (generates random value by default).
#' @param print.mod Print modularity output from clustering algorithm (default FALSE).
#' @param print.align.summary Print summary of clusters which did not align normally (default FALSE).
#'
#' @return \code{liger} object with H.norm and cluster slots set.
#' @export
#'
#' @examples
#' \dontrun{
#' # liger object, factorization complete
#' ligerex
#' # do basic quantile alignment
#' ligerex <- quantileAlignSNF(ligerex)
#' # higher resolution for more clusters (note that SNF is conserved)
#' ligerex <- quantileAlignSNF(ligerex, resolution = 1.2)
#' # change knn_k for more fine-grained local clustering
#' ligerex <- quantileAlignSNF(ligerex, knn_k = 15, resolution = 1.2)
#' }
#'
quantileAlignSNF <- function(
                             object,
                             knn_k = 20,
                             k2 = 500,
                             prune.thresh = 0.2,
                             ref_dataset = NULL,
                             min_cells = 20,
                             quantiles = 50,
                             nstart = 10,
                             resolution = 1,
                             dims.use = 1:ncol(x = object@H[[1]]),
                             dist.use = "CR",
                             center = FALSE,
                             small.clust.thresh = 0,
                             id.number = NULL,
                             print.mod = FALSE,
                             print.align.summary = FALSE) {
  .Deprecated(
    msg = paste(
      "This is a deprecated function. Calling 'quantile_norm' instead.",
      "Note that not all parameters can be passed to 'quantile_norm'.",
      "It's suggested to run 'louvainCluster' subsequently as well."
    )
  )
  return(quantile_norm(object,
    quantiles = quantiles,
    ref_dataset = ref_dataset,
    min_cells = min_cells,
    knn_k = knn_k,
    dims.use = dims.use,
    do.center = center,
    max_sample = 1000,
    eps = 0.9,
    refine.knn = T
  ))
}
