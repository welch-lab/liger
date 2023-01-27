#' Louvain algorithm for community detection
#' @description After quantile normalization, users can additionally run the
#' Louvain algorithm for community detection, which is widely used in
#' single-cell analysis and excels at merging small clusters into broad cell
#' classes.
#'
#' While using quantile normalized factor loadings (result from
#' \code{\link{quantile_norm}}) is recommended, this function looks for
#' unnormalized factor loadings (result from \code{\link{optimizeALS}} or
#' \code{\link{online_iNMF}}) when the former is not available.
#' @param object \linkS4class{liger} object. Should run
#' \code{\link{quantile_norm}} before calling.
#' @param k Integer, the maximum number of nearest neighbors to compute. Default
#' \code{20}.
#' @param resolution Numeric, value of the resolution parameter, a larger value
#' results in a larger number of communities with smaller sizes. Default
#' \code{1.0}.
#' @param prune Numeric. Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the stringency of pruning. \code{0} for no pruning,
#' while \code{1} prunes everything. Default \code{1/15}.
#' @param eps Numeric, the error bound of the nearest neighbor search. Default
#' \code{0.1}.
#' @param nRandomStarts Integer number of random starts. Default \code{10}.
#' @param nIterations Integer, maximal number of iterations per random start.
#' Default \code{100}.
#' @param dims.use Indices of factors to use for Louvain clustering. Default
#' \code{NULL} uses all available factors.
#' @param random.seed Seed of the random number generator. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress.
#' Default \code{TRUE}.
#' @return \code{object} with refined cluster assignment updated in
#' \code{louvain_cluster} variable in \code{cell.meta} slot. Can be fetched
#' with \code{object$louvain_cluster}
#' @export
louvainCluster <- function(
        object,
        resolution = 1.0,
        k = 20,
        prune = 1 / 15,
        eps = 0.1,
        nRandomStarts = 10,
        nIterations = 100,
        dims.use = NULL,
        random.seed = 1,
        verbose = TRUE
) {
    H.norm <- getMatrix(object, "H.norm")
    if (is.null(H.norm)) {
        type <- " unnormalized "
        H.norm <- Reduce(cbind, getMatrix(object, "H"))
    } else type <- " quantile normalized "
    if (is.null(H.norm))
        stop("No factor loading ('H.norm' or 'H') found in `object`.")
    # Not transposing when cbind'ing becausing `t(NULL)` causes error
    if (type == " unnormalized ") H.norm <- t(H.norm)
    if (nrow(H.norm) != ncol(object))
        stop("Invalid factor loading matrix. ",
             "Matrix dimensionality does not match to number of cells.")


    edgeOutPath <- paste0('edge_', sub('\\s', '_', Sys.time()), '.txt')
    edgeOutPath <- gsub("-", "", edgeOutPath)
    edgeOutPath <- gsub(":", "", edgeOutPath)

    if (is.null(dims.use)) {
        useFactors <- seq(ncol(H.norm))
    } else {
        useFactors <- dims.use
    }
    if (isTRUE(verbose))
        .log("Louvain clustering on", type, "cell factor loadings...")
    knn <- RANN::nn2(H.norm[, useFactors], k = k, eps = eps)
    snn <- ComputeSNN(knn$nn.idx, prune = prune)
    WriteEdgeFile(snn, edgeOutPath, display_progress = FALSE)
    clusts <- RunModularityClusteringCpp(
        snn,
        modularityFunction = 1,
        resolution = resolution,
        nRandomStarts = nRandomStarts,
        nIterations = nIterations,
        algorithm = 1,
        randomSeed = random.seed,
        printOutput = FALSE,
        edgefilename = edgeOutPath
    )
    names(clusts) <- colnames(object)
    rownames(snn) <- colnames(object)
    colnames(snn) <- colnames(object)
    clusts <-
        GroupSingletons(ids = clusts,
                        SNN = snn,
                        verbose = FALSE)
    cell.meta(object)$louvain_cluster <- factor(clusts)
    unlink(edgeOutPath)
    return(object)
}
