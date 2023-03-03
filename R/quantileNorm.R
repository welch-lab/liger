#' @rdname quantileNorm
#' @export
setGeneric(
    "quantileNorm",
    function(
        object,
        quantiles = 50,
        reference = NULL,
        minCells = 20,
        nNeighbors = 20,
        useDims = NULL,
        center = FALSE,
        maxSample = 1000,
        eps = 0.9,
        refineKNN = TRUE,
        clusterName = "H.norm_cluster",
        seed = 1,
        verbose = TRUE,
        ref_dataset = reference,
        min_cells = minCells,
        knn_k = nNeighbors,
        dims.use = useDims,
        do.center = center,
        max_sample = maxSample,
        refine.knn = refineKNN,
        rand.seed = seed
    ) standardGeneric("quantileNorm")
)

#' Quantile Align (Normalize) Factor Loadings
#' @description This process builds a shared factor neighborhood graph to
#' jointly cluster cells, then quantile normalizes corresponding clusters.
#'
#' The first step, building the shared factor neighborhood graph, is performed
#' in SNF(), and produces a graph representation where edge weights between
#' cells (across all datasets) correspond to their similarity in the shared
#' factor neighborhood space. An important parameter here is \code{knn_k}, the
#' number of neighbors used to build the shared factor space.
#'
#' Next we perform quantile alignment for each dataset, factor, and cluster (by
#' stretching/compressing datasets' quantiles to better match those of the
#' reference dataset). These aligned factor loadings are combined into a single
#' matrix and returned as "H.norm".
#'
#' This function is a replacement of \code{quantile_norm()}, which will be
#' deprecated. See \code{help("quantile_norm-deprecated")}.
#' @param object A \linkS4class{liger} object with \code{\link{optimizeALS}} run
#' in advance. Or a list object of \eqn{H} matrices, all in \eqn{k * nCell}
#' dimensionality, the \code{H} entry of the output of \code{\link{optimizeALS}}
#' list method is ready for use (i.e. use \code{out$H}).
#' @param nNeighbors Number of nearest neighbors for within-dataset knn graph.
#' Default \code{20}.
#' @param reference Character, numeric or logical selection of one dataset, out
#' of all available datasets in \code{object}, to use as a "reference" for
#' normalization. Default \code{NULL} use the dataset with the largest number of
#' cells.
#' @param minCells Minimum number of cells to consider a cluster shared across
#' datasets. Default \code{20}.
#' @param quantiles Number of quantiles to use for quantile normalization.
#' Default \code{50}.
#' @param eps The error bound of the nearest neighbor search. Lower values give
#' more accurate nearest neighbor graphs but take much longer to compute.
#' Default \code{0.9}.
#' @param useDims Indices of factors to use for shared nearest factor
#' determination. Default \code{NULL} uses all factors.
#' @param center Whether to center the data when scaling factors (centering
#' could be useful for less sparse modalities like methylation data). Default
#' \code{FALSE}.
#' @param maxSample Maximum number of cells used for quantile normalization of
#' each cluster and factor. Default \code{1000}.
#' @param refineKNN whether to increase robustness of cluster assignments using
#' KNN graph. Default \code{TRUE}.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{TRUE}.
#' @return For method on \linkS4class{liger} object, the object will be returned
#' with the \code{H.norm} slot and \code{"H.norm_cluster"} variable in
#' \code{cell.meta} slot updated. For method on list, a result list with entries
#' of \code{H.norm} and \code{clusters} will be returned.
#' @rdname quantileNorm
#' @export
setMethod(
    "quantileNorm",
    signature(object = "liger"),
    function(
        object,
        quantiles = 50,
        reference = NULL,
        minCells = 20,
        nNeighbors = 20,
        useDims = NULL,
        center = FALSE,
        maxSample = 1000,
        eps = 0.9,
        refineKNN = TRUE,
        clusterName = "H.norm_cluster",
        seed = 1,
        verbose = TRUE
    ) {
        .checkValidFactorResult(object, useDatasets = names(object))
        if (is.null(reference)) {
            # If ref_dataset not given, set the one with the largest number of
            # cells as reference.
            # Should not produce intermediate variable here because it'll be
            # recorded as a environment parameter in object@commands
            reference <- names(which.max(sapply(datasets(object), ncol)))
        } else {
            reference <- .checkUseDatasets(object, useDatasets = reference)
            if (length(reference) != 1)
                stop("Should specify only one reference dataset.")
        }
        object <- recordCommand(object, dependencies = "RANN")
        out <- quantileNorm(
            object = getMatrix(object, "H"),
            quantiles = quantiles,
            reference = reference,
            minCells = minCells,
            nNeighbors = nNeighbors,
            useDims = useDims,
            center = center,
            maxSample = maxSample,
            eps = eps,
            refineKNN = refineKNN,
            seed = seed
        )
        object@H.norm <- out$H.norm
        cell.meta(object, clusterName, check = FALSE) <- out$clusters
        return(object)
    })

#' @rdname quantileNorm
#' @export
setMethod(
    "quantileNorm",
    signature(object = "list"),
    function(
        object,
        quantiles = 50,
        reference = NULL,
        minCells = 20,
        nNeighbors = 20,
        useDims = NULL,
        center = FALSE,
        maxSample = 1000,
        eps = 0.9,
        refineKNN = TRUE,
        seed = 1,
        verbose = TRUE
    ) {
        set.seed(seed)
        if (!is.list(object))
            stop("`object` should be a list of matrices")
        if (is.null(names(object))) {
            stop("`object` should be a named list of matrices")
        }
        if (!all(sapply(object, is.matrix)))
            stop("All values in `object` must be a matrix")
        if (!all(sapply(object, nrow) == nrow(object[[1]])))
            stop("All matrices must have the same number of rows (k factors).")
        if (is.character(reference)) {
            if (length(reference) != 1 || !reference %in% names(object))
                stop("Should specify one existing dataset as reference. ",
                     "(character `reference` wrong length or not found)")
        } else if (is.numeric(reference)) {
            if (length(reference) != 1 || reference > length(object))
                stop("Should specify one dataset within the range. ",
                     "(numeric `reference` wrong length or out of bound)")
        } else if (is.logical(reference)) {
            if (length(reference) != length(object) || sum(reference) != 1)
                stop("Should specify one dataset within the range. ",
                     "(logical `reference` wrong length or ",
                     "too many selection)")
        } else {
            stop("Unable to understand `reference`. See `?quantile_norm`.")
        }

        if (is.null(useDims)) {
            useDims <- seq(nrow(object[[1]]))
        }
        # Transposing all H to cell x k
        Hs <- lapply(object, t)
        # fast max factor assignment with Rcpp code
        clusters <- lapply(Hs, max_factor,
                         dims_use = useDims, center_cols = center)
        # TODO: Dumb to mess up factor with characters of numbers
        # Change to numeric in the future. Need to reproduce result for now.
        # clusterAssign <- factor(unlist(clusters))
        clusterAssign <- as.factor(unlist(lapply(clusters, as.character)))
        names(clusterAssign) <- unlist(lapply(Hs, rownames))
        # increase robustness of cluster assignments using knn graph
        if (isTRUE(refineKNN)) {
            for (i in seq_along(Hs)) {
                clustsH <- clusterAssign[rownames(Hs[[i]])]
                H_knn <- RANN::nn2(Hs[[i]], eps = eps, k = nNeighbors,
                                   searchtype = "standard")
                #print(clustsH[1:10])
                # Rcpp method cluster_vote
                newClusts <- cluster_vote(H_knn$nn.idx, clustsH)
                clusterAssign[rownames(Hs[[i]])] <- newClusts
            }
        }
        # TODO: This line need to be removed
        clusters <- lapply(Hs, function(H) clusterAssign[rownames(H)])
        dims <- ncol(Hs[[reference]])
        nClusters <- dims
        for (d in seq_along(Hs)) {
            for (c in seq(nClusters)) {
                # cells 2: cells in d'th dataset belonging to cluster c
                cellIdx2 <- clusters[[d]] == c
                # cells 1: cells in ref dataset belong to cluster c
                cellIdx1 <- clusters[[reference]] == c
                nCells2 <- sum(cellIdx2)
                nCells1 <- sum(cellIdx1)
                if (nCells1 < minCells || nCells2 < minCells) next
                for (k in seq(dims)) {
                    if (nCells2 == 1) {
                        Hs[[d]][cellIdx2, k] <-
                            mean(Hs[[reference]][cellIdx1, k])
                        next
                    }
                    q2 <- stats::quantile(sample(Hs[[d]][cellIdx2, k],
                                          min(nCells2, maxSample)),
                                   seq(0, 1, by = 1 / quantiles))
                    q1 <- stats::quantile(sample(Hs[[reference]][cellIdx1, k],
                                          min(nCells1, maxSample)),
                                   seq(0, 1, by = 1 / quantiles))
                    if (sum(q1) == 0 | sum(q2) == 0 |
                        length(unique(q1)) < 2 | length(unique(q2)) < 2) {
                        newValue <- rep(0, nCells2)
                    } else {
                        warp_func <- withCallingHandlers(
                                stats::approxfun(q2, q1, rule = 2),
                                warning = function(w) {
                                    invokeRestart("muffleWarning")
                                }
                        )
                        newValue <- warp_func(Hs[[d]][cellIdx2, k])
                    }
                    Hs[[d]][cellIdx2, k] <- newValue
                }
            }
        }
        return(list('H.norm' = Reduce(rbind, Hs),
                    'clusters' = clusterAssign))
    })






#' [Deprecated] Quantile Align (Normalize) Factor Loadings
#' @description This process builds a shared factor neighborhood graph to
#' jointly cluster cells, then quantile normalizes corresponding clusters.
#'
#' The first step, building the shared factor neighborhood graph, is performed
#' in SNF(), and produces a graph representation where edge weights between
#' cells (across all datasets) correspond to their similarity in the shared
#' factor neighborhood space. An important parameter here is \code{knn_k}, the
#' number of neighbors used to build the shared factor space.
#'
#' Next we perform quantile alignment for each dataset, factor, and cluster (by
#' stretching/compressing datasets' quantiles to better match those of the
#' reference dataset). These aligned factor loadings are combined into a single
#' matrix and returned as "H.norm".
#'
#' This function will be deprecated. Please use \code{\link{quantileNorm}}
#' instead.
#' @param object A \linkS4class{liger} object with \code{\link{optimizeALS}} run
#' in advance. Or a list object of \eqn{H} matrices, all in \eqn{k * nCell}
#' dimensionality, the \code{H} entry of the output of \code{\link{optimizeALS}}
#' list method is ready for use (i.e. use \code{out$H}).
#' @param knn_k Number of nearest neighbors for within-dataset knn graph.
#' Default \code{20}.
#' @param ref_dataset Character, numeric or logical selection of one dataset,
#' out of all available datasets in \code{object}, to use as a "reference" for
#' normalization. Default \code{NULL} use the dataset with the largest number of
#' cells.
#' @param min_cells Minimum number of cells to consider a cluster shared across
#' datasets. Default \code{20}.
#' @param quantiles Number of quantiles to use for quantile normalization.
#' Default \code{50}.
#' @param eps The error bound of the nearest neighbor search. Lower values give
#' more accurate nearest neighbor graphs but take much longer to compute.
#' Default \code{0.9}.
#' @param dims.use Indices of factors to use for shared nearest factor
#' determination. Default \code{NULL} uses all factors.
#' @param do.center Whether to center the data when scaling factors (centering
#' could be useful for less sparse modalities like methylation data). Default
#' \code{FALSE}.
#' @param max_sample Maximum number of cells used for quantile normalization of
#' each cluster and factor. Default \code{1000}.
#' @param refine.knn whether to increase robustness of cluster assignments using
#' KNN graph. Default \code{TRUE}.
#' @param rand.seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{TRUE}.
#' @return For method on \linkS4class{liger} object, the object will be returned
#' with the \code{H.norm} slot and \code{"H.norm_cluster"} variable in
#' \code{cell.meta} slot updated. For method on list, a result list with entries
#' of \code{H.norm} and \code{clusters} will be returned.
#' @name quantile_norm-deprecated
#' @seealso \code{\link{rliger-deprecated}}
NULL

#' @rdname rliger-deprecated
#' @section \code{quantile_norm}:
#' For \code{quantile_norm}, use \code{\link{quantileNorm}}.
#' @export
setGeneric(
    "quantile_norm",
    function(
        object,
        quantiles = 50,
        ref_dataset = NULL,
        min_cells = 20,
        knn_k = 20,
        dims.use = NULL,
        do.center = FALSE,
        max_sample = 1000,
        eps = 0.9,
        refine.knn = TRUE,
        clusterName = "H.norm_cluster",
        rand.seed = 1,
        verbose = TRUE
    ) standardGeneric("quantile_norm")
)

#' @rdname rliger-deprecated
#' @section \code{quantile_norm}:
#' For \code{quantile_norm}, use \code{\link{quantileNorm}}.
#' @export
setMethod(
    "quantile_norm",
    signature(object = "ANY"),
    definition = function(
        object,
        quantiles = 50,
        ref_dataset = NULL,
        min_cells = 20,
        knn_k = 20,
        dims.use = NULL,
        do.center = FALSE,
        max_sample = 1000,
        eps = 0.9,
        refine.knn = TRUE,
        rand.seed = 1,
        verbose = TRUE
    ) {
        lifecycle::deprecate_warn("1.2.0", "quantile_norm()", "quantileNorm()")
        quantileNorm(
            object = object,
            quantiles = quantiles,
            reference = ref_dataset,
            minCells = min_cells,
            nNeighbors = knn_k,
            useDims = dims.use,
            center = do.center,
            maxSample = max_sample,
            eps = eps,
            refineKNN = refine.knn,
            seed = rand.seed,
            verbose = verbose
        )
    }
)
