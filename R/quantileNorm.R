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
        rand.seed = 1
    ) standardGeneric("quantile_norm"))

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
#' @param object A \linkS4class{liger} object with \code{\link{optimizeALS}} run
#' in advance. Or a list object of \eqn{H} matrices, all in \eqn{k * nCell}
#' dimensionality, the \code{H} entry of the output of \code{\link{optimizeALS}}
#' list method is ready for use (i.e. use \code{out$H}).
#' @param knn_k Number of nearest neighbors for within-dataset knn graph.
#' Default \code{20}.
#' @param ref_dataset Name or numeric index of the dataset to use as a
#' "reference" for normalization. Default \code{NULL} use the dataset with the
#' largest number of cells.
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
#' @return For method on \linkS4class{liger} object, the object will be returned
#' with the \code{H.norm} slot and \code{"H.norm_cluster"} variable in
#' \code{cell.meta} slot updated. For method on list, a result list with entries
#' of \code{H.norm} and \code{clusters} will be returned.
#' @rdname quantile_norm
#' @export
setMethod(
    "quantile_norm",
    signature(object = "liger"),
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
        rand.seed = 1
    ) {
        .checkValidFactorResult(object)
        if (is.null(ref_dataset)) {
            nCells <- sapply(datasets(object), ncol)
            ref_dataset <- names(object)[which.max(nCells)]
        }
        out <- quantile_norm(
            object = getMatrix(object, "H"),
            quantiles = quantiles,
            ref_dataset = ref_dataset,
            min_cells = min_cells,
            knn_k = knn_k,
            dims.use = dims.use,
            do.center = do.center,
            max_sample = max_sample,
            eps = eps,
            refine.knn = refine.knn,
            rand.seed = rand.seed
        )
        object@H.norm <- out$H.norm
        cell.meta(object)$H.norm_cluster <- out$clusters
        return(object)
    })

#' @rdname quantile_norm
#' @export
setMethod(
    "quantile_norm",
    signature(object = "list"),
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
        rand.seed = 1
    ) {
        set.seed(rand.seed)
        if (!is.list(object))
            stop("`object` should be a list of matrices")
        if (is.null(names(object))) {
            stop("`object` should be a named list of matrices")
        }
        if (!all(sapply(object, is.matrix)))
            stop("All values in `object` must be a matrix")
        if (!all(sapply(object, nrow) == nrow(object[[1]])))
            stop("All matrices must have the same number of rows (k factors).")
        if (is.character(ref_dataset) &&
            !ref_dataset %in% names(object)) {
            stop("Specified `ref_dataset` is not found.")
        } else if (is.numeric(ref_dataset) &&
            ref_dataset > length(object)) {
            stop("Specified `ref_dataset` is out of bound")
        } else if (!inherits(ref_dataset, c('character', 'numeric'))) {
            stop("`ref_dataset` should be a character or integer ",
                "specifying which dataset is the reference")
        }
        if (is.null(dims.use)) {
            dims.use <- seq(nrow(object[[1]]))
        }
        # Transposing all H to cell x k
        Hs <- lapply(object, t)
        # fast max factor assignment with Rcpp code
        clusters <- lapply(Hs, max_factor,
                         dims_use = dims.use, center_cols = do.center)
        # TODO: Dumb to mess up factor with characters of numbers
        # Change to numeric in the future. Need to reproduce result for now.
        # clusterAssign <- factor(unlist(clusters))
        clusterAssign <- as.factor(unlist(lapply(clusters, as.character)))
        names(clusterAssign) <- unlist(lapply(Hs, rownames))
        # increase robustness of cluster assignments using knn graph
        if (isTRUE(refine.knn)) {
            for (i in seq_along(Hs)) {
                clustsH <- clusterAssign[rownames(Hs[[i]])]
                H_knn <- RANN::nn2(Hs[[i]], eps = eps, k = knn_k,
                                   searchtype = "standard")
                #print(clustsH[1:10])
                # Rcpp method cluster_vote
                newClusts <- cluster_vote(H_knn$nn.idx, clustsH)
                clusterAssign[rownames(Hs[[i]])] <- newClusts
            }
        }
        # TODO: This line need to be removed
        clusters <- lapply(Hs, function(H) clusterAssign[rownames(H)])
        dims <- ncol(Hs[[ref_dataset]])
        nClusters <- dims
        for (d in seq_along(Hs)) {
            for (c in seq(nClusters)) {
                # cells 2: cells in d'th dataset belonging to cluster c
                cellIdx2 <- clusters[[d]] == c
                # cells 1: cells in ref dataset belong to cluster c
                cellIdx1 <- clusters[[ref_dataset]] == c
                nCells2 <- sum(cellIdx2)
                nCells1 <- sum(cellIdx1)
                if (nCells1 < min_cells || nCells2 < min_cells) next
                for (k in seq(dims)) {
                    if (nCells2 == 1) {
                        Hs[[d]][cellIdx2, k] <-
                            mean(Hs[[ref_dataset]][cellIdx1, k])
                        next
                    }
                    q2 <- quantile(sample(Hs[[d]][cellIdx2, k],
                                          min(nCells2, max_sample)),
                                   seq(0, 1, by = 1 / quantiles))
                    q1 <- quantile(sample(Hs[[ref_dataset]][cellIdx1, k],
                                          min(nCells1, max_sample)),
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

