#' Quantile align (normalize) factor loadings
#'
#' This process builds a shared factor neighborhood graph to jointly cluster cells, then quantile
#' normalizes corresponding clusters.
#'
#' The first step, building the shared factor neighborhood graph, is performed in SNF(), and
#' produces a graph representation where edge weights between cells (across all datasets)
#' correspond to their similarity in the shared factor neighborhood space. An important parameter
#' here is knn_k, the number of neighbors used to build the shared factor space.
#'
#' Next we perform quantile alignment for each dataset, factor, and cluster (by
#' stretching/compressing datasets' quantiles to better match those of the reference dataset). These
#' aligned factor loadings are combined into a single matrix and returned as H.norm.
#'
#' @param object \code{liger} object. Should run optimizeALS before calling.
#' @param knn_k Number of nearest neighbors for within-dataset knn graph (default 20).
#' @param ref_dataset Name of dataset to use as a "reference" for normalization. By default,
#'   the dataset with the largest number of cells is used.
#' @param min_cells Minimum number of cells to consider a cluster shared across datasets (default 20)
#' @param quantiles Number of quantiles to use for quantile normalization (default 50).
#' @param eps  The error bound of the nearest neighbor search. (default 0.9) Lower values give more
#' accurate nearest neighbor graphs but take much longer to computer.
#' @param dims.use Indices of factors to use for shared nearest factor determination (default
#'   1:ncol(H[[1]])).
#' @param do.center Centers the data when scaling factors (useful for less sparse modalities like
#'   methylation data). (default FALSE)
#' @param max_sample Maximum number of cells used for quantile normalization of each cluster
#' and factor. (default 1000)
#' @param refine.knn whether to increase robustness of cluster assignments using KNN graph.(default TRUE)
#' @param rand.seed Random seed to allow reproducible results (default 1)
#' @param ... Arguments passed to other methods
#'
#' @return \code{liger} object with 'H.norm' and 'clusters' slot set.
#'
#' @importFrom stats approxfun
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' # do basic quantile alignment
#' ligerex <- quantile_norm(ligerex)
#' # higher resolution for more clusters (note that SNF is conserved)
#' ligerex <- quantile_norm(ligerex, resolution = 1.2)
#' # change knn_k for more fine-grained local clustering
#' ligerex <- quantile_norm(ligerex, knn_k = 15, resolution = 1.2)
#' }

quantile_normOld <- function(object,
                          ...) {
    UseMethod(generic = 'quantile_normOld', object = object)
}

#' @rdname quantile_normOld
#' @export
#' @method quantile_normOld list
#'
quantile_normOld.list <- function(object,
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
                               ...) {
    set.seed(rand.seed)
    if (!all(sapply(X = object, FUN = is.matrix))) {
        stop("All values in 'object' must be a matrix")
    }
    if (is.null(x = names(x = object))) {
        stop("'object' must be a named list of matrices")
    }
    if (is.character(x = ref_dataset) &&
        !ref_dataset %in% names(x = object)) {
        stop("Cannot find reference dataset")
    } else if (!inherits(x = ref_dataset, what = c('character', 'numeric'))) {
        stop(
            "'ref_dataset' must be a character or integer specifying which dataset is the reference"
        )
    }
    labels <- list()
    if (is.null(dims.use)) {
        use_these_factors <- 1:ncol(object[[1]])
    } else {
        use_these_factors <- dims.use
    }
    # fast max factor assignment with Rcpp code
    labels <-
        lapply(object,
               max_factor,
               dims_use = use_these_factors,
               center_cols = do.center)
    clusters <- as.factor(unlist(lapply(labels, as.character)))
    names(clusters) <- unlist(lapply(object, rownames))

    # increase robustness of cluster assignments using knn graph
    if (refine.knn) {
        clusters <-
            refine_clusts_knn(object, clusters, k = knn_k, eps = eps)
    }
    cluster_assignments <- clusters
    clusters <- lapply(object, function(x) {
        clusters[rownames(x)]
    })
    names(clusters) <- names(object)
    dims <- ncol(object[[ref_dataset]])

    dataset <- unlist(lapply(1:length(object), function(i) {
        rep(names(object)[i], nrow(object[[i]]))
    }))
    Hs <- object
    num_clusters <- dims
    for (k in 1:length(object)) {
        for (j in 1:num_clusters) {
            cells2 <- which(clusters[[k]] == j)
            cells1 <- which(clusters[[ref_dataset]] == j)
            for (i in 1:dims) {
                num_cells2 <- length(cells2)
                num_cells1 <- length(cells1)
                if (num_cells1 < min_cells |
                    num_cells2 < min_cells) {
                    next
                }
                if (num_cells2 == 1) {
                    Hs[[k]][cells2, i] <- mean(Hs[[ref_dataset]][cells1, i])
                    next
                }
                q2 <- quantile(sample(Hs[[k]][cells2, i],
                                      min(num_cells2, max_sample)),
                               seq(0, 1, by = 1 / quantiles))
                q1 <- quantile(sample(Hs[[ref_dataset]][cells1, i],
                                      min(num_cells1, max_sample)),
                               seq(0, 1, by = 1 / quantiles))
                if (sum(q1) == 0 |
                    sum(q2) == 0 | length(unique(q1)) <
                    2 | length(unique(q2)) < 2) {
                    new_vals <- rep(0, num_cells2)
                }
                else {
                    warp_func <-
                        withCallingHandlers(
                            stats::approxfun(q2, q1, rule = 2),
                            warning = function(w) {
                                invokeRestart("muffleWarning")
                            }
                        )
                    new_vals <- warp_func(Hs[[k]][cells2, i])
                }
                Hs[[k]][cells2, i] <- new_vals
            }
        }
    }
    out <- list('H.norm' = Reduce(rbind, Hs),
                'clusters' = cluster_assignments)
    return(out)
}

#' @rdname quantile_normOld
#' @export
#' @method quantile_normOld liger
#'
quantile_normOld.liger <- function(object,
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
                                ...) {
    if (is.null(x = ref_dataset)) {
        ns <- sapply(X = object@H, FUN = nrow)
        ref_dataset <- names(x = object@H)[which.max(x = ns)]
    }
    out <- quantile_norm(
        object = object@H,
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
    for (i in names(x = out)) {
        slot(object = object, name = i) <- out[[i]]
    }
    return(object)
}
