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

#' Quantile Normalize
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
        print(length(out$clusters))
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
        labels <- list()
        if (is.null(dims.use)) {
            dims.use <- seq(nrow(object[[1]]))
        }
        # Transposing all H to cell x k
        Hs <- lapply(object, t)
        # fast max factor assignment with Rcpp code
        clusters <- lapply(Hs, max_factor,
                         dims_use = dims.use, center_cols = do.center)
        cluster_assignments <- factor(unlist(clusters))
        names(cluster_assignments) <- unlist(lapply(Hs, rownames))

        # increase robustness of cluster assignments using knn graph
        if (isTRUE(refine.knn)) {
            for (i in seq_along(Hs)) {
                clusts_H <- cluster_assignments[rownames(Hs[[i]])]
                H_knn <- RANN::nn2(Hs[[i]], eps = eps, k = knn_k,
                                   searchtype = "standard")
                # Rcpp method cluster_vote
                new_clusts <- cluster_vote(H_knn$nn.idx, clusts_H)
                cluster_assignments[rownames(Hs[[i]])] <- new_clusts
            }
        }
        dims <- ncol(Hs[[ref_dataset]])
        num_clusters <- dims
        for (k in seq_along(Hs)) {
            for (j in seq(num_clusters)) {
                cells2 <- which(clusters[[k]] == j)
                cells1 <- which(clusters[[ref_dataset]] == j)
                for (i in seq(dims)) {
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
                        sum(q2) == 0 |
                        length(unique(q1)) < 2 |
                        length(unique(q2)) < 2) {
                        new_vals <- rep(0, num_cells2)
                    } else {
                        warp_func <- withCallingHandlers(
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
    })

