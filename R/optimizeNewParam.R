#' Perform factorization for new value of k
#' @description This uses an efficient strategy for updating that takes
#' advantage of the information in the existing factorization. It is most
#' recommended for values of \code{kNew} smaller than current value (\code{k},
#' which is set when running \code{\link{runINMF}}), where it is more likely to
#' speed up the factorization.
#' @param object A \linkS4class{liger} object. Should have integrative
#' factorization performed e.g. (\code{\link{runINMF}}) in advance.
#' @param kNew Number of factors of factorization.
#' @param lambda Numeric regularization parameter. By default \code{NULL}, this
#' will use the lambda value used in the latest factorization.
#' @param nIteration Number of block coordinate descent iterations to
#' perform. Default \code{30}.
#' @param seed Random seed to allow reproducible results. Default \code{1}. Used
#' by \code{\link{runINMF}} factorization and initialization only when if
#' \code{kNew} is greater than \code{k}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param k.new,max.iters,rand.seed These arguments are now replaced by others
#' and will be removed in the future. Please see usage for replacement.
#' @param thresh \bold{Deprecated}. New implementation of iNMF does not require
#' a threshold for convergence detection. Setting a large enough
#' \code{nIteration} will bring it to convergence.
#' @return \code{object} with \code{W} slot updated with the new \eqn{W}
#' matrix, and the \code{H} and \code{V} slots of each
#' \linkS4class{ligerDataset} object in the \code{datasets} slot updated with
#' the new dataset specific \eqn{H} and \eqn{V} matrix, respectively.
#' @export
#' @seealso \code{\link{runINMF}}, \code{\link{optimizeNewLambda}},
#' \code{\link{optimizeNewData}}
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' # Only running a few iterations for fast examples
#' if (requireNamespace("RcppPlanc", quietly = TRUE)) {
#'     pbmc <- runINMF(pbmc, k = 20, nIteration = 2)
#'     pbmc <- optimizeNewK(pbmc, kNew = 25, nIteration = 2)
#' }
optimizeNewK <- function(
        object,
        kNew,
        lambda = NULL,
        nIteration = 30,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        k.new = kNew,
        max.iters = nIteration,
        rand.seed = seed,
        thresh = NULL
) {
    ## TODO: decide whether to move that initialization to C++ as well,
    ## depending on whether this function is commonly used enough
    .checkValidFactorResult(object)
    lambda <- lambda %||% object@uns$factorization$lambda
    .deprecateArgs(replace = list(k.new = "kNew",
                                  max.iters = "nIteration",
                                  rand.seed = "seed"),
                   defunct = "thresh")
    object <- recordCommand(object)
    k <- object@uns$factorization$k
    if (kNew == k) {
        return(object)
    }
    # g x c
    E <- getMatrix(object, "scaleData", returnList = TRUE)
    # E <- lapply(E, as.matrix)
    # g x k
    W <- getMatrix(object, "W")
    # g x k
    V <- getMatrix(object, "V", returnList = TRUE)
    # k x c
    H <- getMatrix(object, "H", returnList = TRUE)
    nGenes <- length(varFeatures(object))
    nDatasets <- length(object)
    if (isTRUE(verbose)) .log("Initializing with new k...")
    if (kNew > k) {
        set.seed(seed)
        # sqrtLambda <- sqrt(lambda)
        # Initialize W_new g x k_diff
        # TODO: Directly initialize with nGene x kNew-k instead of transposing
        # Doing it now because need to reproduce old result
        W_new <- t(matrix(stats::runif(nGenes*(kNew - k), 0, 2),
                          kNew - k, nGenes))
        # Initialize V_new g x k_diff
        V_new <- lapply(seq(nDatasets), function(i)
            t(matrix(stats::runif(nGenes*(kNew - k), 0, 2),
                     kNew - k, nGenes)))
        # H_new k_diff x c
        H_new <- lapply(seq(nDatasets), function(i) {
            # High level idea is to mathematically directly derive the result
            # of CtC and CtB without creating the giant rbind'ed 'g x c' or
            # even '2g x c' matrices
            # C <- rbind(W_new + V_new[[i]], sqrtLambda * V_new[[i]])
            # B <- as.matrix(rbind(E[[i]] - (W + V[[i]]) %*% H[[i]],
            #                      matrix(0, nGenes, nCells[i])))
            # RcppPlanc::bppnnls(C, B)
            # RcppPlanc::bppnnls_prod(t(C) %*% C, t(C) %*% B)
            L_new <- W_new + V_new[[i]]
            L <- W + V[[i]]
            CtC <- t(L_new) %*% (L_new) +
                lambda * t(V_new[[i]]) %*% V_new[[i]]
            CtB <- -1 * as.matrix(t(L_new) %*% L %*% H[[i]] - t(L_new) %*% E[[i]])
            RcppPlanc::bppnnls_prod(CtC, CtB)
        })
        V_new <- lapply(seq(nDatasets), function(i) {
            # Similarly as what we did with H_new
            # C <- t(cbind(H_new[[i]], sqrtLambda * H_new[[i]]))
            # B <- t(cbind(
            #     E[[i]] - (W + V[[i]]) %*% H[[i]] - W_new %*% H_new[[i]],
            #     matrix(0, nGenes, nCells[[i]])
            # ))
            # t(RcppPlanc::bppnnls(C, B))
            CtC <- H_new[[i]] %*% t(H_new[[i]]) * (1 + lambda)
            CtB <- H_new[[i]] %*% t(E[[i]]) -
                H_new[[i]] %*% t(H[[i]]) %*% t(W + V[[i]]) -
                H_new[[i]] %*% t(H_new[[i]]) %*% t(W_new)
            t(RcppPlanc::bppnnls_prod(CtC, as.matrix(CtB)))
        })
        # Similarly for solving W_new
        # C <- t(Reduce(cbind, H_new))
        # B <- t(Reduce(cbind,  lapply(seq(nDatasets), function(i)
        #     E[[i]] - (W + V[[i]]) %*% H[[i]] - V_new[[i]] %*% H_new[[i]])))
        # W_new <- t(RcppPlanc::bppnnls(C, B))
        CtC <- matrix(0, kNew - k, kNew - k)
        for (i in seq_along(E)) CtC <- CtC + H_new[[i]] %*% t(H_new[[i]])
        CtB <- matrix(0, kNew - k, nGenes)
        for (i in seq_along(E)) {
            CtB <- CtB +
                H_new[[i]] %*% t(E[[i]]) -
                H_new[[i]] %*% t(H[[i]]) %*% t(W + V[[i]]) -
                H_new[[i]] %*% t(H_new[[i]]) %*% t(V_new[[i]])
        }
        W_new <- t(RcppPlanc::bppnnls_prod(CtC, as.matrix(CtB)))
        # H kNew x c
        H <- lapply(seq(nDatasets), function(i) t(rbind(H[[i]], H_new[[i]])))
        # V&W g x kNew
        V <- lapply(seq(nDatasets), function(i) cbind(V[[i]], V_new[[i]]))
        W <- cbind(W, W_new)
    } else {
        deltas <- rep(0, k)
        for (i in seq(nDatasets))
            deltas <- deltas + sapply(seq(k), function(ki)
                norm((W[, ki] + V[[i]][, ki]) %*% t(H[[i]][ki,]), "F")
            )
        k.use <- order(deltas, decreasing = TRUE)[seq(kNew)]
        W <- W[, k.use]
        V <- lapply(V, function(x) x[, k.use])
        H <- lapply(H, function(x) t(x[k.use, ]))
    }
    object <- runINMF(
        object = object,
        k = kNew, lambda = lambda,
        nIteration = nIteration,
        nRandomStarts = 1,
        HInit = H,
        WInit = W,
        VInit = V,
        seed = seed,
        verbose = verbose
    )
    return(object)
}

#' Perform factorization for new data
#'
#' @description Uses an efficient strategy for updating that takes advantage of
#' the information in the existing factorization. Assumes that variable features
#' are presented in the new datasets. Two modes are supported (controlled by
#' \code{merge}):
#' \itemize{
#' \item{Append new data to existing datasets specified by \code{useDatasets}.
#' Here the existing \eqn{V} matrices for the target datasets will directly be
#' used as initialization, and new \eqn{H} matrices for the merged matrices will
#' be initialized accordingly.}
#' \item{Set new data as new datasets. Initial \eqn{V} matrices for them will
#' be copied from datasets specified by \code{useDatasets}, and new \eqn{H}
#' matrices will be initialized accordingly.}
#' }
#' @param object A \linkS4class{liger} object. Should have integrative
#' factorization performed e.g. (\code{\link{runINMF}}) in advance.
#' @param dataNew Named list of \bold{raw count} matrices, genes by cells.
#' @param useDatasets Selection of datasets to append new data to if
#' \code{merge = TRUE}, or the datasets to inherit \eqn{V} matrices from and
#' initialize the optimization when \code{merge = FALSE}. Should match the
#' length and order of \code{dataNew}.
#' @param merge Logical, whether to add the new data to existing
#' datasets or treat as totally new datasets (i.e. calculate new \eqn{V}
#' matrices). Default \code{TRUE}.
#' @param lambda Numeric regularization parameter. By default \code{NULL}, this
#' will use the lambda value used in the latest factorization.
#' @param nIteration Number of block coordinate descent iterations to perform.
#' Default \code{30}.
#' @param seed Random seed to allow reproducible results. Default \code{1}. Used
#' by \code{\link{runINMF}} factorization.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param new.data,which.datasets,add.to.existing,max.iters These arguments are
#' now replaced by others and will be removed in the future. Please see usage
#' for replacement.
#' @param thresh \bold{Deprecated}. New implementation of iNMF does not require
#' a threshold for convergence detection. Setting a large enough
#' \code{nIteration} will bring it to convergence.
#'
#'
#' @return \code{object} with \code{W} slot updated with the new \eqn{W}
#' matrix, and the \code{H} and \code{V} slots of each
#' \linkS4class{ligerDataset} object in the \code{datasets} slot updated with
#' the new dataset specific \eqn{H} and \eqn{V} matrix, respectively.
#' @export
#' @seealso \code{\link{runINMF}}, \code{\link{optimizeNewK}},
#' \code{\link{optimizeNewLambda}}
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' # Only running a few iterations for fast examples
#' if (requireNamespace("RcppPlanc", quietly = TRUE)) {
#'     pbmc <- runINMF(pbmc, k = 20, nIteration = 2)
#'     # Create fake new data by increasing all non-zero count in "ctrl" by 1,
#'     # and make unique cell identifiers
#'     ctrl2 <- rawData(dataset(pbmc, "ctrl"))
#'     ctrl2@x <- ctrl2@x + 1
#'     colnames(ctrl2) <- paste0(colnames(ctrl2), 2)
#'     pbmcNew <- optimizeNewData(pbmc, dataNew = list(ctrl2 = ctrl2),
#'                                useDatasets = "ctrl", nIteration = 2)
#' }
optimizeNewData <- function(
        object,
        dataNew,
        useDatasets,
        merge = TRUE,
        lambda = NULL,
        nIteration = 30,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        new.data = dataNew,
        which.datasets = useDatasets,
        add.to.existing = merge,
        max.iters = nIteration,
        thresh = NULL
) {
    .checkValidFactorResult(object)
    if (length(which.datasets) != length(dataNew)) {
        stop("Length and order of `which.datasets` should match with
             `dataNew`.")
    }
    .deprecateArgs(list(new.data = "dataNew", which.datasets = "useDatasets",
                        max.iters = "nIteration"), "thresh")
    if (is.null(names(dataNew))) {
        stop("`dataNew` has to be a named list.")
    }
    useDatasets <- .checkUseDatasets(object, useDatasets = useDatasets)
    object <- recordCommand(object)
    lambda <- lambda %||% object@uns$factorization$lambda
    k <- object@uns$factorization$k
    # W: g x k
    W <- getMatrix(object, "W")

    if (isTRUE(merge)) {
        if (isTRUE(verbose)) {
            .log("Initializing with new data merged to existing datasets...")
        }
        H.orig <- getMatrix(object, "H")
        for (i in seq_along(useDatasets)) {
            rawOld <- rawData(object, dataset = useDatasets[i])
            rawNew <- mergeSparseAll(list(rawOld, dataNew[[i]]))
            ld <- createLigerDataset(rawData = rawNew,
                                     modal = modalOf(object)[useDatasets[i]],
                                     V = getMatrix(object, "V",
                                                   dataset = useDatasets[i],
                                                   returnList = FALSE))
            dataset(object, useDatasets[i]) <- ld
        }
        object <- normalize(object, useDatasets = useDatasets, verbose = verbose)
        object <- scaleNotCenter(object, useDatasets = useDatasets, verbose = verbose)
        # scaleData: g x c
        E <- getMatrix(object, "scaleData", dataset = useDatasets, returnList = TRUE)
        # V: g x k
        V <- getMatrix(object, "V", dataset = useDatasets, returnList = TRUE)
        # H: k x c
        H_new <- lapply(seq_along(dataNew), function(i) {
            idx <- useDatasets[i]
            # C <- rbind(W + V[[idx]], sqrtLambda * V[[idx]])
            # B <- rbind(E[[idx]][,colnames(dataNew[[i]])],
            #            matrix(0, nGenes, ncol(dataNew[[i]])))
            # RcppPlanc::bppnnls(C, B)
            CtC <- t(W + V[[idx]]) %*% (W + V[[idx]]) +
                lambda * t(V[[idx]]) %*% V[[idx]]
            CtB <- t(W + V[[idx]]) %*% E[[idx]][,colnames(dataNew[[i]])]
            cbind(H.orig[[idx]], RcppPlanc::bppnnls_prod(CtC, as.matrix(CtB)))
        })
        names(H_new) <- useDatasets
        for (n in useDatasets) {
            object@datasets[[n]]@H <- H_new[[n]]
            # ld <- dataset(object, n)
            # ld@H <- H_new[[n]]
            # datasets(object, check = FALSE)[[n]] <- ld
        }
    } else {
        if (isTRUE(verbose)) {
            .log("Initializing with new data added as new datasets...")
        }
        new.names <- names(dataNew)
        if (any(new.names %in% names(object))) {
            stop("Names of `dataNew` must be unique and different from ",
                 "exsiting datasets.")
        }
        for (i in seq_along(new.names)) {
            ld <- createLigerDataset(
                dataNew[[i]],
                V = getMatrix(object, "V", dataset = useDatasets[i])
            )
            dataset(object, new.names[i]) <- ld
        }
        object <- normalize(object, useDatasets = new.names, verbose = verbose)
        object <- scaleNotCenter(object, useDatasets = new.names, verbose = verbose)
        # scaleData: g x c
        E <- getMatrix(object, "scaleData", dataset = new.names, returnList = TRUE)
        # V: g x k
        V <- getMatrix(object, "V", dataset = new.names, returnList = TRUE)
        # H: k x c
        H_new <- lapply(new.names, function(n) {
            # C <- rbind(W + V[[n]], sqrtLambda*V[[n]])
            # B <- rbind(E[[n]], matrix(0, nGenes, nCells[[n]]))
            # RcppPlanc::bppnnls(C, B)
            CtC <- t(W + V[[n]]) %*% (W + V[[n]]) + lambda * (t(V[[n]]) %*% V[[n]])
            CtB <- as.matrix(t(W + V[[n]]) %*% E[[n]])
            RcppPlanc::bppnnls_prod(CtC, CtB)
        })
        names(H_new) <- new.names
        for (n in new.names) {
            object@datasets[[n]]@H <- H_new[[n]]
            # ld <- dataset(object, n)
            # ld@H <- H_new[[n]]
            # datasets(object, check = FALSE)[[n]] <- ld
        }
    }
    object <- runINMF(
        object,
        k = k,
        lambda = lambda,
        nIteration = nIteration,
        HInit = lapply(getMatrix(object, "H"), t),
        WInit = getMatrix(object, "W"),
        VInit = getMatrix(object, "V"),
        verbose = verbose,
        seed = seed
    )
    return(object)
}

#' Perform factorization for new lambda value
#' @description Uses an efficient strategy for updating that takes advantage of
#' the information in the existing factorization; always uses previous k.
#' Recommended mainly when re-optimizing for higher lambda and when new lambda
#' value is significantly different; otherwise may not return optimal results.
#' @param object \linkS4class{liger} object. Should have integrative
#' factorization (e.g. \code{\link{runINMF}}) performed in advance.
#' @param lambdaNew Numeric regularization parameter. Larger values penalize
#' dataset-specific effects more strongly.
#' @param nIteration Number of block coordinate descent iterations to
#' perform. Default \code{30}.
#' @param seed Random seed to allow reproducible results. Default \code{1}. Used
#' by \code{\link{runINMF}} factorization.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param new.lambda,max.iters,rand.seed These arguments are now replaced by
#' others and will be removed in the future. Please see usage for replacement.
#' @param thresh \bold{Deprecated}. New implementation of iNMF does not require
#' a threshold for convergence detection. Setting a large enough
#' \code{nIteration} will bring it to convergence.
#' @return Input \code{object} with optimized factorization values updated.
#' including the \code{W} matrix in \linkS4class{liger} object, and \code{H} and
#' \code{V} matrices in each \linkS4class{ligerDataset} object in the
#' \code{datasets} slot.
#' @export
#' @seealso \code{\link{runINMF}}, \code{\link{optimizeNewK}},
#' \code{\link{optimizeNewData}}
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' if (requireNamespace("RcppPlanc", quietly = TRUE)) {
#'     # Only running a few iterations for fast examples
#'     pbmc <- runINMF(pbmc, k = 20, nIteration = 2)
#'     pbmc <- optimizeNewLambda(pbmc, lambdaNew = 5.5, nIteration = 2)
#' }
optimizeNewLambda <- function(
        object,
        lambdaNew,
        nIteration = 30,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        new.lambda = lambdaNew,
        max.iters = nIteration,
        rand.seed = seed,
        thresh = NULL
) {
    .checkValidFactorResult(object, names(object))
    .deprecateArgs(list(new.lambda = "lambdaNew", max.iters = "nIteration",
                        rand.seed = "seed"), "thresh")
    object <- recordCommand(object)
    if (lambdaNew < object@uns$factorization$lambda && isTRUE(verbose))
        .log("New lambda less than current lambda; new factorization may not ",
             "be optimal. Re-optimization with optimizeAlS recommended ",
             "instead.")

    object <- runINMF(
        object,
        k = object@uns$factorization$k,
        lambda = lambdaNew,
        nIteration = nIteration,
        HInit = getMatrix(object, "H"),
        WInit = getMatrix(object, "W"),
        seed = seed,
        verbose = verbose
    )
    return(object)
}

#' Perform factorization for subset of data
#' @description Uses an efficient strategy for updating that takes advantage of
#' the information in the existing factorization.
#' @param object \linkS4class{liger} object. Should have integrative
#' factorization (e.g. \code{\link{runINMF}}) performed in advance.
#' @param clusterVar,useClusters Together select the clusters to subset the
#' object conveniently. \code{clusterVar} is the name of variable in
#' \code{cellMeta(object)} and \code{useClusters} should be vector of names of
#' clusters in the variable. \code{clusterVar} is by default the default
#' cluster (See \code{\link{runCluster}}, or \code{\link{defaultCluster}} at
#' "Cell metadata access"). Users can otherwise select cells explicitly with
#' \code{cellIdx} for complex conditions. \code{useClusters} overrides
#' \code{cellIdx}.
#' @param lambda Numeric regularization parameter. By default \code{NULL}, this
#' will use the lambda value used in the latest factorization.
#' @param nIteration Maximum number of block coordinate descent iterations to
#' perform. Default \code{30}.
#' @param cellIdx Valid index vector that applies to the whole object. See
#' \code{\link{subsetLiger}} for requirement. Default \code{NULL}.
#' @param scaleDatasets Names of datasets to re-scale after subsetting.
#' Default \code{NULL} does not re-scale.
#' @param seed Random seed to allow reproducible results. Default \code{1}. Used
#' by \code{\link{runINMF}} factorization.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param cell.subset,cluster.subset,max.iters,datasets.scale These arguments
#' are now replaced by others and will be removed in the future. Please see
#' usage for replacement.
#' @param thresh \bold{Deprecated}. New implementation of iNMF does not require
#' a threshold for convergence detection. Setting a large enough
#' \code{nIteration} will bring it to convergence.
#' @return Subset \code{object} with factorization matrices optimized, including
#' the \code{W} matrix in \linkS4class{liger} object, and \code{W} and \code{V}
#' matrices in each \linkS4class{ligerDataset} object in the \code{datasets}
#' slot. \code{scaleData} in the \linkS4class{ligerDataset} objects of
#' datasets specified by \code{scaleDatasets} will also be updated to reflect
#' the subset.
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' if (requireNamespace("RcppPlanc", quietly = TRUE)) {
#'     # Only running a few iterations for fast examples
#'     pbmc <- runINMF(pbmc, k = 20, nIteration = 2)
#'     pbmc <- optimizeSubset(pbmc, cellIdx = sort(sample(ncol(pbmc), 200)),
#'                            nIteration = 2)
#' }
optimizeSubset <- function(
        object,
        clusterVar = NULL,
        useClusters = NULL,
        lambda = NULL,
        nIteration = 30,
        cellIdx = NULL,
        scaleDatasets = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        # Deprecated
        cell.subset = cellIdx,
        cluster.subset = useClusters,
        max.iters = nIteration,
        datasets.scale = scaleDatasets,
        thresh = NULL
) {
    .deprecateArgs(list(cell.subset = "cellIdx", cluster.subset = "useClusters",
                        max.iters = "nIteration",
                        datasets.scale = "scaleDatasets"),
                   defunct = "thresh")
    .checkValidFactorResult(object, useDatasets = names(object))
    lambda <- lambda %||% object@uns$factorization$lambda
    clusterVar <- clusterVar %||% object@uns$defaultCluster
    if (!is.null(useClusters)) {
        if (is.null(clusterVar)) {
            stop("No `clusterVar` specified. Specify variable name with
                 `clusterVar`, or see `?optimizeSubset`.")
        }
        clusterVar <- .fetchCellMetaVar(object, clusterVar,
                                        checkCategorical = TRUE)
        if (all(!useClusters %in% levels(clusterVar))) {
            stop("`useCluster` must contain existing levels in ",
                 "`object[[clusterVar]]`.")
        }
        cellIdx <- which(clusterVar %in% useClusters)
    }
    if (length(cellIdx) == 0) {
        stop("No subset specified. Use either a combination of `clusterVar` ",
             "and `useClusters` or explicit `cellIdx`.")
    }
    cellIdx <- .idxCheck(object, cellIdx, orient = "cell")
    object <- recordCommand(object)
    object <- subsetLiger(object, cellIdx = cellIdx, verbose = verbose)
    if (!is.null(scaleDatasets))
        object <- scaleNotCenter(object, useDatasets = scaleDatasets,
                                 verbose = verbose)
    object <- runINMF(
        object,
        k = object@uns$factorization$k,
        lambda = lambda,
        nIteration = nIteration,
        HInit = lapply(getMatrix(object, "H"), t),
        WInit = getMatrix(object, "W"),
        VInit = getMatrix(object, "V"),
        seed = seed,
        verbose = verbose
    )
    return(object)
}
