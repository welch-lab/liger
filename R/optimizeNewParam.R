#' Perform factorization for new value of k
#' @description This uses an efficient strategy for updating that takes
#' advantage of the information in the existing factorization. It is most
#' recommended for values of \code{k.new} smaller than current value (\code{k},
#' which is set when running \code{\link{optimizeALS}}), where it
#' is more likely to speed up the factorization.
#' @param object \linkS4class{liger} object. Should call
#' \code{\link{optimizeALS}} in advance.
#' @param k.new Number of factors of factorization.
#' @param lambda Regularization parameter. By default, this will use the lambda
#' last used with \code{\link{optimizeALS}}.
#' @param thresh Convergence threshold. Convergence occurs when
#' \eqn{|obj_0-obj|/(mean(obj0,obj)) < thresh}. Default \code{1e-4}.
#' @param max.iters Maximum number of block coordinate descent iterations to
#' perform. Default \code{100}.
#' @param rand.seed Random seed to set. Only relevant if \code{k.new} is greater
#' than \code{k}. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return \code{object} with \code{W} slot updated with the new \eqn{W}
#' matrix, and the \code{H} and \code{V} slots of each
#' \linkS4class{ligerDataset} object in the \code{datasets} slot updated with
#' the new dataset specific \eqn{H} and \eqn{V} matrix, respectively.
#' @export
#' @seealso \code{\link{optimizeALS}}
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' # Only running a few iterations for fast examples
#' pbmc <- optimizeALS(pbmc, k = 20, maxIter = 2)
#' pbmc <- optimizeNewK(pbmc, k.new = 15, max.iters = 2)
optimizeNewK <- function(
        object,
        k.new,
        lambda = NULL,
        thresh = 1e-4,
        max.iters = 100,
        rand.seed = 1,
        verbose = getOption("ligerVerbose")
) {
    .checkValidFactorResult(object)
    if (is.null(lambda)) lambda <- object@uns$factorization$lambda
    object <- recordCommand(object)
    k <- object@uns$factorization$k
    if (k.new == k) {
        return(object)
    }
    # g x c
    E <- getMatrix(object, "scaleData", returnList = TRUE)
    # g x k
    W <- getMatrix(object, "W")
    # g x k
    V <- getMatrix(object, "V", returnList = TRUE)
    # k x c
    H <- getMatrix(object, "H", returnList = TRUE)
    nGenes <- length(varFeatures(object))
    nDatasets <- length(object)
    nCells <- sapply(datasets(object), ncol)
    if (isTRUE(verbose)) .log("Initializing with new k...")
    if (k.new > k) {
        set.seed(rand.seed)
        sqrtLambda <- sqrt(lambda)
        # Initialize W_new g x k_diff
        # TODO: Directly initialize with nGene x k.new-k instead of transposing
        # Doing it now because need to reproduce old result
        W_new <- t(matrix(stats::runif(nGenes*(k.new - k), 0, 2),
                          k.new - k, nGenes))
        # Initialize V_new g x k_diff
        V_new <- lapply(seq(nDatasets), function(i)
            t(matrix(stats::runif(nGenes*(k.new - k), 0, 2),
                     k.new - k, nGenes)))
        # H_new k_diff x c
        H_new <- lapply(seq(nDatasets), function(i) {
            solveNNLS(rbind(W_new + V_new[[i]], sqrtLambda * V_new[[i]]),
                      rbind(E[[i]] - (W + V[[i]]) %*% H[[i]],
                            matrix(0, nGenes, nCells[i])))
        })
        V_new <- lapply(seq(nDatasets), function(i) {
            t(solveNNLS(
                t(cbind(H_new[[i]], sqrtLambda * H_new[[i]])),
                t(cbind(
                    E[[i]] - (W + V[[i]]) %*% H[[i]] - W_new %*% H_new[[i]],
                    matrix(0, nGenes, nCells[[i]])))
            ))
        })
        W_new <- t(solveNNLS(
            t(Reduce(cbind, H_new)),
            t(Reduce(cbind,  lapply(seq(nDatasets), function(i)
                E[[i]] - (W + V[[i]]) %*% H[[i]] - V_new[[i]] %*% H_new[[i]])))
        ))
        # H k.new x c
        H <- lapply(seq(nDatasets), function(i) rbind(H[[i]], H_new[[i]]))
        # V&W g x k.new
        V <- lapply(seq(nDatasets), function(i) cbind(V[[i]], V_new[[i]]))
        W <- cbind(W, W_new)
    } else {
        deltas <- rep(0, k)
        for (i in seq(nDatasets))
            deltas <- deltas + sapply(seq(k), function(ki)
                norm((W[, ki] + V[[i]][, ki]) %*% t(H[[i]][ki,]), "F")
            )
        k.use <- order(deltas, decreasing = TRUE)[seq(k.new)]
        W <- W[, k.use]
        V <- lapply(V, function(x) x[, k.use])
        H <- lapply(H, function(x) x[k.use, ])
    }
    object <- optimizeALS(
        object,
        k.new,
        lambda = lambda,
        thresh = thresh,
        maxIter = max.iters,
        H.init = H,
        W.init = W,
        V.init = V,
        seed = rand.seed,
        verbose = verbose
    )
    return(object)
}

#' Perform factorization for new data
#'
#' @description Uses an efficient strategy for updating that takes advantage of
#' the information in the existing factorization. Assumes that variable featuers
#' are represented in the new datasets.
#' @param object \linkS4class{liger} object. Should call
#' \code{\link{optimizeALS}} in advance.
#' @param new.data Named list of raw count matrices (one or more).
#' @param add.to.existing Logical, whether to add the new data to existing
#' datasets or treat as totally new datasets (i.e. calculate new \eqn{V}
#' matrices). Default \code{TRUE}.
#' @param which.datasets Names of datasets to append new data to if
#' \code{add.to.existing = TRUE}, or the names of datasets to inherit \eqn{V}
#' matrices from and initialize the optimization when \code{add.to.existing =
#' FALSE}. Should match the length and order of \code{new.data}.
#' @param lambda Regularization parameter. Default \code{NULL} uses the lambda
#' last used for factorization, stored at
#' \code{object@uns$factorization$lambda}.
#' @param thresh Convergence threshold. Convergence occurs when
#' \eqn{|obj_0-obj|/(mean(obj0,obj)) < thresh}. Default \code{1e-4}.
#' @param max.iters Maximum number of block coordinate descent iterations to
#' perform. Default \code{100}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return \code{object} with \code{W} slot updated with the new \eqn{W}
#' matrix, and the \code{H} and \code{V} slots of each
#' \linkS4class{ligerDataset} object in the \code{datasets} slot updated with
#' the new dataset specific \eqn{H} and \eqn{V} matrix, respectively.
#' @export
#' @seealso \code{\link{optimizeALS}}
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' # Only running a few iterations for fast examples
#' pbmc <- optimizeALS(pbmc, k = 20, maxIter = 2)
#' # Create fake new data by increasing all non-zero count in "ctrl" by 1,
#' # and make unique cell identiciers
#' ctrl2 <- rawData(dataset(pbmc, "ctrl"))
#' ctrl2@x <- ctrl2@x + 1
#' colnames(ctrl2) <- paste0(colnames(ctrl2), 2)
#' pbmcNew <- optimizeNewData(pbmc, new.data = list(ctrl2 = ctrl2),
#'                            which.datasets = "ctrl", max.iters = 2)
#' pbmcNew
optimizeNewData <- function(
        object,
        new.data,
        which.datasets,
        add.to.existing = TRUE,
        lambda = NULL,
        thresh = 1e-4,
        max.iters = 100,
        verbose = getOption("ligerVerbose")
) {
    .checkValidFactorResult(object)
    if (length(which.datasets) != length(new.data)) {
        stop("Length and order of `which.datasets` should match with
             `new.data`.")
    }
    which.datasets <- .checkUseDatasets(object, useDatasets = which.datasets)
    object <- recordCommand(object)
    if (is.null(lambda)) lambda <- object@uns$factorization$lambda
    sqrtLambda <- sqrt(lambda)
    nGenes <- length(varFeatures(object))
    k <- object@uns$factorization$k
    # W: g x k
    W <- getMatrix(object, "W")
    if (isTRUE(verbose)) .log("Initializing with new data...")
    if (isTRUE(add.to.existing)) {
        H.orig <- getMatrix(object, "H")
        # TODO Establish dataset merging/extending functionality first
        for (i in seq_along(which.datasets)) {
            rawOld <- rawData(dataset(object, which.datasets[i]))
            rawNew <- mergeSparseAll(list(rawOld, new.data[[i]]))
            ld <- createLigerDataset(rawData = rawNew,
                                     V = getMatrix(object, "V",
                                                   dataset = which.datasets[i],
                                                   returnList = FALSE))
            dataset(object, which.datasets[i]) <- ld
        }
        object <- normalize(object, useDatasets = which.datasets)
        object <- scaleNotCenter(object, useDatasets = which.datasets)
        # scaleData: g x c
        E <- getMatrix(object, "scaleData")
        # V: g x k
        V <- getMatrix(object, "V")
        # H: k x c
        H_new <- lapply(1:length(new.data), function(i) {
            solveNNLS(
                rbind(
                    W + V[[which.datasets[i]]],
                    sqrtLambda * V[[which.datasets[i]]]
                ),
                rbind(
                    E[[which.datasets[i]]][,colnames(new.data[[i]])],
                    matrix(0, nGenes, ncol(new.data[[i]]))
                )
            )
        })
        names(H_new) <- which.datasets
        for (n in which.datasets) {
            ld <- dataset(object, n)
            ld@H <- cbind(H.orig[[n]], H_new[[n]])
            datasets(object, check = FALSE)[[n]] <- ld
        }
    } else {
        new.names <- names(new.data)
        if (is.null(new.names)) {
            stop("`new.names` has to be a named list when ",
                 "`add.to.existing` = FALSE.")
        }
        for (i in seq_along(new.names)) {
            ld <- as.ligerDataset(new.data[[i]])
            ld@V <- getMatrix(object, "V", dataset = which.datasets[i])
            dataset(object, new.names[i]) <- ld
        }
        object <- normalize(object, useDatasets = new.names)
        object <- scaleNotCenter(object, useDatasets = new.names)
        nCells <- lapply(datasets(object), ncol)
        # scaleData: g x c
        E <- getMatrix(object, "scaleData")
        # V: g x k
        V <- getMatrix(object, "V")
        # H: k x c
        H_new <- lapply(new.names, function(n) {
            solveNNLS(rbind(W + V[[n]], sqrtLambda*V[[n]]),
                      rbind(E[[n]], matrix(0, nGenes, nCells[[n]])))
        })
        names(H_new) <- new.names
        for (n in new.names) {
            ld <- dataset(object, n)
            ld@H <- H_new[[n]]
            datasets(object, check = FALSE)[[n]] <- ld
        }
    }
    object <- optimizeALS(
        object,
        k,
        lambda,
        thresh,
        max.iters,
        H.init = getMatrix(object, "H"),
        W.init = getMatrix(object, "W"),
        V.init = getMatrix(object, "V"),
        verbose = verbose
    )
    return(object)
}

#' Perform factorization for new lambda value
#' @description Uses an efficient strategy for updating that takes advantage of
#' the information in the existing factorization; uses previous k. Recommended
#' mainly when re-optimizing for higher lambda and when new lambda value is
#' significantly different; otherwise may not return optimal results.
#' @param object \linkS4class{liger} object. Should call
#' \code{\link{optimizeALS}} in advance.
#' @param new.lambda Numeric regularization parameter. Larger values penalize
#' dataset-specific effects more strongly.
#' @param thresh Convergence threshold. Convergence occurs when
#' \eqn{|obj_0-obj|/(mean(obj_0,obj)) < thresh}. Default \code{1e-4}.
#' @param max.iters Maximum number of block coordinate descent iterations to
#' perform. Default \code{100}.
#' @param rand.seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return Input \code{object} with optimized factorization values updated.
#' including the \code{W} matrix in \linkS4class{liger} object, and \code{W} and
#' \code{V} matrices in each \linkS4class{ligerDataset} object in the
#' \code{datasets} slot.
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' # Only running a few iterations for fast examples
#' pbmc <- optimizeALS(pbmc, k = 20, maxIter = 2)
#' pbmc <- optimizeNewLambda(pbmc, new = 5.5, max.iters = 2)
optimizeNewLambda <- function(
        object,
        new.lambda,
        thresh = 1e-4,
        max.iters = 100,
        rand.seed = 1,
        verbose = getOption("ligerVerbose")
) {
    .checkValidFactorResult(object, names(object))
    if (new.lambda < object@uns$factorization$lambda && isTRUE(verbose))
        .log("New lambda less than current lambda; new factorization may not ",
             "be optimal. Re-optimization with optimizeAlS recommended ",
             "instead.")
    object <- optimizeALS(
        object,
        k = object@uns$factorization$k,
        lambda = new.lambda,
        thresh = thresh,
        maxIter = max.iters,
        H.init = getMatrix(object, "H"),
        W.init = getMatrix(object, "W"),
        seed = rand.seed,
        verbose = verbose
    )
    return(object)
}

#' Perform factorization for subset of data
#' @description Uses an efficient strategy for updating that takes advantage of
#' the information in the existing factorization.
#' @param object \linkS4class{liger} object. Should call
#' \code{\link{optimizeALS}} in advance.
#' @param cellIdx Valid index vector that applies to the whole object. See
#' \code{\link{subsetLiger}} for requirement.
#' @param lambda Numeric regularization parameter. Default \code{5}.
#' @param thresh Convergence threshold. Convergence occurs when
#' \eqn{|obj_0-obj|/(mean(obj_0,obj)) < thresh}. Default \code{1e-4}.
#' @param max.iters Maximum number of block coordinate descent iterations to
#' perform. Default \code{100}.
#' @param datasets.scale Names of datasets to rescale after subsetting.
#' Default \code{NULL} does not rescale.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param cell.subset,cluster.subset \bold{Deprecated}. Please use
#' \code{cellIdx} to explicitly specify.
#' @return Subset \code{object} with factorization matrices reset, including
#' the \code{W} matrix in \linkS4class{liger} object, and \code{W} and \code{V}
#' matrices in each \linkS4class{ligerDataset} object in the \code{datasets}
#' slot. \code{scaleData} in the \linkS4class{ligerDataset} objects of
#' datasets specified by \code{datasets.scale} will also be updated to reflect
#' the subset.
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' # Only running a few iterations for fast examples
#' pbmc <- optimizeALS(pbmc, k = 20, maxIter = 2)
#' pbmc <- optimizeSubset(pbmc, cellIdx = sort(sample(ncol(pbmc), 200)),
#'                        max.iters = 2)
optimizeSubset <- function(
        object,
        cellIdx,
        lambda = NULL,
        thresh = 1e-4,
        max.iters = 100,
        datasets.scale = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        # Deprecated
        cell.subset = NULL,
        cluster.subset = NULL
) {
    .deprecateArgs(list(cell.subset = "cellIdx"), defunct = "cluster.subset")
    .checkValidFactorResult(object, useDatasets = names(object))
    if (is.null(lambda)) lambda <- object@uns$factorization$lambda
    cellIdx <- .idxCheck(object, cellIdx, orient = "cell")
    if (!is.null(datasets.scale))
        datasets.scale <- .checkUseDatasets(object, datasets.scale)
    object <- recordCommand(object)

    if (!is.null(datasets.scale))
        object <- scaleNotCenter(object, useDatasets = datasets.scale)
    object <- subsetLiger(object, cellIdx = cellIdx)
    object <- optimizeALS(
        object,
        k = object@uns$factorization$k,
        lambda = lambda,
        thresh = thresh,
        maxIter = max.iters,
        H.init = getMatrix(object, "H"),
        W.init = getMatrix(object, "W"),
        V.init = getMatrix(object, "V"),
        seed = seed,
        verbose = verbose,
    )
    return(object)
}
