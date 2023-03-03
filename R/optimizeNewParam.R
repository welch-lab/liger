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
#' \code{TRUE}.
#' @return \code{object} with \code{W} slot updated with the new \eqn{W}
#' matrix, and the \code{H} and \code{V} slots of each
#' \linkS4class{ligerDataset} object in the \code{datasets} slot updated with
#' the new dataset specific \eqn{H} and \eqn{V} matrix, respectively.
#' @export
#' @seealso \code{\link{optimizeALS}}
optimizeNewK <- function(
        object,
        k.new,
        lambda = NULL,
        thresh = 1e-4,
        max.iters = 100,
        rand.seed = 1,
        verbose = TRUE
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
        max.iters = max.iters,
        H.init = H,
        W.init = W,
        V.init = V,
        rand.seed = rand.seed,
        verbose = verbose
    )
    return(object)
}

#' Perform factorization for new data
#'
#' Uses an efficient strategy for updating that takes advantage of the information in the existing
#' factorization. Assumes that selected genes (var.genes) are represented in the new datasets.
#'
#' @param object \code{liger} object. Should call optimizeALS before calling.
#' @param new.data List of rawData matrices (one or more). Each list entry should be named.
#' @param which.datasets List of datasets to append new.data to if add.to.existing is true.
#'   Otherwise, the most similar existing datasets for each entry in new.data.
#' @param add.to.existing Add the new data to existing datasets or treat as totally new datasets
#'   (calculate new Vs?) (default TRUE)
#' @param lambda Regularization parameter. By default, this will use the lambda last used with
#'   optimizeALS.
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#'   (default 1e-4).
#' @param max.iters Maximum number of block coordinate descent iterations to perform (default 100).
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return \code{liger} object with H, W, and V slots reset. Raw.data, normData, and scaleData will
#'   also be updated to include the new data.
#'
#' @export
#' @examples
#' \dontrun{
#' # Given preprocessed liger object: ligerex (contains two datasets Y and Z)
#' # get factorization using three restarts and 20 factors
#' ligerex <- optimizeALS(ligerex, k = 20, lambda = 5, nrep = 3)
#' # acquire new data (Y_new, Z_new) from the same cell type, let's add it to existing datasets
#' new_data <- list(Y_set = Y_new, Z_set = Z_new)
#' ligerex2 <- optimizeNewData(ligerex, new.data = new_data, which.datasets = list('y_set', 'z_set'))
#' # acquire new data from different cell type (X), we'll just add another dataset
#' # it's probably most similar to y_set
#' ligerex <- optimizeNewData(ligerex, new.data = list(x_set = X), which.datasets = list('y_set'),
#'                            add.to.existing = FALSE)
#' }
optimizeNewData <- function(
        object,
        new.data,
        which.datasets,
        add.to.existing = TRUE,
        lambda = NULL,
        thresh = 1e-4,
        max.iters = 100,
        verbose = TRUE
) {
    .checkValidFactorResult(object)
    object <- recordCommand(object)
    if (is.null(lambda)) lambda <- object@uns$factorization$lambda
    sqrtLambda <- sqrt(lambda)
    nGenes <- length(varFeatures(object))
    k <- object@uns$factorization$k
    # W: g x k
    W <- getMatrix(object, "W")
    # V: g x k
    V <- getMatrix(object, "V")
    if (isTRUE(verbose)) .log("Initializing with new data...")
    if (isTRUE(add.to.existing)) {
        # TODO Establish dataset merging/extending functionality first
        for (i in 1:length(new.data)) {
            if (verbose) {
                message(dim(object@rawData[[which.datasets[[i]]]]))
            }
            object@rawData[[which.datasets[[i]]]] <-
                cbind(object@rawData[[which.datasets[[i]]]],
                      new.data[[i]])
            if (verbose) {
                message(dim(object@rawData[[which.datasets[[i]]]]))
            }
        }
        object <- normalize(object)
        object <- scaleNotCenter(object)
        H_new <- lapply(1:length(new.data), function(i) {
            t(solveNNLS(rbind(
                t(W + V[[which.datasets[[i]]]]),
                sqrtLambda * t(V[[which.datasets[[i]]]])
            ),
            rbind(
                t(object@scaleData[[which.datasets[[i]]]][colnames(new.data[[i]]), ]),
                matrix(0, nGenes, ncol(new.data[[i]]))
            )))
        })
        for (i in 1:length(new.data)) {
            object@H[[which.datasets[[i]]]] <-
                rbind(object@H[[which.datasets[[i]]]], H_new[[i]])
        }
    } else {
        new.names <- names(new.data)
        for (i in seq_along(new.names)) {
            ld <- createLigerDataset(rawData = new.data[[i]],
                                     V = dataset(object, which.datasets[i])@V)
            dataset(object, new.names[i]) <- ld
        }
        object <- normalize(object)
        object <- scaleNotCenter(object)
        nCells <- lapply(datasets(object), ncol)
        # scaleData: g x c
        E <- lapply(datasets(object), function(ld) scaleData(ld))
        # H: k x c
        H_new <- lapply(new.names, function(n) {
            solveNNLS(rbind(W + V[[n]], sqrtLambda*V[[n]]),
                      rbind(E[[n]], matrix(0, nGenes, nCells[[n]])))
        })
        for (n in new.names) {
            ld <- dataset(object, n)
            ld@H <- H_new[[n]]
            datasets(object, check = FALSE)[[n]] <- ld
        }
    }
    # H: k x c
    H <- lapply(datasets(object), function(ld) ld@H)
    object <- optimizeALS(
        object,
        k,
        lambda,
        thresh,
        max.iters,
        H.init = H,
        W.init = W,
        V.init = V,
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
#' \code{TRUE}.
#' @return Input \code{object} with optimized factorization values updated.
#' including the \code{W} matrix in \linkS4class{liger} object, and \code{W} and
#' \code{V} matrices in each \linkS4class{ligerDataset} object in the
#' \code{datasets} slot.
#' @export
optimizeNewLambda <- function(
        object,
        new.lambda,
        thresh = 1e-4,
        max.iters = 100,
        rand.seed = 1,
        verbose = TRUE
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
        rand.seed = rand.seed,
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
#' @return Subset \code{object} with factorization matrices reset, including
#' the \code{W} matrix in \linkS4class{liger} object, and \code{W} and \code{V}
#' matrices in each \linkS4class{ligerDataset} object in the \code{datasets}
#' slot. \code{scaleData} in the \linkS4class{ligerDataset} objects of
#' datasets specified by \code{datasets.scale} will also be updated to reflect
#' the subset.
#' @export
optimizeSubset <- function(
        object,
        cellIdx,
        lambda = NULL,
        thresh = 1e-4,
        max.iters = 100,
        datasets.scale = NULL,
        ...,
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
        max.iters = max.iters,
        H.init = getMatrix(object, "H"),
        W.init = getMatrix(object, "W"),
        V.init = getMatrix(object, "V"),
        ...
    )
    return(object)
}
