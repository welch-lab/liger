#' Perform iNMF on scaled datasets
#' @description
#' Performs integrative non-negative matrix (iNMF) factorization to return
#' factorized \eqn{H}, \eqn{W}, and \eqn{V} matrices. It optimizes the iNMF
#' objective function using block coordinate descent (alternating non-negative
#' least squares), where the number of factors is set by \code{k}. TODO: include
#' objective function equation here in documentation (using deqn)
#'
#' For each dataset, this factorization produces an \eqn{H} matrix (cells by k),
#' a \eqn{V} matrix (k by genes), and a shared \eqn{W} matrix (k by genes). The
#' \eqn{H} matrices represent the cell factor loadings. \eqn{W} is held
#' consistent among all datasets, as it represents the shared components of the
#' metagenes across datasets. The \eqn{V} matrices represent the
#' dataset-specific components of the metagenes.
#' @param object A \linkS4class{liger} object or a named list of matrix object,
#' where the names represents dataset names and matrices are scaled on the same
#' set of variable features, with rows as features and columns as cells.
#' @param k Inner dimension of factorization (number of factors). Run
#' \code{\link{suggestK}} to determine appropriate value; a general rule of
#' thumb is that a higher \code{k} will be needed for datasets with more
#' sub-structure.
#' @param lambda Regularization parameter. Larger values penalize
#' dataset-specific effects more strongly (i.e. alignment should increase as
#' \code{lambda} increases). Default \code{5}.
#' @param thresh Convergence threshold. Convergence occurs when
#' \eqn{|obj_0-obj|/(mean(obj_0,obj)) < thresh}. Default \code{1e-6}.
#' @param maxIter Maximum number of block coordinate descent iterations to
#' perform. Default \code{30}.
#' @param nrep Number of restarts to perform (iNMF objective function is
#' non-convex, so taking the best objective from multiple successive
#' initialization is recommended). For easier reproducibility, this increments
#' the random seed by 1 for each consecutive restart, so future factorization
#' of the same dataset can be run with one rep if necessary. Default \code{1}.
#' @param H.init Initial values to use for \eqn{H} matrices. A list object where
#' each element is the initial \eqn{H} matrix of each dataset. Default
#' \code{NULL}.
#' @param W.init Initial values to use for \eqn{W} matrix. A matrix object.
#' Default \code{NULL}.
#' @param V.init Initial values to use for \eqn{V} matrices. A list object where
#' each element is the initial \eqn{V} matrix of each dataset. Default
#' \code{NULL}.
#' @param method NNLS subproblem solver. Choose from \code{"liger"} (default
#' original implementation), \code{"planc"} or \code{"rcppml"}.
#' @param useUnshared Logical, whether to include unshared variable features and
#' run optimizeUANLS algorithm. Defaul \code{FALSE}. Running
#' \code{\link{selectGenes}} with \code{unshared = TRUE} and then running
#' \code{\link{scaleNotCenter}} is required.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param readH5 \code{TRUE} to force reading H5 based data into memory and
#' conduct factorization. \code{"auto"} reads H5 dataset with less than 8000
#' cells. \code{FALSE} will stop users from running if H5 data presents.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param max.iters,use.unshared,rand.seed \bold{Deprecated}. See Usage section
#' for replacement.
#' @param print.obj \bold{Defunct}. Whether to print objective function values
#' after convergence when \code{verbose = TRUE}. Now always print when verbose.
#' @return \code{object} with \code{W} slot updated with the result \eqn{W}
#' matrix, and the \code{H} and \code{V} slots of each
#' \linkS4class{ligerDataset} object in the \code{datasets} slot updated with
#' the dataset specific \eqn{H} and \eqn{V} matrix, respectively.
#' @rdname runINMF_R
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' # Only running a few iterations for fast examples
#' pbmc <- runINMF(pbmc, k = 20, maxIter = 2)
setGeneric(
    "runINMF_R",
    function(
        object,
        k,
        lambda = 5.0,
        thresh = 1e-6,
        maxIter = 30,
        nrep = 1,
        H.init = NULL,
        W.init = NULL,
        V.init = NULL,
        method = c("planc", "liger", "rcppml"),
        useUnshared = FALSE,
        seed = 1,
        readH5 = "auto",
        verbose = getOption("ligerVerbose")
    ) standardGeneric("runINMF_R")
)

#' @rdname runINMF_R
#' @export
setMethod(
    "runINMF_R",
    signature(object = "liger"),
    function(
        object,
        k,
        lambda = 5.0,
        thresh = 1e-6,
        maxIter = 30,
        nrep = 1,
        H.init = NULL,
        W.init = NULL,
        V.init = NULL,
        useUnshared = FALSE,
        seed = 1,
        readH5 = "auto",
        verbose = getOption("ligerVerbose")
    ) {
        .checkObjVersion(object)
        object <- recordCommand(object)
        if (isFALSE(useUnshared)) {
            object <- removeMissing(object, orient = "cell",
                                    verbose = verbose)
            data <- lapply(datasets(object), function(ld) {
                if (is.null(scaleData(ld)))
                    stop("Scaled data not available. ",
                         "Run `scaleNotCenter(object)` first")
                if (isH5Liger(ld)) {
                    if (!isFALSE(readH5)) {
                        h5d <- scaleData(ld)
                        if (readH5 == "auto") {
                            if (h5d$dims[2] <= 8000) {
                                warning("Automatically reading H5 based ",
                                        "scaled dense matrix into memory. ",
                                        "Dim: ", h5d$dims[1], "x", h5d$dims[2],
                                        immediate. = verbose)
                                return(h5d[,])
                            } else {
                                stop("Scaled data in H5 based dataset with ",
                                     "more than 8000 cells will not be ",
                                     "automatically read into memory. Use ",
                                     "`readH5 = TRUE` to force reading, or ",
                                     "try `online_iNMF()` instead.")
                            }
                        } else if (isTRUE(readH5)) {
                            return(h5d[,])
                        } else {
                            stop("Can only set `readH5` to TRUE, FALSE, ",
                                 "or 'auto'.")
                        }
                    } else {
                        stop("H5 based dataset detected while `readH5` is ",
                             "set to FALSE.")
                    }
                } else {
                    return(scaleData(ld))
                }
            })
            out <- runINMF_R(
                object = data,
                k = k,
                lambda = lambda,
                thresh = thresh,
                maxIter = maxIter,
                nrep = nrep,
                H.init = H.init,
                W.init = W.init,
                V.init = V.init,
                useUnshared = FALSE,
                seed = seed,
                verbose = verbose
            )
            object@W <- out$W
            for (d in names(object)) {
                ld <- dataset(object, d)
                ld@H <- out$H[[d]]
                ld@V <- out$V[[d]]
                datasets(object, check = FALSE)[[d]] <- ld
            }
            object@uns$factorization$k <- k
            object@uns$factorization$lambda <- lambda
        } else {
            object <- runUINMF(
                object = object,
                k = k,
                lambda = lambda,
                thresh = thresh,
                maxIter = maxIter,
                nrep = nrep,
                seed = seed,
                verbose = verbose
            )
        }
        return(object)
    }
)

#' @rdname runINMF_R
#' @export
setMethod(
    "runINMF_R",
    signature(object = "list"),
    function(
        object,
        k,
        lambda = 5.0,
        maxIter = 30,
        nrep = 1,
        H.init = NULL,
        W.init = NULL,
        V.init = NULL,
        method = c("planc", "liger", "rcppml"),
        useUnshared = FALSE,
        seed = 1,
        readH5 = "auto",
        verbose = getOption("ligerVerbose")
    ) {
        # E ==> cell x gene scaled matrices
        E <- object
        nDatasets <- length(E)
        nCells <- sapply(E, ncol)
        nGenes <- nrow(E[[1]])
        if (k >= nGenes) {
            stop("Select k lower than the number of variable genes: ", nGenes)
        }
        Wm <- matrix(0, nGenes, k)
        Vm <- rep(list(matrix(0, nGenes, k)), nDatasets)
        Hm <- lapply(nCells, function(n) matrix(0, n, k))

        bestObj <- Inf
        bestSeed <- seed
        for (i in seq(nrep)) {
            set.seed(seed = seed + i - 1)
            startTime <- Sys.time()
            if (!is.null(W.init))
                W <- .checkInit(W.init, nCells, nGenes, k, "W")
            else W <- matrix(stats::runif(nGenes * k, 0, 2), nGenes, k)

            if (!is.null(V.init)) {
                V <- .checkInit(V.init, nCells, nGenes, k, "V")
            } else
                V <- lapply(seq(nDatasets), function(i) {
                    matrix(stats::runif(nGenes * k, 0, 2), nGenes, k)})

            if (!is.null(H.init)) {
                H <- .checkInit(H.init, nCells, nGenes, k, "H")
                H <- lapply(H, t)
            } else
                H <- lapply(nCells, function(n) {
                    matrix(stats::runif(n * k, 0, 2), n, k)
                })

            if (isTRUE(verbose)) {
                .log("Start iNMF with seed: ", seed + i - 1, "...")
                if (maxIter > 0)
                    pb <- utils::txtProgressBar(0, maxIter, style = 3)
            }
            iter <- 0
            while (iter < maxIter) {
                H <- inmfSolveH(W = W, V = V, E = E, lambda = lambda)
                V <- inmfSolveV(W = W, H = H, E = E, lambda = lambda)
                W <- inmfSolveW(H = H, V = V, E = E, lambda = lambda)
                iter <- iter + 1
                if (isTRUE(verbose) && maxIter > 0)
                    utils::setTxtProgressBar(pb, value = iter)

            }
            if (isTRUE(verbose) && maxIter > 0) {
                utils::setTxtProgressBar(pb, value = maxIter)
                cat("\n")
            }
            obj <- inmf_calcObj(E, H, W, V, lambda)
            if (obj < bestObj) {
                Wm <- W
                Hm <- H
                Vm <- V
                bestObj <- obj
                bestSeed <- seed + i - 1
            }
            endTime <- difftime(time1 = Sys.time(), time2 = startTime,
                                units = "auto")
            if (isTRUE(verbose)) {
                .log("Finished in ", endTime, " ", units(endTime),
                     "\nObjective error: ", bestObj)
                .log("Objective: ", obj)
                .log("Best results with seed ", bestSeed)
            }
        }
        out <- list(H = lapply(Hm, t), V = Vm, W = Wm)
        factorNames <- paste0("Factor_", seq(k))
        for (i in seq(nDatasets)) {
            dimnames(out$H[[i]]) <- list(factorNames, colnames(object[[i]]))
            dimnames(out$V[[i]]) <- list(rownames(object[[i]]), factorNames)
        }
        names(out$V) <- names(out$H) <- names(object)
        dimnames(out$W) <- list(rownames(object[[1]]), factorNames)
        return(out)
    }
)

inmf_calcObj <- function(E, H, W, V, lambda) {
    # E - dgCMatrix
    # H, W, V - matrix
    obj <- 0
    for (i in seq_along(E)) {
        obj <- obj +
            Matrix::norm(E[[i]] - (W + V[[i]]) %*% t(H[[i]]), "F") ^ 2 +
            lambda*norm(V[[i]] %*% t(H[[i]]), "F") ^ 2
    }
    return(obj)
}

inmfSolveH <- function(W, V, E, lambda) {
    H <- list()
    for (i in seq_along(E)) {
        CtC <- t(W + V[[i]]) %*% (W + V[[i]]) + lambda*(t(V[[i]]) %*% V[[i]])
        CtB <- as.matrix(t(W + V[[i]]) %*% E[[i]])
        H[[i]] <- t(RcppPlanc::bppnnls_prod(CtC, CtB))
    }
    return(H)
}

inmfSolveV <- function(W, H, E, lambda) {
    V <- list()
    for (i in seq_along(E)) {
        CtC <- (1 + lambda)*(t(H[[i]]) %*% H[[i]])
        CtB <- as.matrix(t(H[[i]]) %*% t(E[[i]]))
        CtB <- CtB - t(H[[i]]) %*% H[[i]] %*% t(W)
        V[[i]] <- t(RcppPlanc::bppnnls_prod(CtC, CtB))
    }
    return(V)
}

inmfSolveW <- function(H, V, E, lambda) {
    m <- nrow(E[[1]])
    k <- ncol(H[[1]])
    CtC <- matrix(0, k, k)
    CtB <- matrix(0, k, m)
    for (i in seq_along(E)) {
        CtC <- CtC + t(H[[i]]) %*% H[[i]]
        CtB <- CtB + as.matrix(t(H[[i]]) %*% t(E[[i]])) -
            t(H[[i]]) %*% H[[i]] %*% t(V[[i]])
    }
    return(t(RcppPlanc::bppnnls_prod(CtC, CtB)))
}
