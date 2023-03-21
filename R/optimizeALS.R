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
#' perform. Default \code{10}.
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
#' @rdname optimizeALS
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' # Only running a few iterations for fast examples
#' pbmc <- optimizeALS(pbmc, k = 20, maxIter = 2)
setGeneric(
    "optimizeALS",
    function(
        object,
        k,
        lambda = 5.0,
        thresh = 1e-6,
        maxIter = 10,
        nrep = 1,
        H.init = NULL,
        W.init = NULL,
        V.init = NULL,
        useUnshared = FALSE,
        seed = 1,
        readH5 = "auto",
        verbose = getOption("ligerVerbose"),
        # Deprecated coding style
        max.iters = maxIter,
        use.unshared = useUnshared,
        rand.seed = seed,
        # Deprecated functionality
        print.obj = NULL
    ) standardGeneric("optimizeALS")
)

#' @rdname optimizeALS
#' @export
setMethod(
    "optimizeALS",
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
        verbose = getOption("ligerVerbose"),
        # Deprecated coding style
        max.iters = maxIter,
        use.unshared = useUnshared,
        rand.seed = seed,
        # Deprecated functionality
        print.obj = NULL
    ) {
        .deprecateArgs(list(max.iters = "maxIter", use.unshared = "useUnshared",
                            rand.seed = "seed"), defunct = "print.obj")
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
            out <- optimizeALS(
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
            rownames(object@W) <- varFeatures(object)
            for (d in names(object)) {
                ld <- dataset(object, d)
                ld@H <- out$H[[d]]
                colnames(ld@H) <- colnames(ld)
                ld@V <- out$V[[d]]
                rownames(ld@V) <- varFeatures(object)
                datasets(object, check = FALSE)[[d]] <- ld
            }
            object@uns$factorization$k <- k
            object@uns$factorization$lambda <- lambda
        } else {
            object <- optimizeUANLS(
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

#' @rdname optimizeALS
#' @export
setMethod(
    "optimizeALS",
    signature(object = "list"),
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
        verbose = getOption("ligerVerbose"),
        # Deprecated coding style
        max.iters = maxIter,
        use.unshared = useUnshared,
        rand.seed = seed,
        # Deprecated functionality
        print.obj = NULL
    ) {
        .deprecateArgs(list(max.iters = "maxIter", use.unshared = "useUnshared",
                            rand.seed = "seed"), defunct = "print.obj")
        if (!all(sapply(object, is.matrix))) {
            stop("All values in 'object' must be a matrix")
        }
        # E ==> cell x gene scaled matrices
        E <- lapply(object, t)
        nDatasets <- length(E)
        nCells <- sapply(E, nrow)
        tmp <- gc() # nolint
        nGenes <- ncol(E[[1]])
        if (k >= nGenes) {
            stop("Select k lower than the number of variable genes: ", nGenes)
        }
        Wm <- matrix(0, k, nGenes)
        Vm <- rep(list(matrix(0, k, nGenes)), nDatasets)
        Hm <- lapply(nCells, function(n) matrix(0, n, k))
        tmp <- gc()
        bestObj <- Inf
        bestSeed <- seed
        runStats <- matrix(0, nrep, 2)
        for (i in seq(nrep)) {
            set.seed(seed = seed + i - 1)
            startTime <- Sys.time()
            if (!is.null(W.init))
                W <- t(.checkInit(W.init, nCells, nGenes, k, "W"))
            else W <- matrix(stats::runif(nGenes * k, 0, 2), k, nGenes)
            if (!is.null(V.init)) {
                V <- .checkInit(V.init, nCells, nGenes, k, "V")
                V <- lapply(V, t)
            } else
                V <- lapply(seq(nDatasets), function(i) {
                    matrix(stats::runif(nGenes * k, 0, 2), k, nGenes)})
            if (!is.null(H.init)) {
                H <- .checkInit(H.init, nCells, nGenes, k, "H")
                H <- lapply(H, t)
            } else
                H <- lapply(nCells, function(n) {
                    matrix(stats::runif(n * k, 0, 2), n, k)
                })
            tmp <- gc()
            delta <- 1
            iters <- 0
            sqrtLambda <- sqrt(lambda)

            obj0 <- sum(sapply(
                seq(nDatasets),
                function(i) norm(E[[i]] - H[[i]] %*% (W + V[[i]]), "F") ^ 2
            )) +
                sum(sapply(
                    seq(nDatasets),
                    function(i) lambda*norm(H[[i]] %*% V[[i]], "F") ^ 2
                ))
            tmp <- gc()
            if (isTRUE(verbose)) {
                .log("Start iNMF with seed: ", seed + i - 1, "...")
                pb <- utils::txtProgressBar(0, maxIter, style = 3)
            }
            while (delta > thresh & iters < maxIter) {
                H <- lapply(
                    seq(nDatasets),
                    function(i)
                        t(solveNNLS(
                            C = rbind(t(W + V[[i]]), sqrtLambda*t(V[[i]])),
                            B = rbind(t(E[[i]]), matrix(0, nGenes, nCells[i]))
                        ))
                )
                tmp <- gc()
                V <- lapply(
                    seq(nDatasets),
                    function(i)
                        solveNNLS(C = rbind(H[[i]], sqrtLambda*H[[i]]),
                            B = rbind(E[[i]] - H[[i]] %*% W,
                                      matrix(0, nCells[[i]], nGenes)))
                )
                tmp <- gc()
                W <- solveNNLS(C = rbindlist(H),
                               B = rbindlist(lapply(seq(nDatasets),
                                   function(i) E[[i]] - H[[i]] %*% V[[i]]
                               )))
                tmp <- gc()
                obj <- sum(sapply(
                    seq(nDatasets),
                    function(i) norm(E[[i]] - H[[i]] %*% (W + V[[i]]), "F") ^ 2
                )) +
                    sum(sapply(
                        seq(nDatasets),
                        function(i)
                            lambda*norm(H[[i]] %*% V[[i]], "F") ^ 2
                    ))
                tmp <- gc()
                delta <- abs(obj0 - obj) / (mean(obj0, obj))
                obj0 <- obj
                iters <- iters + 1
                if (isTRUE(verbose))
                    utils::setTxtProgressBar(pb, value = iters)
            }
            if (isTRUE(verbose)) {
                utils::setTxtProgressBar(pb, value = maxIter)
                cat("\n")
            }
            # if (iters == maxIter) {
            #   print("Warning: failed to converge within the allowed number of iterations.
            #         Re-running with a higher maxIter is recommended.")
            # }
            if (obj < bestObj) {
                Wm <- W
                Hm <- H
                Vm <- V
                bestObj <- obj
                bestSeed <- seed + i - 1
            }
            endTime <- difftime(time1 = Sys.time(), time2 = startTime,
                                 units = "auto")
            runStats[i, 1] <- as.double(endTime)
            runStats[i, 2] <- iters
            if (isTRUE(verbose)) {
                .log("Finished in ", runStats[i, 1], " ", units(endTime),
                     ", ", iters, " iterations. \nMax iterations set: ",
                     maxIter, "\nFinal objective delta: ", delta)
                .log("Objective: ", obj)
                .log("Best results with seed ", bestSeed)
            }
        }
        out <- list(H = Hm, V = Vm, W = t(Wm))
        factorNames <- paste0("Factor_", seq(k))
        for (i in seq(nDatasets)) {
            out$H[[i]] <- t(out$H[[i]])
            colnames(out$H[[i]]) <- colnames(object[[i]])
            rownames(out$H[[i]]) <- factorNames
            out$V[[i]] <- t(out$V[[i]])
            rownames(out$V[[i]]) <- rownames(object[[i]])
            colnames(out$V[[i]]) <- factorNames
        }
        names(out$V) <- names(out$H) <- names(object)
        rownames(out$W) <- rownames(object[[1]])
        colnames(out$W) <- factorNames
        return(out)
    }
)

# Binds list of matrices row-wise (vertical stack)
rbindlist <- function(mat_list) do.call(rbind, mat_list)
