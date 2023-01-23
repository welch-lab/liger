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
#' @param max.iters Maximum number of block coordinate descent iterations to
#' perform. Default \code{30}.
#' @param nrep Number of restarts to perform (iNMF objective function is
#' non-convex, so taking the best objective from multiple successive
#' initializations is recommended). For easier reproducibility, this increments
#' the random seed by 1 for each consecutive restart, so future factorizations
#' of the same dataset can be run with one rep if necessary. Default \code{1}.
#' @param H.init Initial values to use for \eqn{H} matrices. A list object where
#' each element is the initial \eqn{H} matrix of each dataset. Default
#' \code{NULL}.
#' @param W.init Initial values to use for \eqn{W} matrix. A matrix object.
#' Default \code{NULL}.
#' @param V.init Initial values to use for \eqn{V} matrices. A list object where
#' each element is the initial \eqn{V} matrix of each dataset. Default
#' \code{NULL}.
#' @param rand.seed Random seed to allow reproducible results. Default \code{1}.
#' @param print.obj Print objective function values after convergence when
#' \code{verbose = TRUE}. Default \code{FALSE}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{TRUE}.
#' @rdname optimizeALS
#' @export
setGeneric(
    "optimizeALS",
    function(
        object,
        k,
        lambda = 5.0,
        thresh = 1e-6,
        max.iters = 30,
        nrep = 1,
        H.init = NULL,
        W.init = NULL,
        V.init = NULL,
        use.unshared = FALSE,
        lamda.u = NULL,
        rand.seed = 1,
        print.obj = FALSE,
        readH5 = "auto",
        verbose = TRUE
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
        max.iters = 30,
        nrep = 1,
        H.init = NULL,
        W.init = NULL,
        V.init = NULL,
        use.unshared = FALSE,
        lamda.u = NULL,
        rand.seed = 1,
        print.obj = FALSE,
        readH5 = "auto",
        verbose = TRUE
    ) {
        if (isFALSE(use.unshared)) {
            object <- removeMissing(object, orient = "cell",
                                    verbose = verbose)
            data <- lapply(datasets(object), function(ld) {
                if (isH5Liger(ld)) {
                    if (!isFALSE(readH5)) {
                        h5d <- scale.data(ld)
                        if (readH5 == "auto") {
                            if (h5d$dims[2] <= 8000) {
                                readH5 <- TRUE
                                warning("Automatically reading H5 based ",
                                        "scaled dense matrix into memory. ",
                                        "Dim: ", h5d$dims[1], "x", h5d$dims[2],
                                        immediate. = verbose)
                            }
                            else {
                                stop("Scaled data in H5 based dataset with ",
                                     "more than 8000 cells will not be ",
                                     "automatically read into memory.")
                            }
                        }
                        if (isTRUE(readH5)) {
                            return(h5d[,])
                        }
                    } else {
                        stop("H5 based dataset detected while `readH5` is ",
                             "set to FALSE.")
                    }
                } else {
                    return(scale.data(ld))
                }
            })
            out <- optimizeALS(
                object = data,
                k = k,
                lambda = lambda,
                thresh = thresh,
                max.iters = max.iters,
                nrep = nrep,
                H.init = H.init,
                W.init = W.init,
                V.init = V.init,
                use.unshared = FALSE,
                rand.seed = rand.seed,
                print.obj = print.obj,
                verbose = verbose
            )
            colnames(out$W) <- var.features(object)
            object@W <- out$W
            names(out$H) <- names(out$V) <- names(object)
            for (d in names(object)) {
                ld <- dataset(object, d)
                rownames(out$H[[d]]) <- colnames(ld)
                ld@H <- out$H[[d]]
                colnames(out$V[[d]]) <- rownames(scale.data(ld))
                ld@V <- out$V[[d]]
                datasets(object, check = FALSE)[[d]] <- ld
            }
            #object@parameters$lambda <- lambda
            return(object)
        } else {
            object <- optimize_UANLS(
                object = object,
                k = k,
                lambda = lambda,
                thresh = thresh,
                max.iters = max.iters,
                nrep = nrep,
                rand.seed = rand.seed,
                print.obj = print.obj
            )
        }
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
        max.iters = 30,
        nrep = 1,
        H.init = NULL,
        W.init = NULL,
        V.init = NULL,
        use.unshared = FALSE,
        lamda.u = NULL,
        rand.seed = 1,
        print.obj = FALSE,
        verbose = TRUE
    ) {
        if (!all(sapply(object, is.matrix))) {
            stop("All values in 'object' must be a matrix")
        }
        E <- object
        E <- lapply(E, t)
        # N ==> nDataset
        nDatasets <- length(E)
        # TODO Fix the direction of scaled data. g x c or c x g?
        # ns ==> nCells
        nCells <- sapply(E, nrow)
        tmp <- gc()
        # g ==> nGenes
        nGenes <- ncol(E[[1]])
        if (k >= nGenes) {
            stop('Select k lower than the number of variable genes: ', g)
        }
        W_m <- matrix(0, k, nGenes)
        V_m <- lapply(seq(nDatasets), function(i) matrix(0, k, nGenes))
        H_m <- lapply(nCells, function(n) matrix(0, n, k))
        tmp <- gc()
        best_obj <- Inf
        run_stats <- matrix(0, nrep, 2)
        for (i in seq(nrep)) {
            set.seed(seed = rand.seed + i - 1)
            start_time <- Sys.time()
            if (!is.null(W.init)) W <- W.init
            else W <- matrix(runif(nGenes*k, 0, 2), k, nGenes)
            if (!is.null(V.init)) V <- V.init
            else
                V <- rep(list(matrix(runif(nGenes*k, 0, 2), k, nGenes)),
                         nDatasets)
            if (!is.null(H.init)) H <- H.init
            else H <- lapply(nCells, function(n) matrix(runif(n*k, 0, 2), n, k))
            tmp <- gc()
            delta <- 1
            iters <- 0
            sqrt_lambda <- sqrt(lambda)
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
                .log("Start iNMF...")
                pb <- utils::txtProgressBar(0, max.iters, style = 3)
            }
            while (delta > thresh & iters < max.iters) {
                H <- lapply(
                    seq(nDatasets),
                    function(i)
                        t(solveNNLS(
                            C = rbind(t(W + V[[i]]), sqrt_lambda*t(V[[i]])),
                            B = rbind(t(E[[i]]), matrix(0, nGenes, nCells[i]))
                        ))
                )
                tmp <- gc()
                V <- lapply(
                    seq(nDatasets),
                    function(i)
                        solveNNLS(C = rbind(H[[i]], sqrt_lambda*H[[i]]),
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
                utils::setTxtProgressBar(pb, value = max.iters)
                cat("\n")
            }
            # if (iters == max.iters) {
            #   print("Warning: failed to converge within the allowed number of iterations.
            #         Re-running with a higher max.iters is recommended.")
            # }
            if (obj < best_obj) {
                W_m <- W
                H_m <- H
                V_m <- V
                best_obj <- obj
                best_seed <- rand.seed + i - 1
            }
            end_time <- difftime(time1 = Sys.time(), time2 = start_time,
                                 units = "auto")
            run_stats[i, 1] <- as.double(end_time)
            run_stats[i, 2] <- iters
            if (isTRUE(verbose)) {
                .log("Finished in ", run_stats[i, 1], " ", units(end_time),
                     ", ", iters, " iterations. \nMax iterations set: ",
                     max.iters, "\nFinal objective delta: ", delta)
                if (isTRUE(print.obj)) .log("Objective: ", obj)
                .log("Best results with seed ", best_seed)
            }
        }
        out <- list()
        out$H <- H_m
        for (i in seq(nDatasets))
            rownames(out$H[[i]]) <- rownames(object[[i]])
        out$V <- V_m
        names(out$V) <- names(out$H) <- names(object)
        out$W <- W_m
        return(out)
    }
)
