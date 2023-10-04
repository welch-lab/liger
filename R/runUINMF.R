#' Perform iNMF on scaled datasets, and include unshared features
#' @param object \linkS4class{liger} object. Should run
#' \code{\link{selectGenes}} with \code{unshared = TRUE} and then run
#' \code{\link{scaleNotCenter}} in advance.
#' @param k Integer, inner dimension of factorization (number of factors).
#' Default \code{20}.
#' @param lambda Numeric, the regularization parameter. Default \code{5}.
#' @param nIteration Total number of iterations to perform. Default \code{30}.
#' @param nRandomStarts Number of restarts to perform (iNMF objective function
#' is non-convex, so taking the best objective from multiple successive
#' initializations is recommended). For easier reproducibility, this increments
#' the random seed by 1 for each consecutive restart, so future factorizations
#' of the same dataset can be run with one rep if necessary. Default \code{1}.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @export
#' @rdname runUINMF
runUINMF <- function(
        object,
        k = 20,
        lambda = 5,
        nIteration = 30,
        nRandomStarts = 1,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    UseMethod("runUINMF", object)
}

#' @export
#' @rdname runUINMF
#' @method runUINMF liger
runUINMF.liger <- function(
        object,
        k = 20,
        lambda = 5,
        nIteration = 30,
        nRandomStarts = 1,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    .checkObjVersion(object)
    object <- recordCommand(object, dependencies = "RcppPlanc")
    object <- removeMissing(object, orient = "cell", verbose = verbose)
    # Elist <- getMatrix(object, "scaleData", returnList = TRUE)

    Elist <- lapply(datasets(object), function(ld) {
        if (is.null(scaleData(ld)))
            stop("Scaled data not available. ",
                 "Run `scaleNotCenter(object)` first")
        return(scaleData(ld))
    })
    Ulist <- getMatrix(object, "scaleUnsharedData", returnList = TRUE)
    if (all(sapply(Ulist, is.null))) {
        stop("No scaled data for unshared feature found. Run `selectGenes()` ",
             "with `unshared = TRUE` and then `scaleNotCenter()`.")
    }
    res <- runUINMF.list(Elist, Ulist, k = k, lambda = lambda,
                         nIteration = nIteration, nRandomStarts = nRandomStarts,
                         seed = seed, verbose = verbose, ...)
    for (i in seq_along(object)) {
        ld <- dataset(object, i)
        ld@H <- res$H[[i]]
        ld@V <- res$V[[i]]
        if (!is.null(ld@scaleUnsharedData)) ld@U <- res$U[[i]]
        datasets(object, check = FALSE)[[i]] <- ld
    }
    object@W <- res$W
    object@uns$factorization <- list(k = k, lambda = lambda)
    return(object)
}

#' @export
#' @rdname runUINMF
#' @method runUINMF list
#' @param unsharedList List of matrices for unshared features
runUINMF.list <- function(
        object,
        unsharedList,
        k = 20,
        lambda = 5,
        nIteration = 30,
        nRandomStarts = 1,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    bestObj <- Inf
    bestRes <- NULL
    bestSeed <- NULL
    for (i in seq(nRandomStarts)) {
        seed <- seed + i - 1
        set.seed(seed)
        res <- RcppPlanc::uinmf(object, unsharedList, k = k, lambda = lambda,
                                niter = nIteration, verbose = verbose)
        if (res$objErr < bestObj) {
            bestRes <- res
            bestObj <- res$objErr
            bestSeed <- seed
        }
    }
    if (isTRUE(verbose) && nRandomStarts > 1) {
        .log("Best objective error: ", bestObj, "\nBest seed: ", bestSeed)
    }
    rm(res)
    features <- rownames(object[[1]])
    unsharedFeatures <- lapply(unsharedList, rownames)
    factorNames <- paste0("Factor_", seq(k))
    barcodes <- lapply(object, colnames)
    for (i in seq_along(object)) {
        bestRes$H[[i]] <- t(bestRes$H[[i]])
        dimnames(bestRes$H[[i]]) <- list(factorNames, barcodes[[i]])
        dimnames(bestRes$V[[i]]) <- list(features, factorNames)
        dimnames(bestRes$U[[i]]) <- list(unsharedFeatures[[i]], factorNames)
        rownames(bestRes$U[[i]]) <- unsharedFeatures[[i]]
        colnames(bestRes$U[[i]]) <- factorNames
    }
    dimnames(bestRes$W) <- list(features, factorNames)
    return(bestRes)
}

#' @export
#' @rdname runUINMF
#' @method runUINMF Seurat
#' @param datasetVar Variable for dataset belonging.
runUINMF.Seurat <- function(
        object,
        unsharedList,
        datasetVar,
        k = 20,
        lambda = 5,
        nIteration = 30,
        nRandomStarts = 1,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat installation required. Please run\n",
             "install.packages(\"Seurat\")")
    }
    EBind <- Seurat::GetAssayData(object, "scale.data")
    if (any(EBind < 0)) {
        stop("Non-negative Matrix Factorization requires non-negative data. ",
             "Please scale the library-size-normalized data without centering.")
    }
    if (is.character(datasetVar) && length(datasetVar) == 1) {
        datasetVar <- object[[datasetVar]][[1]]
    }
    if (!is.factor(datasetVar) || length(datasetVar) != ncol(EBind)) {
        stop("Invalid `datasetVar`. Please see `?runUINMF` for instruction.")
    }
    datasetVar <- droplevels(datasetVar)
    Es <- lapply(levels(datasetVar), function(d) {
        as(EBind[, datasetVar == d], "CsparseMatrix")
    })
    names(Es) <- levels(datasetVar)
    runUINMF.list(
        object = Es,
        unsharedList = unsharedList,
        k = k,
        lambda = lambda,
        nIteration = nIteration,
        nRandomStarts = nRandomStarts,
        seed = seed,
        verbose = verbose,
        ...
    )
}

