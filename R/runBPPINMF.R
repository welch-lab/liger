#' Perform iNMF on scaled datasets
#' @description
#' Performs integrative non-negative matrix (iNMF) factorization to return
#' factorized \eqn{H}, \eqn{W}, and \eqn{V} matrices, using highly optimized
#' fast and memory efficient implementation extended from Planc (Kannan, 2016).
#' Pre-installation of extension package \code{RcppPlanc} is required. The
#' underlying algorithm adopts the identical ANLS strategy as
#' \code{\link{optimizeALS}} in the old version of LIGER.
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
#' @param HInit Initial values to use for \eqn{H} matrices. A list object where
#' each element is the initial \eqn{H} matrix of each dataset. Default
#' \code{NULL}.
#' @param WInit Initial values to use for \eqn{W} matrix. A matrix object.
#' Default \code{NULL}.
#' @param VInit Initial values to use for \eqn{V} matrices. A list object where
#' each element is the initial \eqn{V} matrix of each dataset. Default
#' \code{NULL}.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param ... Arguments passed to methods.
#' @return \code{object} with \code{W} slot updated with the result \eqn{W}
#' matrix, and the \code{H} and \code{V} slots of each
#' \linkS4class{ligerDataset} object in the \code{datasets} slot updated with
#' the dataset specific \eqn{H} and \eqn{V} matrix, respectively.
#' @rdname runBPPINMF
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' # Only running a few iterations for fast examples
#' pbmc <- runINMF(pbmc, k = 20, maxIter = 2)
runBPPINMF <- function(
        object,
        k,
        lambda = 5.0,
        thresh = 1e-6,
        maxIter = 30,
        nrep = 1,
        HInit = NULL,
        WInit = NULL,
        VInit = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    UseMethod("runBPPINMF", object)
}

#' @rdname runBPPINMF
#' @export
#' @param readH5 \code{TRUE} to force reading H5 based data into memory and
#' conduct factorization. \code{"auto"} reads H5 dataset with less than 8000
#' cells. \code{FALSE} will stop users from running if H5 data presents.
#' @method runBPPINMF liger
runBPPINMF.liger <- function(
        object,
        k,
        lambda = 5.0,
        thresh = 1e-6,
        maxIter = 30,
        nrep = 1,
        HInit = NULL,
        WInit = NULL,
        VInit = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        readH5 = "auto",
        ...
) {
    .checkObjVersion(object)
    object <- recordCommand(object)
    object <- removeMissing(object, orient = "cell", verbose = verbose)
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
    out <- runBPPINMF.list(
        object = data,
        k = k,
        lambda = lambda,
        thresh = thresh,
        maxIter = maxIter,
        nrep = nrep,
        HInit = HInit,
        WInit = WInit,
        VInit = VInit,
        seed = seed,
        verbose = verbose,
        barcodeList = lapply(datasets(object), colnames),
        features = varFeatures(object)
    )
    # return(out)
    object@W <- out$W
    for (d in names(object)) {
        ld <- dataset(object, d)
        ld@H <- out$H[[d]]
        ld@V <- out$V[[d]]
        datasets(object, check = FALSE)[[d]] <- ld
    }
    object@uns$factorization$k <- k
    object@uns$factorization$lambda <- lambda
    return(object)
}

#' @rdname runBPPINMF
#' @export
#' @param barcodeList List object of barcodes for each datasets, for setting
#' dimnames of output \eqn{H} matrices. Default \code{NULL} uses \code{colnames}
#' of matrices in the \code{object}.
#' @param features Character vector of feature names, for setting dimnames of
#' output \eqn{V} and \eqn{W} matrices. Default \code{NULL} uses \code{rownames}
#' of matrices in the \code{object}.
#' @method runBPPINMF list
runBPPINMF.list <- function(
        object,
        k,
        lambda = 5.0,
        thresh = 1e-6,
        maxIter = 30,
        nrep = 1,
        HInit = NULL,
        WInit = NULL,
        VInit = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        barcodeList = NULL,
        features = NULL,
        ...
) {
    if (!requireNamespace("RcppPlanc", quietly = TRUE))
        stop("RcppPlanc installation required. Currently, please get the ",
             "GitHub private repository access from the lab and run: \n",
             "devtools::install_github(\"welch-lab/RcppPlanc\")")

    if (is.null(barcodeList)) {
        barcodeList <- lapply(object, colnames)
    } else if (!identical(sapply(barcodeList, length), sapply(object, ncol))) {
        stop("Given `barcodeList` cannot match to columns of all matrices.")
    }

    if (is.null(features)) features <- rownames(object[[1]])
    else {
        if (!identical(length(features), unique(sapply(object, nrow)))) {
            stop("Either the length of `features` does not match with the ",
                 "rows of all matrices, or the matrices have varied numbers ",
                 "of rows.")
        }
    }

    nCells <- sapply(object, ncol)
    nGenes <- nrow(object[[1]])
    nDatasets <- length(object)
    if (k >= nGenes) {
        stop("Select k lower than the number of variable genes: ", nGenes)
    }
    Wm <- matrix(0, k, nGenes)
    Vm <- rep(list(matrix(0, k, nGenes)), nDatasets)
    Hm <- lapply(nCells, function(n) matrix(0, n, k))

    bestObj <- Inf
    bestSeed <- seed
    for (i in seq(nrep)) {
        if (isTRUE(verbose)) {
            .log("Replicate run ", i, "...")
        }
        set.seed(seed = seed + i - 1)
        if (!is.null(WInit))
            W <- t(.checkInit(WInit, nCells, nGenes, k, "W"))
        else W <- matrix(stats::runif(nGenes * k, 0, 2), k, nGenes)
        if (!is.null(VInit)) {
            V <- .checkInit(VInit, nCells, nGenes, k, "V")
            V <- lapply(V, t)
        } else
            V <- lapply(seq(nDatasets), function(i) {
                matrix(stats::runif(nGenes * k, 0, 2), k, nGenes)})
        if (!is.null(HInit)) {
            H <- .checkInit(HInit, nCells, nGenes, k, "H")
            H <- lapply(H, t)
        } else
            H <- lapply(nCells, function(n) {
                matrix(stats::runif(n * k, 0, 2), n, k)
            })

        out <- RcppPlanc::bppinmf_sparse(
            objectList = object, k = k, lambda = lambda, maxIter = maxIter,
            thresh = thresh, verbose = verbose, Hinit = HInit, Winit = WInit,
            Vinit = VInit
        )
        # return(out)
        if (out$objErr < bestObj) {
            Wm <- out$W
            Hm <- out$H
            Vm <- out$V
            bestObj <- out$objErr
            bestSeed <- seed + i - 1
        }
    }
    if (isTRUE(verbose)) {
        .log("Best objective error: ", bestObj, "\nBest seed: ", bestSeed)
    }
    factorNames <- paste0("Factor_", seq(k))
    for (i in seq(nDatasets)) {
        Hm[[i]] <- t(Hm[[i]])
        colnames(Hm[[i]]) <- barcodeList[[i]]
        rownames(Hm[[i]]) <- factorNames
        rownames(Vm[[i]]) <- features
        colnames(Vm[[i]]) <- factorNames
    }
    names(Vm) <- names(Hm) <- names(object)
    rownames(Wm) <- features
    colnames(Wm) <- factorNames
    return(out)
}
