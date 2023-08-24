#' Perform iNMF on scaled datasets
#' @description
#' Performs integrative non-negative matrix factorization (iNMF) (J.D. Welch,
#' 2019) to return factorized \eqn{H}, \eqn{W}, and \eqn{V} matrices. The
#' objective function is stated as
#'
#' \deqn{\arg\min_{H\ge0,W\ge0,V\ge0}\sum_{i}^{d}||E_i-(W+V_i)Hi||^2_F+\lambda\sum_{i}^{d}||V_iH_i||_F^2}
#'
#' where \eqn{E_i} is the input non-negative matrix of the i'th dataset, \eqn{d}
#' is the total number of datasets.
#'
#' The factorization produces a shared \eqn{W} matrix (genes by k), and for each
#' dataset, an \eqn{H} matrix (k by cells) and a \eqn{V} matrix (genes by k).
#' The \eqn{H} matrices represent the cell factor loadings. \eqn{W} is held
#' consistent among all datasets, as it represents the shared components of the
#' metagenes across datasets. The \eqn{V} matrices represent the
#' dataset-specific components of the metagenes.
#'
#' This function adopts highly optimized fast and memory efficient
#' implementation extended from Planc (Kannan, 2016). Pre-installation of
#' extension package \code{RcppPlanc} is required. The underlying algorithm
#' adopts the identical ANLS strategy as \code{\link{optimizeALS}} in the old
#' version of LIGER.
#' @param object A \linkS4class{liger} object, a Seurat object or a named list
#' of matrix, dgCMatrix, H5D objects, where the names represents dataset names
#' and matrices are scaled on the same set of variable features, with rows as
#' features and columns as cells.
#' @param k Inner dimension of factorization (number of factors). Run
#' \code{\link{suggestK}} to determine appropriate value; a general rule of
#' thumb is that a higher \code{k} will be needed for datasets with more
#' sub-structure.
#' @param lambda Regularization parameter. Larger values penalize
#' dataset-specific effects more strongly (i.e. alignment should increase as
#' \code{lambda} increases). Default \code{5}.
#' @param nIteration Total number of block coordinate descent iterations to
#' perform. Default \code{30}.
#' @param nRandomStarts Number of restarts to perform (iNMF objective function
#' is non-convex, so taking the best objective from multiple successive
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
#' @rdname runBPPINMF
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' pbmc <- runBPPINMF(pbmc, k = 20)
runBPPINMF <- function(
        object,
        k,
        lambda = 5.0,
        nIteration = 30,
        nRandomStarts = 1,
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
#' @return The liger method returns the input \linkS4class{liger} object with
#' factorization result updated. A list of all \eqn{H} matrices can be accessed
#' with \code{getMatrix(object, "H")}, a list of all \eqn{V} matrices can be
#' accessed with \code{getMatrix(object, "V")}, and the \eqn{W} matrix can be
#' accessed with \code{getMatrix(object, "W")}.
runBPPINMF.liger <- function(
        object,
        k,
        lambda = 5.0,
        nIteration = 30,
        nRandomStarts = 1,
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
        return(scaleData(ld))
    })
    dataClasses <- sapply(data, function(x) class(x)[1])
    if (!all(dataClasses == dataClasses[1])) {
        stop("Currently the scaledData of all datasets have to be of the same class.")
    }
    out <- runBPPINMF.list(
        object = data,
        k = k,
        lambda = lambda,
        nIteration = nIteration,
        nRandomStarts = nRandomStarts,
        HInit = HInit,
        WInit = WInit,
        VInit = VInit,
        seed = seed,
        verbose = verbose,
        barcodeList = lapply(datasets(object), colnames),
        features = varFeatures(object)
    )

    return(out)
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
#' @return The list method returns a list of entries \code{H}, \code{V} and
#' \code{W}. \code{H} is a list of \eqn{H} matrices for each dataset. \code{V}
#' is a list of \eqn{V} matrices for each dataset. \code{W} is the shared
#' \eqn{W} matrix.
#' @method runBPPINMF list
runBPPINMF.list <- function(
        object,
        k,
        lambda = 5.0,
        nIteration = 30,
        nRandomStarts = 1,
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
    datasetNames <- names(object)
    if (k >= nGenes) {
        stop("Select k lower than the number of variable genes: ", nGenes)
    }
    Wm <- matrix(0, k, nGenes)
    Vm <- rep(list(matrix(0, k, nGenes)), nDatasets)
    Hm <- lapply(nCells, function(n) matrix(0, n, k))

    bestObj <- Inf
    bestSeed <- seed
    for (i in seq(nRandomStarts)) {
        if (isTRUE(verbose)) {
            .log("Replicate run ", i, "...")
        }
        set.seed(seed = seed + i - 1)
        if (!is.null(WInit))
            W <- t(.checkInit(WInit, nCells, nGenes, k, "W"))
        else W <- t(matrix(stats::runif(nGenes * k, 0, 2), k, nGenes))
        if (!is.null(VInit)) {
            V <- .checkInit(VInit, nCells, nGenes, k, "V")
            V <- lapply(V, t)
        } else {
            V <- lapply(seq(nDatasets), function(i) {
                t(matrix(stats::runif(nGenes * k, 0, 2), k, nGenes))})
        }
        if (!is.null(HInit)) {
            H <- .checkInit(HInit, nCells, nGenes, k, "H")
            H <- lapply(H, t)
        } else {
            H <- lapply(nCells, function(n) {
                matrix(stats::runif(n * k, 0, 2), n, k)
            })
        }
        if (inherits(object[[1]], "H5D")) {
            # RcppPlanc::bppinmf_h5dense()
            stop("TODO: Push Yichen to test bppinmf_h5sparse/bppinmf_h5dense!")
        } else {
            out <- RcppPlanc::bppinmf(
                objectList = object, k = k, lambda = lambda, niter = nIteration,
                verbose = verbose)
        }

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
    return(list(H = Hm, V = Vm, W = Wm))
}

#' @rdname runBPPINMF
#' @export
#' @param datasetVar Variable name in metadata indicating a factor of dataset
#' belonging, or directly a factor that match with the number of cells.
#' @return The Seurat method returns a list of entries \code{H}, \code{V} and
#' \code{W}. \code{H} is a list of \eqn{H} matrices for each dataset. \code{V}
#' is a list of \eqn{V} matrices for each dataset. \code{W} is the shared
#' \eqn{W} matrix.
#' @method runBPPINMF Seurat
runBPPINMF.Seurat <- function(
        object,
        datasetVar,
        k,
        lambda = 5.0,
        nIteration = 30,
        nRandomStarts = 1,
        HInit = NULL,
        WInit = NULL,
        VInit = NULL,
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
        stop("Invalid `datasetVar`. Please see `?runBPPINMF` for instruction.")
    }
    datasetVar <- droplevels(datasetVar)
    Es <- lapply(levels(datasetVar), function(d) {
        as(EBind[, datasetVar == d], "CsparseMatrix")
    })
    names(Es) <- levels(datasetVar)
    runBPPINMF.list(
        object = Es,
        k = k,
        lambda = lambda,
        nIteration = nIteration,
        nRandomStarts = nRandomStarts,
        HInit = HInit,
        WInit = WInit,
        VInit = VInit,
        seed = seed,
        verbose = verbose
    )
}
