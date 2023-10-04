#' Perform online iNMF on scaled datasets
#' @description Perform online integrative non-negative matrix factorization to
#' represent multiple single-cell datasets in terms of \eqn{H}, \eqn{W}, and
#' \eqn{V} matrices. It optimizes the iNMF objective function using online
#' learning (non-negative least squares for H matrix, hierarchical alternating
#' least squares for \eqn{W} and \eqn{V} matrices), where the number of factors
#' is set by \code{k}. The function allows online learning in 3 scenarios:
#'
#' \enumerate{
#'  \item Fully observed datasets;
#'  \item Iterative refinement using continually arriving datasets;
#'  \item Projection of new datasets without updating the existing factorization
#' }
#'
#' All three scenarios require fixed memory independent of the number of cells.
#'
#' For each dataset, this factorization produces an \eqn{H} matrix (k by cell),
#' a \eqn{V} matrix (genes by k) \eqn{C^\mathsf{T}C}, and a shared \eqn{W} matrix (genes by k). The
#' \eqn{H} matrices represent the cell factor loadings. \eqn{W} is identical
#' among all datasets, as it represents the shared components of the metagenes
#' across datasets. The \eqn{V} matrices represent the dataset-specific
#' components of the metagenes.
#' @details
#' For optional initialization, \code{W.init} must be a matrix object with
#' number of rows equal to number of variable genes (denoted as \code{g}) and
#' number of columns equal to \code{k}. Any of \code{V.init}, \code{A} and
#' \code{B} must be a list object of n matrices where n is the number of
#' datasets in \code{object}. For \code{V.init}, each matrix should be of size g
#' x k. For \code{A.init}, each matrix should be k x k and for \code{B.init},
#' each matrix should be g x k.
#'
#' Minibatch iterations is performed on small subset of cells. The exact
#' minibatch size applied on each dataset is \code{miniBatch_size} multiplied by
#' the proportion of cells in this dataset out of all cells. The setting of
#' \code{miniBatch_size} is by default \code{5000}, which is reasonable.
#' However, a smaller value such as \code{1000} may be necessary for analyzing
#' very small datasets. In general, \code{miniBatch_size} should be no larger
#' than the number of cells in the smallest dataset. An epoch is one completion
#' of calculation on all cells after a number of iterations of minibatches.
#' Therefore, the total number of iterations is determined by the setting of
#' \code{max.epochs}, total number of cells, and \code{miniBatch_size}.
#'
#' @param object \linkS4class{liger} object. Scaled data required.
#' @param newDatasets New datasets for scenario 2 or scenario 3.
#' @param projection Perform data integration with scenario 3. See description.
#' Default \code{FALSE}.
#' @param WInit Optional initialization for W. See detail. Default \code{NULL}.
#' @param VInit Optional initialization for V. See detail. Default \code{NULL}.
#' @param AInit Optional initialization for A. See detail. Default \code{NULL}.
#' @param BInit Optional initialization for B. See detail. Default \code{NULL}.
#' @param k Inner dimension of factorization--number of metagenes. A value in
#' the range 20-50 works well for most analyses. Default \code{20}.
#' @param lambda Regularization parameter. Larger values penalize
#' dataset-specific effects more strongly (i.e. alignment should increase as
#' lambda increases). We recommend always using the default value except
#' possibly for analyses with relatively small differences (biological
#' replicates, male/female comparisons, etc.) in which case a lower value such
#' as 1.0 may improve reconstruction quality. Default \code{5.0}.
#' @param maxEpochs The number of epochs to iterate through. See detail.
#' Default \code{5}.
#' @param HALSiters Maximum number of block coordinate descent (HALS
#' algorithm) iterations to perform for each update of \eqn{W} and \eqn{V}.
#' Default \code{1}. Changing this parameter is not recommended.
#' @param miniBatchSize Total number of cells in each minibatch. See detail.
#' Default \code{5000}.
#' @param seed Random seed to allow reproducible results. Default \code{123}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return \code{object} with \code{W} slot updated with resulting \eqn{W}
#' matrix; the \code{H}, \code{V}, \code{A} and \code{B} slots of each
#' \linkS4class{ligerDataset} object in \code{datasets} slot is updated with the
#' corresponding result matrices.
#' @export
#' @rdname runOnlineINMF
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' # Minibatch size has to be less than number of cell in the smallest dataset
#' # Scenario 1
#' pbmc <- online_iNMF(pbmc, miniBatch_size = 100)
#' # Scenario 2
#' # Fake new dataset by increasing all non-zero value in "ctrl" by 1
#' ctrl2 <- rawData(dataset(pbmc, "ctrl"))
#' ctrl2@x <- ctrl2@x + 1
#' colnames(ctrl2) <- paste0(colnames(ctrl2), 2)
#' pbmc2 <- online_iNMF(pbmc, k = 20, X_new = list(ctrl2 = ctrl2),
#'                      miniBatch_size = 100)
#' # Scenario 3
#' pbmc3 <- online_iNMF(pbmc, k = 20, X_new = list(ctrl2 = ctrl2),
#'                      miniBatch_size = 100, projection = TRUE)
runOnlineINMF <- function(
        object,
        k = 20,
        lambda = 5,
        newDatasets = NULL,
        projection = FALSE,
        maxEpochs = 5,
        HALSiter = 1,
        miniBatchSize = 5000,
        seed = 123,
        verbose = getOption("ligerVerbose"),
        ...
) {
    UseMethod("runOnlineINMF", object)
}

#' @export
#' @rdname runOnlineINMF
#' @method runOnlineINMF liger
runOnlineINMF.liger <- function(
        object,
        k = 20,
        lambda = 5,
        newDatasets = NULL,
        projection = FALSE,
        maxEpochs = 5,
        HALSiter = 1,
        miniBatchSize = 5000,
        seed = 123,
        verbose = getOption("ligerVerbose"),
        ...
) {
    .checkObjVersion(object)
    object <- recordCommand(object, dependencies = c("hdf5r", "RcppPlanc"))
    Es <- getMatrix(object, "scaleData", returnList = TRUE)
    Es <- lapply(datasets(object), function(ld) {
        sd <- scaleData(ld)
        if (is.null(sd))
            stop("Scaled data not available. Run `scaleNotCenter()` first")
        if (inherits(sd, "H5D")) return(.H5DToH5Mat(sd))
        else if (inherits(sd, "H5Group"))
            return(.H5GroupToH5SpMat(sd, c(length(varFeatures(object)), ncol(ld))))
        else return(scaleData(ld))
    })
    WInit <- VInit <- AInit <- BInit <- NULL
    if (!is.null(newDatasets)) {
        WInit <- getMatrix(object, "W")
        VInit <- getMatrix(object, "V")
        AInit <- getMatrix(object, "A")
        BInit <- getMatrix(object, "B")
        if (is.null(WInit) || any(sapply(VInit, is.null)) ||
            any(sapply(AInit, is.null)) || any(sapply(BInit, is.null))) {
            stop("Cannot find complete online iNMF result for current ",
                 "datasets. Please run `runOnlineINMF()` without `newDataset` ",
                 "first.")
        }
        # Put new datasets into input liger object
        if (is.list(newDatasets)) {
            allType <- sapply(newDatasets, function(x) class(x)[1])
            if (!all(allType == "dgCMatrix")) {
                stop("`newDatasets` must be a list of dgCMatrix for now.")
            }
            if (is.null(names(newDatasets))) {
                stop("`newDatasets` must be a named list.")
            }
            allNewNames <- character()
            for (i in seq_along(newDatasets)) {
                if (names(newDatasets)[i] %in% names(object)) {
                    newName <- paste0(names(newDatasets)[i], ".1")
                } else newName <- names(newDatasets)[i]
                dataset(object, newName) <- newDatasets[[i]]
                allNewNames <- c(allNewNames, newName)
            }
            object <- normalize(object, useDatasets = allNewNames,
                                verbose = verbose)
            object <- scaleNotCenter(object, useDatasets = allNewNames,
                                     verbose = verbose)
            newDatasets <- getMatrix(object, slot = "scaleData",
                                     dataset = allNewNames, returnList = TRUE)
        }
    }
    object <- closeAllH5(object)
    res <- runOnlineINMF.list(Es, newDatasets = newDatasets,
                              projection = projection, k = k, lambda = lambda,
                              maxEpochs = maxEpochs,
                              miniBatchSize = miniBatchSize,
                              HALSiter = HALSiter, verbose = verbose,
                              WInit = WInit, VInit = VInit, AInit = AInit,
                              BInit = BInit, seed = seed)

    for (i in seq_along(object)) {
        ld <- dataset(object, i)
        ld@H <- res$H[[i]]
        ld@V <- res$V[[i]]
        ld@A <- res$A[[i]]
        ld@B <- res$B[[i]]
        datasets(object, check = FALSE)[[i]] <- ld
    }
    object@W <- res$W
    object@uns$factorization <- list(k = k, lambda = lambda)
    suppressMessages({object <- restoreH5Liger(object)})
    return(object)
}

#' @export
#' @rdname runOnlineINMF
#' @method runOnlineINMF list
runOnlineINMF.list <- function(
        object,
        k = 20,
        lambda = 5,
        newDatasets = NULL,
        projection = FALSE,
        maxEpochs = 5,
        WInit = NULL,
        VInit = NULL,
        AInit = NULL,
        BInit = NULL,
        HALSiter = 1,
        miniBatchSize = 5000,
        seed = 123,
        verbose = getOption("ligerVerbose"),
        ...
) {
    if (!requireNamespace("RcppPlanc", quietly = TRUE))
        stop("RcppPlanc installation required. Currently, please get the ",
             "GitHub private repository access from the lab and run: \n",
             "devtools::install_github(\"welch-lab/RcppPlanc\")")
    nDatasets <- length(object) + length(newDatasets)
    barcodeList <- c(lapply(object, colnames), lapply(newDatasets, colnames))
    features <- rownames(object[[1]])
    if (!is.null(seed)) set.seed(seed)
    res <- RcppPlanc::onlineINMF(objectList = object, newDatasets = newDatasets,
                                 project = projection, k = k, lambda = lambda,
                                 maxEpoch = maxEpochs,
                                 minibatchSize = miniBatchSize,
                                 maxHALSIter = HALSiter, Vinit = VInit,
                                 Winit = WInit, Ainit = AInit, Binit = BInit,
                                 verbose = verbose)
    factorNames <- paste0("Factor_", seq(k))
    for (i in seq(nDatasets)) {
        res$H[[i]] <- t(res$H[[i]])
        dimnames(res$H[[i]]) <- list(factorNames, barcodeList[[i]])
        dimnames(res$V[[i]]) <- list(features, factorNames)
        dimnames(res$A[[i]]) <- list(factorNames, factorNames)
        dimnames(res$B[[i]]) <- list(features, factorNames)
    }
    names(res$A) <- names(res$B) <- names(res$V) <-
        names(res$H) <- c(names(object), names(newDatasets))
    dimnames(res$W) <- list(features, factorNames)
    return(res)
}

#' @export
#' @rdname runOnlineINMF
#' @method runOnlineINMF Seurat
runOnlineINMF.Seurat <- function(
        object,
        k = 20,
        lambda = 5,
        newDatasets = NULL,
        projection = FALSE,
        maxEpochs = 5,
        HALSiter = 1,
        miniBatchSize = 5000,
        seed = 123,
        verbose = getOption("ligerVerbose"),
        ...
) {

}
