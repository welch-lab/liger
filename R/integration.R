#' Integrate scaled datasets with iNMF or variant methods
#' @description
#' LIGER provides dataset integration methods based on iNMF (integrative
#' Non-negative Matrix Factorization [1]) and its variants (online iNMF [2] and
#' UINMF [3]). This function wraps \code{\link{runINMF}},
#' \code{\link{runOnlineINMF}} and \code{\link{runUINMF}}, of which the help
#' pages have more detailed description.
#' @param object A \linkS4class{liger} object or a Seurat object with
#' non-negative scaled data of variable features (Done with
#' \code{\link{scaleNotCenter}}).
#' @param k Inner dimension of factorization (number of factors). Generally, a
#' higher \code{k} will be needed for datasets with more sub-structure. Default
#' \code{20}.
#' @param lambda Regularization parameter. Larger values penalize
#' dataset-specific effects more strongly (i.e. alignment should increase as
#' \code{lambda} increases). Default \code{5}.
#' @param method iNMF variant algorithm to use for integration. Choose from
#' \code{"iNMF"}, \code{"onlineINMF"}, \code{"UINMF"}. Default \code{"iNMF"}.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param ... Arguments passed to other methods and wrapped functions.
#' @export
#' @rdname runIntegration
#' @return Updated input object. For detail, please refer to the refered method
#' linked in Description.
#' @references
#' \enumerate{
#'  \item{Joshua D. Welch and et al., Single-Cell Multi-omic Integration Compares
#' and Contrasts Features of Brain Cell Identity, Cell, 2019}
#'  \item{Chao Gao and et al., Iterative single-cell multi-omic integration using
#' online learning, Nat Biotechnol., 2021}
#'  \item{April R. Kriebel and Joshua D. Welch, UINMF performs mosaic integration
#' of single-cell multi-omic datasets using nonnegative matrix factorization,
#' Nat. Comm., 2022}
#' }
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' if (requireNamespace("RcppPlanc", quietly = TRUE)) {
#'     pbmc <- runIntegration(pbmc)
#' }
runIntegration <- function(
        object,
        k = 20,
        lambda = 5.0,
        method = c("iNMF", "onlineINMF", "UINMF"),
        ...
) {
    UseMethod("runIntegration", object)
}

#' @rdname runIntegration
#' @export
#' @method runIntegration liger
runIntegration.liger <- function(
        object,
        k = 20,
        lambda = 5.0,
        method = c("iNMF", "onlineINMF", "UINMF"),
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    method <- match.arg(method)
    object <- switch(
        method,
        iNMF = runINMF(object, k = k, lambda = lambda, seed = seed,
                       verbose = verbose, ...),
        onlineINMF = runOnlineINMF(object, k = k, lambda = lambda, seed = seed,
                                   verbose = verbose, ...),
        UINMF = runUINMF(object = object, k = k, lambda = lambda, seed = seed,
                         verbose = verbose, ...)
    )
    return(object)
}

#' @rdname runIntegration
#' @export
#' @method runIntegration Seurat
#' @param datasetVar Metadata variable name that stores the dataset source
#' annotation. Default \code{"orig.ident"}.
#' @param useLayer For Seurat>=4.9.9, the name of layer to retrieve input
#' non-negative scaled data. Default \code{"ligerScaleData"}. For older Seurat,
#' always retrieve from \code{scale.data} slot.
#' @param assay Name of assay to use. Default \code{NULL} uses current active
#' assay.
runIntegration.Seurat <- function(
        object,
        k = 20,
        lambda = 5.0,
        method = c("iNMF", "onlineINMF"),
        datasetVar = "orig.ident",
        useLayer = "ligerScaleData",
        assay = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    method <- match.arg(method)
    object <- switch(
        method,
        iNMF = runINMF(object, k = k, lambda = lambda, seed = seed,
                       useLayer = useLayer, assay = assay,
                       datasetVar = datasetVar, verbose = verbose, ...),
        onlineINMF = runOnlineINMF(object, k = k, lambda = lambda, seed = seed,
                                   useLayer = useLayer, assay = assay,
                                   datasetVar = datasetVar,
                                   verbose = verbose, ...)
    )
    return(object)
}


############################### regular INMF ###################################

#' Perform iNMF on scaled datasets
#' @description
#' Performs integrative non-negative matrix factorization (iNMF) (J.D. Welch,
#' 2019) using block coordinate descent (alternating non-negative
#' least squares, ANLS) to return factorized \eqn{H}, \eqn{W}, and \eqn{V}
#' matrices. The objective function is stated as
#'
#' \deqn{\arg\min_{H\ge0,W\ge0,V\ge0}\sum_{i}^{d}||E_i-(W+V_i)Hi||^2_F+
#' \lambda\sum_{i}^{d}||V_iH_i||_F^2}
#'
#' where \eqn{E_i} is the input non-negative matrix of the i'th dataset, \eqn{d}
#' is the total number of datasets. \eqn{E_i} is of size \eqn{m \times n_i} for
#' \eqn{m} variable genes and \eqn{n_i} cells, \eqn{H_i} is of size
#' \eqn{n_i \times k}, \eqn{V_i} is of size \eqn{m \times k}, and \eqn{W} is of
#' size \eqn{m \times k}.
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
#' @section Difference from optimizeALS():
#' In the old version implementation, we compute the objective error at the end
#' of each iteration, and then compares if the algorithm is reaching a
#' convergence, using an argument \code{thresh}. Now, since the computation of
#' objective error is indeed expensive, we canceled this feature and directly
#' runs a default of 30 (\code{nIteration}) iterations, which empirically leads
#' to a convergence most of the time. Given that the new version is highly
#' optimized, running this many iteration should be acceptable.
#' @param object A \linkS4class{liger} object or a Seurat object with
#' non-negative scaled data of variable features (Done with
#' \code{\link{scaleNotCenter}}).
#' @param k Inner dimension of factorization (number of factors). Generally, a
#' higher \code{k} will be needed for datasets with more sub-structure. Default
#' \code{20}.
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
#' @rdname runINMF
#' @export
#' @return
#' \itemize{
#'  \item{liger method - Returns updated input \linkS4class{liger} object
#'  \itemize{
#'      \item{A list of all \eqn{H} matrices can be accessed with
#'          \code{getMatrix(object, "H")}}
#'      \item{A list of all \eqn{V} matrices can be accessed with
#'          \code{getMatrix(object, "V")}}
#'      \item{The \eqn{W} matrix can be accessed with
#'          \code{getMatrix(object, "W")}}
#'  }}
#'  \item{Seurat method - Returns updated input Seurat object
#'  \itemize{
#'      \item{\eqn{H} matrices for all datasets will be concatenated and
#'          transposed (all cells by k), and form a DimReduc object in the
#'          \code{reductions} slot named by argument \code{reduction}.}
#'      \item{\eqn{W} matrix will be presented as \code{feature.loadings} in the
#'          same DimReduc object.}
#'      \item{\eqn{V} matrices, an objective error value and the dataset
#'          variable used for the factorization is currently stored in
#'          \code{misc} slot of the same DimReduc object.}
#'  }}
#' }
#' @references Joshua D. Welch and et al., Single-Cell Multi-omic Integration
#' Compares and Contrasts Features of Brain Cell Identity, Cell, 2019
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' if (requireNamespace("RcppPlanc", quietly = TRUE)) {
#'     pbmc <- runINMF(pbmc)
#' }
runINMF <- function(
        object,
        k = 20,
        lambda = 5.0,
        ...
) {
    UseMethod("runINMF", object)
}

#' @rdname runINMF
#' @export
#' @method runINMF liger
runINMF.liger <- function(
        object,
        k = 20,
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
    .checkObjVersion(object)
    object <- recordCommand(object, ..., dependencies = "RcppPlanc")
    object <- removeMissing(object, orient = "cell", verbose = verbose)
    data <- lapply(datasets(object), function(ld) {
        if (is.null(scaleData(ld)))
            stop("Scaled data not available. ",
                 "Run `scaleNotCenter(object)` first")
        return(scaleData(ld))
    })
    dataClasses <- sapply(data, function(x) class(x)[1])
    if (!all(dataClasses == dataClasses[1])) {
        stop("Currently the scaledData of all datasets have to be of the ",
             "same class.")
    }
    out <- .runINMF.list(
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

    object@W <- out$W
    rownames(object@W) <- varFeatures(object)
    for (d in names(object)) {
        ld <- dataset(object, d)
        ld@H <- out$H[[d]]
        ld@V <- out$V[[d]]
        if (isH5Liger(ld)) {
            colnames(ld@H) <- colnames(ld)
            rownames(ld@V) <- varFeatures(object)
        }
        datasets(object, check = FALSE)[[d]] <- ld
    }
    object@uns$factorization <- list(k = k, lambda = lambda)
    return(object)
}

#' @rdname runINMF
#' @export
#' @param datasetVar Metadata variable name that stores the dataset source
#' annotation. Default \code{"orig.ident"}.
#' @param layer For Seurat>=4.9.9, the name of layer to retrieve input
#' non-negative scaled data. Default \code{"ligerScaleData"}. For older Seurat,
#' always retrieve from \code{scale.data} slot.
#' @param assay Name of assay to use. Default \code{NULL} uses current active
#' assay.
#' @param reduction Name of the reduction to store result. Also used as the
#' feature key. Default \code{"inmf"}.
#' @method runINMF Seurat
runINMF.Seurat <- function(
        object,
        k = 20,
        lambda = 5.0,
        datasetVar = "orig.ident",
        layer = "ligerScaleData",
        assay = NULL,
        reduction = "inmf",
        nIteration = 30,
        nRandomStarts = 1,
        HInit = NULL,
        WInit = NULL,
        VInit = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    Es <- .getSeuratData(object, layer = layer, slot = "scale.data",
                              assay = assay)
    # the last [,1] converts data.frame to the vector/factor
    datasetVar <- object[[datasetVar]][,1]
    if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
    datasetVar <- droplevels(datasetVar)

    if (!is.list(Es)) {
        Es <- splitRmMiss(Es, datasetVar)
        Es <- lapply(Es, methods::as, Class = "CsparseMatrix")
    }
    for (i in seq_along(Es)) {
        if (any(Es[[i]]@x < 0)) {
            stop("Negative data encountered for integrative Non-negative ",
                 "Matrix Factorization. Please run `scaleNotCenter()` first.")
        }
    }

    res <- .runINMF.list(
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
    Hconcat <- t(Reduce(cbind, res$H))
    colnames(Hconcat) <- paste0(reduction, "_", seq_len(ncol(Hconcat)))
    object[[reduction]] <- Seurat::CreateDimReducObject(
        embeddings = Hconcat,
        loadings = res$W,
        assay = assay,
        misc = list(V = res$V, objErr = res$objErr, dataset = datasetVar)
    )
    return(object)
}

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
#' @noRd
.runINMF.list <- function(
        object,
        k = 20,
        lambda = 5.0,
        nIteration = 30,
        nRandomStarts = 1,
        HInit = NULL,
        WInit = NULL,
        VInit = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        barcodeList = NULL,
        features = NULL
) {
    if (!requireNamespace("RcppPlanc", quietly = TRUE)) # nocov start
        stop("RcppPlanc installation required. Currently, please get the ",
             "GitHub private repository access from the lab and run: \n",
             "devtools::install_github(\"welch-lab/RcppPlanc\")") # nocov end

    bestResult <- list()
    bestObj <- Inf
    bestSeed <- seed
    for (i in seq(nRandomStarts)) {
        if (isTRUE(verbose) && nRandomStarts > 1) {
            .log("Replicate run ", i, "...")
        }
        set.seed(seed = seed + i - 1)
        if (inherits(object[[1]], "H5D")) {
            # RcppPlanc::bppinmf_h5dense()
            stop("TODO: Push Yichen to test bppinmf_h5sparse/bppinmf_h5dense!")
        } else {
            out <- RcppPlanc::inmf(objectList = object, k = k, lambda = lambda,
                                   niter = nIteration, Hinit = HInit,
                                   Vinit = VInit, Winit = WInit,
                                   verbose = verbose)
        }
        if (out$objErr < bestObj) {
            bestResult <- out
            bestObj <- out$objErr
            bestSeed <- seed + i - 1
        }
    }
    if (isTRUE(verbose) && nRandomStarts > 1) {
        .log("Best objective error: ", bestObj, "\nBest seed: ", bestSeed)
    }
    barcodeList <- lapply(object, colnames)
    features <- rownames(object[[1]])
    factorNames <- paste0("Factor_", seq(k))
    for (i in seq_along(object)) {
        bestResult$H[[i]] <- t(bestResult$H[[i]])
        dimnames(bestResult$H[[i]]) <- list(factorNames, barcodeList[[i]])
        dimnames(bestResult$V[[i]]) <- list(features, factorNames)
    }
    names(bestResult$V) <- names(bestResult$H) <- names(object)
    dimnames(bestResult$W) <- list(features, factorNames)
    return(bestResult)
}

#' [Deprecated] Perform iNMF on scaled datasets
#' @description
#' \bold{Please turn to \code{\link{runINMF}} or \code{\link{runIntegration}}}.
#'
#' Perform integrative non-negative matrix factorization to return factorized H,
#' W, and V matrices. It optimizes the iNMF objective function using block
#' coordinate descent (alternating non-negative least squares), where the number
#' of factors is set by k. TODO: include objective function equation here in
#' documentation (using deqn)
#'
#' For each dataset, this factorization produces an H matrix (cells by k), a V
#' matrix (k by genes), and a shared W matrix (k by genes). The H matrices
#' represent the cell factor loadings. W is held consistent among all datasets,
#' as it represents the shared components of the metagenes across datasets. The
#' V matrices represent the dataset-specific components of the metagenes.
#' @param object \code{liger} object. Should normalize, select genes, and scale
#' before calling.
#' @param k Inner dimension of factorization (number of factors). Run suggestK
#' to determine appropriate value; a general rule of thumb is that a higher k
#' will be needed for datasets with more sub-structure.
#' @param lambda Regularization parameter. Larger values penalize
#' dataset-specific effects more strongly (ie. alignment should increase as
#' lambda increases). Run suggestLambda to determine most appropriate value for
#' balancing dataset alignment and agreement (default 5.0).
#' @param thresh Convergence threshold. Convergence occurs when
#' |obj0-obj|/(mean(obj0,obj)) < thresh. (default 1e-6)
#' @param max.iters Maximum number of block coordinate descent iterations to
#' perform (default 30).
#' @param nrep Number of restarts to perform (iNMF objective function is
#' non-convex, so taking the best objective from multiple successive
#' initializations is recommended). For easier reproducibility, this increments
#' the random seed by 1 for each consecutive restart, so future factorizations
#' of the same dataset can be run with one rep if necessary. (default 1)
#' @param H.init Initial values to use for H matrices. (default NULL)
#' @param W.init Initial values to use for W matrix (default NULL)
#' @param V.init Initial values to use for V matrices (default NULL)
#' @param rand.seed Random seed to allow reproducible results (default 1).
#' @param print.obj Print objective function values after convergence (default
#' FALSE).
#' @param verbose Print progress bar/messages (TRUE by default)
#' @param ... Arguments passed to other methods
#' @return \code{liger} object with H, W, and V slots set.
#' @name optimizeALS-deprecated
#' @seealso \code{\link{rliger2-deprecated}}
NULL

#' @rdname rliger2-deprecated
#' @section \code{optimizeALS}:
#' For \code{optimizeALS}, use \code{\link{runIntegration}} or
#' \code{\link{runINMF}}. For the case of
#' \code{optimizeALS(use.unshared = TRUE)}, use \code{\link{runIntegration}}
#' with \code{method = "UINMF"} or \code{\link{runUINMF}} instead.
#' @export
optimizeALS <- function( # nocov start
        object,
        k,
        lambda = 5.0,
        thresh = NULL,
        max.iters = 30,
        nrep = 1,
        H.init = NULL,
        W.init = NULL,
        V.init = NULL,
        use.unshared = FALSE,
        rand.seed = 1,
        print.obj = NULL,
        verbose = TRUE,
        ...
) {
    if (isTRUE(use.unshared)) {
        lifecycle::deprecate_warn(
            "1.99.0", "optimizeALS(use.unshared = 'TRUE')",
            details = "Please use `runIntegration()` with `method = 'UINMF'` or `runUINMF()` instead."
        )
        # Call UINMF
        object <- runUINMF(object = object, k = k, lambda = lambda,
                           nIteration = max.iters, nRandomStarts = nrep,
                           seed = rand.seed, verbose = verbose)
    } else {
        lifecycle::deprecate_warn(
            "1.99.0", "optimizeALS()",
            "runIntegration()")
        object <- runINMF(object = object, k = k, lambda = lambda,
                          nIteration = max.iters, nRandomStarts = nrep,
                          HInit = H.init, WInit = W.init, VInit = V.init,
                          seed = rand.seed, verbose = verbose)
    }
    return(object)
} # nocov end

############################### online inmf ####################################

#' Perform online iNMF on scaled datasets
#' @description Perform online integrative non-negative matrix factorization to
#' represent multiple single-cell datasets in terms of \eqn{H}, \eqn{W}, and
#' \eqn{V} matrices. It optimizes the iNMF objective function (see
#' \code{\link{runINMF}}) using online learning (non-negative least squares for
#' \eqn{H} matrices, and hierarchical alternating least squares (HALS) for
#' \eqn{V} matrices and \eqn{W}), where the number of factors is set by
#' \code{k}. The function allows online learning in 3 scenarios:
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
#' a \eqn{V} matrix (genes by k), and a shared \eqn{W}
#' matrix (genes by k). The \eqn{H} matrices represent the cell factor loadings.
#' \eqn{W} is identical among all datasets, as it represents the shared
#' components of the metagenes across datasets. The \eqn{V} matrices represent
#' the dataset-specific components of the metagenes.
#'
#' @details
#' For performing scenario 2 or 3, a complete set of factorization result from
#' a run of scenario 1 is required. Given the structure of a \linkS4class{liger}
#' object, all of the required information can be retrieved automatically.
#' Under the circumstance where users need customized information for existing
#' factorization, arguments \code{WInit}, \code{VInit}, \code{AInit} and
#' \code{BInit} are exposed. The requirements for these argument follows:
#' \itemize{
#'  \item{WInit - A matrix object of size \eqn{m \times k}. (see
#'      \code{\link{runINMF}} for notation)}
#'  \item{VInit - A list object of matrices each of size \eqn{m \times k}.
#'      Number of matrices should match with \code{newDatasets}.}
#'  \item{AInit - A list object of matrices each of size \eqn{k \times k}.
#'      Number of matrices should match with \code{newDatasets}.}
#'  \item{BInit - A list object of matrices each of size \eqn{m \times k}.
#'      Number of matrices should match with \code{newDatasets}.}
#' }
#'
#' Minibatch iterations is performed on small subset of cells. The exact
#' minibatch size applied on each dataset is \code{minibatchSize} multiplied by
#' the proportion of cells in this dataset out of all cells. In general,
#' \code{minibatchSize} should be no larger than the number of cells in the
#' smallest dataset (considering both \code{object} and \code{newDatasets}).
#' Therefore, a smaller value may be necessary for analyzing very small
#' datasets.
#'
#' An epoch is one completion of calculation on all cells after a number of
#' iterations of minibatches. Therefore, the total number of iterations is
#' determined by the setting of \code{maxEpochs}, total number of cells, and
#' \code{minibatchSize}.
#'
#' Currently, Seurat S3 method does not support working on Scenario 2 and 3,
#' because there is no simple solution for organizing a number of miscellaneous
#' matrices with a single Seurat object. We strongly recommend that users create
#' a \linkS4class{liger} object which has the specific structure.
#' @param object \linkS4class{liger} object. Scaled data required.
#' @param newDatasets Named list of \linkS4class{dgCMatrix}. New datasets for
#' scenario 2 or scenario 3. Default \code{NULL} triggers scenario 1.
#' @param projection Whether to perform data integration with scenario 3 when
#' \code{newDatasets} is specified. See description. Default \code{FALSE}.
#' @param WInit,VInit,AInit,BInit Optional initialization for \eqn{W}, \eqn{V},
#' \eqn{A}, and \eqn{B} matrices, respectively. Must be presented all together.
#' See detail. Default \code{NULL}.
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
#' @param HALSiter Maximum number of block coordinate descent (HALS
#' algorithm) iterations to perform for each update of \eqn{W} and \eqn{V}.
#' Default \code{1}. Changing this parameter is not recommended.
#' @param minibatchSize Total number of cells in each minibatch. See detail.
#' Default \code{5000}.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param ... Arguments passed to other S3 methods of this function.
#' @return
#' \itemize{
#'  \item{liger method - Returns updated input \linkS4class{liger} object.
#'  \itemize{
#'      \item{A list of all \eqn{H} matrices can be accessed with
#'          \code{getMatrix(object, "H")}}
#'      \item{A list of all \eqn{V} matrices can be accessed with
#'          \code{getMatrix(object, "V")}}
#'      \item{The \eqn{W} matrix can be accessed with
#'          \code{getMatrix(object, "W")}}
#'      \item{Meanwhile, intermediate matrices \eqn{A} and \eqn{B} produced in
#'          HALS update can also be accessed similarly.}
#'  }
#'  }
#'  \item{Seurat method - Returns updated input Seurat object.
#'  \itemize{
#'      \item{\eqn{H} matrices for all datasets will be concatenated and
#'          transposed (all cells by k), and form a DimReduc object in the
#'          \code{reductions} slot named by argument \code{reduction}.}
#'      \item{\eqn{W} matrix will be presented as \code{feature.loadings} in the
#'          same DimReduc object.}
#'      \item{\eqn{V} matrices, \eqn{A} matrices, \eqn{B} matricesm an objective
#'          error value and the dataset variable used for the factorization is
#'          currently stored in \code{misc} slot of the same DimReduc object.}
#'  }}
#' }
#' @references Chao Gao and et al., Iterative single-cell multi-omic integration
#' using online learning, Nat Biotechnol., 2021
#' @export
#' @rdname runOnlineINMF
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' if (requireNamespace("RcppPlanc", quietly = TRUE)) {
#'     # Scenario 1
#'     pbmc <- runOnlineINMF(pbmc, minibatchSize = 200)
#'     # Scenario 2
#'     # Fake new dataset by increasing all non-zero value in "ctrl" by 1
#'     ctrl2 <- rawData(dataset(pbmc, "ctrl"))
#'     ctrl2@x <- ctrl2@x + 1
#'     colnames(ctrl2) <- paste0(colnames(ctrl2), 2)
#'     pbmc2 <- runOnlineINMF(pbmc, k = 20, newDatasets = list(ctrl2 = ctrl2),
#'                            minibatchSize = 100)
#'     # Scenario 3
#'     pbmc3 <- runOnlineINMF(pbmc, k = 20, newDatasets = list(ctrl2 = ctrl2),
#'                            projection = TRUE)
#' }
runOnlineINMF <- function(
        object,
        k = 20,
        lambda = 5,
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
        minibatchSize = 5000,
        WInit = NULL,
        VInit = NULL,
        AInit = NULL,
        BInit = NULL,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    .checkObjVersion(object)
    object <- recordCommand(object, ..., dependencies = c("RcppPlanc"))
    Es <- getMatrix(object, "scaleData", returnList = TRUE)
    Es <- lapply(datasets(object), function(ld) {
        sd <- scaleData(ld)
        if (is.null(sd))
            stop("Scaled data not available. Run `scaleNotCenter()` first")
        # if (inherits(sd, "H5D")) return(.H5DToH5Mat(sd))
        # else
        if (inherits(sd, "H5Group"))
            return(.H5GroupToH5SpMat(sd, c(length(varFeatures(object)),
                                           ncol(ld))))
        else return(sd)
    })
    if (!is.null(newDatasets)) {
        WInit <- WInit %||% getMatrix(object, "W", returnList = FALSE)
        VInit <- VInit %||% getMatrix(object, "V", returnList = TRUE)
        AInit <- AInit %||% getMatrix(object, "A", returnList = TRUE)
        BInit <- BInit %||% getMatrix(object, "B", returnList = TRUE)
        if (is.null(WInit) || any(sapply(VInit, is.null)) ||
            any(sapply(AInit, is.null)) || any(sapply(BInit, is.null))) {
            stop("Cannot find complete online iNMF result for current ",
                 "datasets. Please run `runOnlineINMF()` without `newDataset` ",
                 "first.")
        }

        newNames <- names(newDatasets)
        if (any(newNames %in% names(object))) {
            stop("Names of `newDatasets` overlap with existing datasets.")
        }
        if (is.list(newDatasets)) {
            # A list of raw data
            if (is.null(names(newDatasets))) {
                stop("The list of new datasets must be named.")
            }
            for (i in seq_along(newDatasets)) {
                if (inherits(newDatasets[[i]], "dgCMatrix")) {
                    dataset(object, names(newDatasets)[i]) <- newDatasets[[i]]
                } else if (is.character(newDatasets[[i]])) {
                    # Assume is H5 filename
                    ld <- createH5LigerDataset(newDatasets[[i]])
                    dataset(object, names(newDatasets[i])) <- ld
                } else {
                    stop("Cannot interpret `newDatasets` element ", i)
                }
            }
        } else if (inherits(newDatasets, "liger")) {
            # A liger object with all new datasets
            object <- c(object, newDatasets)
        } else {
            stop("`newDatasets` must be either a named list or a liger object")
        }

        object <- normalize(object, useDatasets = newNames)
        object <- scaleNotCenter(object, useDatasets = newNames)
        newDatasets <- list()
        for (d in newNames) {
            ld <- dataset(object, d)
            sd <- scaleData(ld)
            # if (inherits(sd, "H5D")) {
            #     newDatasets[[d]] <- .H5DToH5Mat(sd)
            # } else
            if (inherits(sd, "H5Group")) {
                newDatasets[[d]] <- .H5GroupToH5SpMat(sd, c(length(varFeatures(object)), ncol(ld)))
            } else {
                newDatasets[[d]] <- sd
            }
        }
    }
    object <- closeAllH5(object)
    res <- .runOnlineINMF.list(Es, newDatasets = newDatasets,
                               projection = projection, k = k, lambda = lambda,
                               maxEpochs = maxEpochs,
                               minibatchSize = minibatchSize,
                               HALSiter = HALSiter, verbose = verbose,
                               WInit = WInit, VInit = VInit, AInit = AInit,
                               BInit = BInit, seed = seed)
    if (!isTRUE(projection)) {
        # Scenario 1&2, everything updated
        for (i in seq_along(object)) {
            ld <- dataset(object, i)
            ld@H <- res$H[[i]]
            ld@V <- res$V[[i]]
            ld@A <- res$A[[i]]
            ld@B <- res$B[[i]]
            if (isH5Liger(ld)) {
                colnames(ld@H) <- colnames(ld)
                rownames(ld@V) <- varFeatures(object)
                rownames(ld@B) <- varFeatures(object)
            }
            datasets(object, check = FALSE)[[i]] <- ld
        }
        object@W <- res$W
        rownames(object@W) <- varFeatures(object)
    } else {
        # Scenario 3, only H of newDatasets returned
        for (i in seq_along(newDatasets)) {
            dname <- names(newDatasets)[i]
            ld <- dataset(object, dname)
            ld@H <- res$H[[i]]
            if (isH5Liger(ld)) {
                colnames(ld@H) <- colnames(ld)
            }
            datasets(object, check = FALSE)[[dname]] <- ld
        }
    }
    object@uns$factorization <- list(k = k, lambda = lambda)
    suppressMessages({object <- restoreH5Liger(object)})
    return(object)
}

.runOnlineINMF.list <- function(
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
        minibatchSize = 5000,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    if (!requireNamespace("RcppPlanc", quietly = TRUE)) # nocov start
        stop("RcppPlanc installation required. Currently, please get the ",
             "GitHub private repository access from the lab and run: \n",
             "devtools::install_github(\"welch-lab/RcppPlanc\")") # nocov end
    nDatasets <- length(object) + length(newDatasets)
    barcodeList <- c(lapply(object, colnames), lapply(newDatasets, colnames))
    names(barcodeList) <- c(names(object), names(newDatasets))
    features <- rownames(object[[1]])
    if (!is.null(seed)) set.seed(seed)

    # # If minibatchSize > smallest invovled dataset, auto reset with warning
    # minibatchSize_min <- min(sapply(object, ncol))
    # if (!is.null(newDatasets)) {
    #     minibatchSize_min2 <- min(sapply(newDatasets, ncol))
    #     minibatchSize_min <- min(minibatchSize_min, minibatchSize_min2)
    # }
    # if (minibatchSize > minibatchSize_min) {
    #     warning("Minibatch size larger than the smallest dataset involved.\n",
    #             "  Setting to the smallest dataset size: ", minibatchSize_min,
    #             immediate. = TRUE)
    #     minibatchSize <- minibatchSize_min
    # }
    res <- RcppPlanc::onlineINMF(objectList = object, newDatasets = newDatasets,
                                 project = projection, k = k, lambda = lambda,
                                 maxEpoch = maxEpochs,
                                 minibatchSize = minibatchSize,
                                 maxHALSIter = HALSiter, Vinit = VInit,
                                 Winit = WInit, Ainit = AInit, Binit = BInit,
                                 verbose = verbose)
    factorNames <- paste0("Factor_", seq(k))
    if (isTRUE(projection)) {
        # Scenario 3 only got H for new datasets
        for (i in seq_along(newDatasets)) {
            dname <- names(newDatasets)[i]
            res$H[[i]] <- t(res$H[[i]])
            dimnames(res$H[[i]]) <- list(factorNames, barcodeList[[dname]])
            names(res$H) <- names(newDatasets)
        }
    } else {
        # Scenario 1&2 got everything
        for (i in seq(nDatasets)) {
            res$H[[i]] <- t(res$H[[i]])
            dimnames(res$H[[i]]) <- list(factorNames, barcodeList[[i]])
            dimnames(res$V[[i]]) <- list(features, factorNames)
            dimnames(res$A[[i]]) <- list(factorNames, factorNames)
            dimnames(res$B[[i]]) <- list(features, factorNames)
        }
        names(res$B) <- names(res$A) <- names(res$V) <- names(res$H) <-
            names(barcodeList)
        dimnames(res$W) <- list(features, factorNames)
    }
    return(res)
}

#' @export
#' @rdname runOnlineINMF
#' @method runOnlineINMF Seurat
#' @param datasetVar Metadata variable name that stores the dataset source
#' annotation. Default \code{"orig.ident"}.
#' @param layer For Seurat>=4.9.9, the name of layer to retrieve input
#' non-negative scaled data. Default \code{"ligerScaleData"}. For older Seurat,
#' always retrieve from \code{scale.data} slot.
#' @param assay Name of assay to use. Default \code{NULL} uses current active
#' assay.
#' @param reduction Name of the reduction to store result. Also used as the
#' feature key. Default \code{"onlineINMF"}.
runOnlineINMF.Seurat <- function(
        object,
        k = 20,
        lambda = 5,
        datasetVar = "orig.ident",
        layer = "ligerScaleData",
        assay = NULL,
        reduction = "onlineINMF",
        maxEpochs = 5,
        HALSiter = 1,
        minibatchSize = 5000,
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    Es <- .getSeuratData(object, layer = layer, slot = "scale.data",
                         assay = assay)
    # the last [,1] converts data.frame to the vector/factor
    datasetVar <- object[[datasetVar]][,1]
    if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
    datasetVar <- droplevels(datasetVar)

    if (!is.list(Es)) {
        Es <- splitRmMiss(Es, datasetVar)
        Es <- lapply(Es, methods::as, Class = "CsparseMatrix")
    }
    for (i in seq_along(Es)) {
        if (any(Es[[i]]@x < 0)) {
            stop("Negative data encountered for integrative Non-negative ",
                 "Matrix Factorization. Please run `scaleNotCenter()` first.")
        }
    }

    res <- .runOnlineINMF.list(
        object = Es, k = k, lambda = lambda,
        newDatasets = NULL, projection = FALSE,
        maxEpochs = maxEpochs, HALSiter = HALSiter,
        minibatchSize = minibatchSize, seed = seed, verbose = verbose,
        WInit = NULL, VInit = NULL, AInit = NULL, BInit = NULL,
    )
    Hconcat <- t(Reduce(cbind, res$H))
    colnames(Hconcat) <- paste0(reduction, "_", seq_len(ncol(Hconcat)))
    object[[reduction]] <- Seurat::CreateDimReducObject(
        embeddings = Hconcat,
        loadings = res$W,
        assay = assay,
        misc = list(V = res$V, A = res$A, B = res$B, objErr = res$objErr,
                    dataset = datasetVar)
    )
    return(object)
}

#' [Deprecated] Perform online iNMF on scaled datasets
#' @description
#' \bold{Please turn to \code{\link{runOnlineINMF}} or
#' \code{\link{runIntegration}}}.
#'
#' Perform online integrative non-negative matrix factorization to represent
#' multiple single-cell datasets in terms of H, W, and V matrices. It optimizes
#' the iNMF objective function using online learning (non-negative least squares
#' for H matrix, hierarchical alternating least squares for W and V matrices),
#' where the number of factors is set by k. The function allows online learning
#' in 3 scenarios: (1) fully observed datasets; (2) iterative refinement using
#' continually arriving datasets; and (3) projection of new datasets without
#' updating the existing factorization. All three scenarios require fixed memory
#' independent of the number of cells.
#'
#' For each dataset, this factorization produces an H matrix (cells by k), a V
#' matrix (k by genes), and a shared W matrix (k by genes). The H matrices
#' represent the cell factor loadings. W is identical among all datasets, as it
#' represents the shared components of the metagenes across datasets. The V
#' matrices represent the dataset-specific components of the metagenes.
#' @param object \code{liger} object with data stored in HDF5 files. Should
#' normalize, select genes, and scale before calling.
#' @param X_new List of new datasets for scenario 2 or scenario 3. Each list
#' element should be the name of an HDF5 file.
#' @param projection Perform data integration by shared metagene (W) projection
#' (scenario 3). (default FALSE)
#' @param W.init Optional initialization for W. (default NULL)
#' @param V.init Optional initialization for V (default NULL)
#' @param H.init Optional initialization for H (default NULL)
#' @param A.init Optional initialization for A (default NULL)
#' @param B.init Optional initialization for B (default NULL)
#' @param k Inner dimension of factorization--number of metagenes (default 20).
#' A value in the range 20-50 works well for most analyses.
#' @param lambda Regularization parameter. Larger values penalize
#' dataset-specific effects more
#'   strongly (ie. alignment should increase as lambda increases). We recommend
#'   always using the default value except
#'   possibly for analyses with relatively small differences (biological
#'   replicates, male/female comparisons, etc.)
#'   in which case a lower value such as 1.0 may improve reconstruction quality.
#'   (default 5.0).
#' @param max.epochs Maximum number of epochs (complete passes through the
#' data). (default 5)
#' @param miniBatch_max_iters Maximum number of block coordinate descent (HALS
#' algorithm) iterations to perform for each update of W and V (default 1).
#' Changing this parameter is not  recommended.
#' @param miniBatch_size Total number of cells in each minibatch (default 5000).
#' This is a reasonable default, but a smaller value such as 1000 may be
#' necessary for analyzing very small datasets. In general, minibatch size
#' should be no larger than the number of cells in the smallest dataset.
#' @param h5_chunk_size Chunk size of input hdf5 files (default 1000). The chunk
#' size should be no larger than the batch size.
#' @param seed Random seed to allow reproducible results (default 123).
#' @param verbose Print progress bar/messages (TRUE by default)
#' @return \code{liger} object with H, W, V, A and B slots set.
#' @name online_iNMF-deprecated
NULL

#' @rdname rliger2-deprecated
#' @section \code{online_iNMF}:
#' For \code{online_iNMF}, use \code{\link{runIntegration}} with
#' \code{method = "online"} or \code{\link{runOnlineINMF}}.
#' @export
online_iNMF <- function( # nocov start
        object,
        X_new = NULL,
        projection = FALSE,
        W.init = NULL,
        V.init = NULL,
        H.init = NULL,
        A.init = NULL,
        B.init = NULL,
        k = 20,
        lambda = 5,
        max.epochs = 5,
        miniBatch_max_iters = 1,
        miniBatch_size = 5000,
        h5_chunk_size = 1000,
        seed = 123,
        verbose = TRUE
) {
    lifecycle::deprecate_warn(
        "1.99.0", "online_iNMF()",
        details = "Please use `runIntegration()` with `method = 'online'`, or `runOnlineINMF()` instead."
    )
    object <- runOnlineINMF.liger(
        object = object, k = k, lambda = lambda, maxEpochs = max.epochs,
        HALSiter = miniBatch_max_iters, minibatchSize = miniBatch_size,
        seed = seed, verbose = verbose, newDatasets = X_new,
        projection = projection, WInit = W.init, VInit = V.init, AInit = A.init,
        BInit = B.init
    )
    return(object)
} # nocov end



################################### UINMF ######################################

#' Perform Mosaic iNMF (UINMF) on scaled datasets with unshared features
#' @description
#' Performs mosaic integrative non-negative matrix factorization (UINMF) (A.R.
#' Kriebel, 2022) using block coordinate descent (alternating non-negative
#' least squares, ANLS) to return factorized \eqn{H}, \eqn{W}, \eqn{V} and
#' \eqn{U} matrices. The objective function is stated as
#'
#' \deqn{\arg\min_{H\ge0,W\ge0,V\ge0,U\ge0}\sum_{i}^{d}
#' ||\begin{bmatrix}E_i \\ P_i \end{bmatrix} -
#' (\begin{bmatrix}W \\ 0 \end{bmatrix}+
#' \begin{bmatrix}V_i \\ U_i \end{bmatrix})Hi||^2_F+
#' \lambda_i\sum_{i}^{d}||\begin{bmatrix}V_i \\ U_i \end{bmatrix}H_i||_F^2}
#'
#' where \eqn{E_i} is the input non-negative matrix of the \eqn{i}'th dataset,
#' \eqn{P_i} is the input non-negative matrix for the unshared features,
#' \eqn{d} is the total number of datasets. \eqn{E_i} is of size
#' \eqn{m \times n_i} for \eqn{m} shared features and \eqn{n_i} cells, \eqn{P_i}
#' is of size \eqn{u_i \times n_i} for \eqn{u_i} unshared feaetures,
#' \eqn{H_i} is of size \eqn{k \times n_i}, \eqn{V_i} is of size
#' \eqn{m \times k}, \eqn{W} is of size \eqn{m \times k} and \eqn{U_i} is of
#' size \eqn{u_i \times k}.
#'
#' The factorization produces a shared \eqn{W} matrix (genes by k). For each
#' dataset, an \eqn{H} matrix (k by cells), a \eqn{V} matrix (genes by k) and
#' a \eqn{U} matrix (unshared genes by k). The \eqn{H} matrices represent the
#' cell factor loadings. \eqn{W} is held consistent among all datasets, as it
#' represents the shared components of the metagenes across datasets. The
#' \eqn{V} matrices represent the dataset-specific components of the metagenes,
#' \eqn{U} matrices are similar to \eqn{V}s but represents the loading
#' contributed by unshared features.
#'
#' This function adopts highly optimized fast and memory efficient
#' implementation extended from Planc (Kannan, 2016). Pre-installation of
#' extension package \code{RcppPlanc} is required. The underlying algorithm
#' adopts the identical ANLS strategy as \code{\link{optimizeALS}(unshared =
#' TRUE)} in the old version of LIGER.
#' @param object \linkS4class{liger} object. Should run
#' \code{\link{selectGenes}} with \code{unshared = TRUE} and then run
#' \code{\link{scaleNotCenter}} in advance.
#' @param k Inner dimension of factorization (number of factors). Generally, a
#' higher \code{k} will be needed for datasets with more sub-structure. Default
#' \code{20}.
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
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param ... Arguments passed to other methods and wrapped functions.
#' @export
#' @references April R. Kriebel and Joshua D. Welch, UINMF performs mosaic
#' integration of single-cell multi-omic datasets using nonnegative matrix
#' factorization, Nat. Comm., 2022
#' @note
#' Currently, Seurat S3 method is not supported for UINMF because there is no
#' simple solution for organizing a number of miscellaneous matrices with a
#' single Seurat object. We strongly recommend that users create a
#' \linkS4class{liger} object which has the specific structure.
#' @return
#' \itemize{
#'  \item{liger method - Returns updated input \linkS4class{liger} object.
#'  \itemize{
#'      \item{A list of all \eqn{H} matrices can be accessed with
#'          \code{getMatrix(object, "H")}}
#'      \item{A list of all \eqn{V} matrices can be accessed with
#'          \code{getMatrix(object, "V")}}
#'      \item{The \eqn{W} matrix can be accessed with
#'          \code{getMatrix(object, "W")}}
#'      \item{A list of all \eqn{U} matrices can be accessed with
#'          \code{getMatrix(object, "U")}}
#'  }
#'  }
#' }
#' @rdname runUINMF
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc, useUnsharedDatasets = c("ctrl", "stim"))
#' pbmc <- scaleNotCenter(pbmc)
#' if (!is.null(getMatrix(pbmc, "scaleUnsharedData", "ctrl")) &&
#'     !is.null(getMatrix(pbmc, "scaleUnsharedData", "stim"))) {
#'     # TODO: unshared variable features cannot be detected from this example
#'     pbmc <- runUINMF(pbmc)
#' }
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
    object <- recordCommand(object, ..., dependencies = "RcppPlanc")
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
             "with `useUnsharedDatasets` specified, ",
             "and then `scaleNotCenter()`.")
    }
    res <- .runUINMF.list(Elist, Ulist, k = k, lambda = lambda,
                          nIteration = nIteration,
                          nRandomStarts = nRandomStarts,
                          seed = seed, verbose = verbose, ...)
    for (d in names(object)) {
        ld <- dataset(object, d)
        ld@H <- res$H[[d]]
        ld@V <- res$V[[d]]
        if (!is.null(ld@scaleUnsharedData)) {
            ld@U <- res$U[[d]]
        }
        datasets(object, check = FALSE)[[d]] <- ld
    }
    object@W <- res$W
    object@uns$factorization <- list(k = k, lambda = lambda)
    return(object)
}

#' @param unsharedList List of matrices for unshared features
#' @noRd
.runUINMF.list <- function(
        object,
        unsharedList,
        k = 20,
        lambda = 5,
        nIteration = 30,
        nRandomStarts = 1,
        seed = 1,
        verbose = getOption("ligerVerbose")
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
    for (d in names(object)) {
        bestRes$H[[d]] <- t(bestRes$H[[d]])
        dimnames(bestRes$H[[d]]) <- list(factorNames, barcodes[[d]])
        dimnames(bestRes$V[[d]]) <- list(features, factorNames)
        if (d %in% names(bestRes$U)) {
            dimnames(bestRes$U[[d]]) <- list(unsharedFeatures[[d]], factorNames)
        }
    }
    dimnames(bestRes$W) <- list(features, factorNames)
    return(bestRes)
}





########################### Quantile Normalization #############################

#' Quantile Align (Normalize) Factor Loadings
#' @description This process builds a shared factor neighborhood graph to
#' jointly cluster cells, then quantile normalizes corresponding clusters.
#'
#' The first step, building the shared factor neighborhood graph, is performed
#' in SNF(), and produces a graph representation where edge weights between
#' cells (across all datasets) correspond to their similarity in the shared
#' factor neighborhood space. An important parameter here is \code{nNeighbors},
#' the number of neighbors used to build the shared factor space.
#'
#' Next we perform quantile alignment for each dataset, factor, and cluster (by
#' stretching/compressing datasets' quantiles to better match those of the
#' reference dataset).
#' @param object A \linkS4class{liger} or Seurat object with valid factorization
#' result available (i.e. \code{\link{runIntegration}} performed in advance).
#' @param quantiles Number of quantiles to use for quantile normalization.
#' Default \code{50}.
#' @param reference Character, numeric or logical selection of one dataset, out
#' of all available datasets in \code{object}, to use as a "reference" for
#' normalization. Default \code{NULL} use the dataset with the largest number of
#' cells.
#' @param minCells Minimum number of cells to consider a cluster shared across
#' datasets. Default \code{20}.
#' @param nNeighbors Number of nearest neighbors for within-dataset knn graph.
#' Default \code{20}.
#' @param useDims Indices of factors to use for shared nearest factor
#' determination. Default \code{NULL} uses all factors.
#' @param center Whether to center the data when scaling factors. Could be
#' useful for less sparse modalities like methylation data. Default
#' \code{FALSE}.
#' @param maxSample Maximum number of cells used for quantile normalization of
#' each cluster and factor. Default \code{1000}.
#' @param eps The error bound of the nearest neighbor search. Lower values give
#' more accurate nearest neighbor graphs but take much longer to compute.
#' Default \code{0.9}.
#' @param refineKNN whether to increase robustness of cluster assignments using
#' KNN graph. Default \code{TRUE}.
#' @param clusterName Variable name that will store the clustering result
#' in metadata of a \linkS4class{liger} object or a \code{Seurat} object.
#' Default \code{"quantileNorm_cluster"}
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param ... Arguments passed to other S3 methods of this function.
#' @return Updated input object
#' \itemize{
#'  \item{liger method
#'  \itemize{
#'      \item{Update the \code{H.norm} slot for the alignment cell factor
#'          loading, ready for running graph based community detection
#'          clustering or dimensionality reduction for visualization.}
#'      \item{Update the \code{cellMata} slot with a cluster assignment basing
#'          on cell factor loading}
#'  }}
#'  \item{Seurat method
#'  \itemize{
#'      \item{Update the \code{reductions} slot with a new \code{DimReduc}
#'          object containing the aligned cell factor loading.}
#'      \item{Update the metadata with a cluster assignment basing on cell
#'          factor loading}
#'  }}
#' }
#' @examples
#' pbmc <- quantileNorm(pbmcPlot)
#' @export
#' @rdname quantileNorm
quantileNorm <- function(
        object,
        ...
) {
    UseMethod("quantileNorm", object)
}

#' @export
#' @rdname quantileNorm
#' @method quantileNorm liger
quantileNorm.liger <- function(
        object,
        quantiles = 50,
        reference = NULL,
        minCells = 20,
        nNeighbors = 20,
        useDims = NULL,
        center = FALSE,
        maxSample = 1000,
        eps = 0.9,
        refineKNN = TRUE,
        clusterName = "quantileNorm_cluster",
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    .checkObjVersion(object)
    .checkValidFactorResult(object, checkV = FALSE)
    reference <- reference %||% names(which.max(sapply(datasets(object), ncol)))
    reference <- .checkUseDatasets(object, useDatasets = reference)
    if (length(reference) != 1)
        stop("Should specify only one reference dataset.")
    object <- recordCommand(object, ..., dependencies = "RANN")
    out <- .quantileNorm.HList(
        object = getMatrix(object, "H"),
        quantiles = quantiles,
        reference = reference,
        minCells = minCells,
        nNeighbors = nNeighbors,
        useDims = useDims,
        center = center,
        maxSample = maxSample,
        eps = eps,
        refineKNN = refineKNN,
        seed = seed
    )
    object@H.norm <- out$H.norm
    cellMeta(object, clusterName, check = FALSE) <- out$clusters
    return(object)
}

#' @export
#' @rdname quantileNorm
#' @method quantileNorm Seurat
#' @param reduction Name of the reduction where LIGER integration result is
#' stored. Default \code{"inmf"}.
quantileNorm.Seurat <- function(
        object,
        reduction = "inmf",
        quantiles = 50,
        reference = NULL,
        minCells = 20,
        nNeighbors = 20,
        useDims = NULL,
        center = FALSE,
        maxSample = 1000,
        eps = 0.9,
        refineKNN = TRUE,
        clusterName = "quantileNorm_cluster",
        seed = 1,
        verbose = getOption("ligerVerbose"),
        ...
) {
    resName <- paste0(reduction, "Norm")
    reduction <- object[[reduction]]
    if (!inherits(reduction, "DimReduc")) {
        stop("Specified `reduction` does not points to a DimReduc.")
    }
    # Retrieve some information. Might have better ways instead of using `@`
    ## Due to proper formatting in Seurat object, Hconcat is already cell x k
    Hconcat <- reduction[[]]
    datasetVar <- reduction@misc$dataset
    W <- reduction[]
    assay <- reduction@assay.used

    reference <- reference %||% names(which.max(table(datasetVar)))
    # Debating whether to implement something without the need to split it
    HList <- lapply(levels(datasetVar), function(d) {
        t(Hconcat[datasetVar == d, , drop = FALSE])
    })
    names(HList) <- levels(datasetVar)
    alignment <- .quantileNorm.HList(
        HList, quantiles = quantiles, reference = reference,
        minCells = minCells, nNeighbors = nNeighbors, useDims = useDims,
        center = center, maxSample = maxSample, eps = eps,
        refineKNN = refineKNN, seed = seed, verbose = verbose
    )
    reddim <- Seurat::CreateDimReducObject(
        embeddings = alignment$H.norm, loadings = W,
        assay = assay, key = paste0(resName, "_")
    )
    object[[resName]] <- reddim
    object[[paste0(resName, ".cluster")]] <- alignment$clusters
    return(object)
}

.quantileNorm.HList <- function(
        object,
        quantiles = 50,
        reference = NULL,
        minCells = 20,
        nNeighbors = 20,
        useDims = NULL,
        center = FALSE,
        maxSample = 1000,
        eps = 0.9,
        refineKNN = TRUE,
        seed = 1,
        verbose = getOption("ligerVerbose")
) {
    set.seed(seed)
    if (is.character(reference)) {
        if (length(reference) != 1 || !reference %in% names(object))
            stop("Should specify one existing dataset as reference. ",
                 "(character `reference` wrong length or not found)")
    } else if (is.numeric(reference)) {
        if (length(reference) != 1 || reference > length(object))
            stop("Should specify one dataset within the range. ",
                 "(numeric `reference` wrong length or out of bound)")
    } else if (is.logical(reference)) {
        if (length(reference) != length(object) || sum(reference) != 1)
            stop("Should specify one dataset within the range. ",
                 "(logical `reference` wrong length or ",
                 "too many selection)")
    } else {
        stop("Unable to understand `reference`. See `?quantileNorm`.")
    }
    useDims <- useDims %||% seq_len(nrow(object[[1]]))
    # Transposing all H to cell x k
    Hs <- lapply(object, t)
    # fast max factor assignment with Rcpp code
    clusters <- lapply(Hs, max_factor_rcpp, dims_use = useDims, center = center)
    # TODO: Dumb to mess up factor with characters of numbers
    # Change to numeric in the future. Need to reproduce result for now.
    # clusterAssign <- factor(unlist(clusters))
    clusterAssign <- as.factor(unlist(lapply(clusters, as.character)))
    names(clusterAssign) <- unlist(lapply(Hs, rownames))
    # increase robustness of cluster assignments using knn graph
    if (isTRUE(refineKNN)) {
        for (i in seq_along(Hs)) {
            clustsH <- clusterAssign[rownames(Hs[[i]])]
            H_knn <- RANN::nn2(Hs[[i]], eps = eps, k = nNeighbors,
                               searchtype = "standard")
            # Rcpp method cluster_vote
            newClusts <- cluster_vote_rcpp(H_knn$nn.idx, clustsH)
            clusterAssign[rownames(Hs[[i]])] <- newClusts
        }
    }
    # TODO: This line need to be removed
    clusters <- lapply(Hs, function(H) clusterAssign[rownames(H)])
    dims <- ncol(Hs[[reference]])
    nClusters <- dims
    for (d in seq_along(Hs)) {
        for (c in seq(nClusters)) {
            # cells 2: cells in d'th dataset belonging to cluster c
            cellIdx2 <- clusters[[d]] == c
            # cells 1: cells in ref dataset belong to cluster c
            cellIdx1 <- clusters[[reference]] == c
            nCells2 <- sum(cellIdx2)
            nCells1 <- sum(cellIdx1)
            if (nCells1 < minCells || nCells2 < minCells) next
            for (k in seq(dims)) {
                if (nCells2 == 1) {
                    Hs[[d]][cellIdx2, k] <-
                        mean(Hs[[reference]][cellIdx1, k])
                    next
                }
                q2 <- stats::quantile(sample(Hs[[d]][cellIdx2, k],
                                             min(nCells2, maxSample)),
                                      seq(0, 1, by = 1 / quantiles))
                q1 <- stats::quantile(sample(Hs[[reference]][cellIdx1, k],
                                             min(nCells1, maxSample)),
                                      seq(0, 1, by = 1 / quantiles))
                if (sum(q1) == 0 | sum(q2) == 0 |
                    length(unique(q1)) < 2 | length(unique(q2)) < 2) {
                    newValue <- rep(0, nCells2)
                } else {
                    warp_func <- withCallingHandlers(
                        stats::approxfun(q2, q1, rule = 2),
                        warning = function(w) {
                            invokeRestart("muffleWarning")
                        }
                    )
                    newValue <- warp_func(Hs[[d]][cellIdx2, k])
                }
                Hs[[d]][cellIdx2, k] <- newValue
            }
        }
    }
    return(list('H.norm' = Reduce(rbind, Hs), 'clusters' = clusterAssign))
}

#' [Deprecated] Quantile align (normalize) factor loading
#' @description
#' \bold{Please turn to \code{\link{quantileNorm}}.}
#'
#' This process builds a shared factor neighborhood graph to jointly cluster
#' cells, then quantile normalizes corresponding clusters.
#'
#' The first step, building the shared factor neighborhood graph, is performed
#' in SNF(), and produces a graph representation where edge weights between
#' cells (across all datasets) correspond to their similarity in the shared
#' factor neighborhood space. An important parameter here is knn_k, the number
#' of neighbors used to build the shared factor space.
#'
#' Next we perform quantile alignment for each dataset, factor, and cluster (by
#' stretching/compressing datasets' quantiles to better match those of the
#' reference dataset). These aligned factor loadings are combined into a single
#' matrix and returned as H.norm.
#'
#' @param object \code{liger} object. Should run optimizeALS before calling.
#' @param knn_k Number of nearest neighbors for within-dataset knn graph
#' (default 20).
#' @param ref_dataset Name of dataset to use as a "reference" for normalization.
#' By default, the dataset with the largest number of cells is used.
#' @param min_cells Minimum number of cells to consider a cluster shared across
#' datasets (default 20)
#' @param quantiles Number of quantiles to use for quantile normalization
#' (default 50).
#' @param eps  The error bound of the nearest neighbor search. (default 0.9)
#' Lower values give more accurate nearest neighbor graphs but take much longer
#' to computer.
#' @param dims.use Indices of factors to use for shared nearest factor
#' determination (default 1:ncol(H[[1]])).
#' @param do.center Centers the data when scaling factors (useful for less
#' sparse modalities like methylation data). (default FALSE)
#' @param max_sample Maximum number of cells used for quantile normalization of
#' each cluster and factor. (default 1000)
#' @param refine.knn whether to increase robustness of cluster assignments using
#' KNN graph.(default TRUE)
#' @param rand.seed Random seed to allow reproducible results (default 1)
#' @return \code{liger} object with 'H.norm' and 'clusters' slot set.
#' @name quantile_norm-deprecated
#' @seealso \code{\link{rliger2-deprecated}}
NULL

#' @rdname rliger2-deprecated
#' @section \code{quantile_norm}:
#' For \code{quantile_norm}, use \code{\link{quantileNorm}}.
#' @export
quantile_norm <- function( # nocov start
        object,
        quantiles = 50,
        ref_dataset = NULL,
        min_cells = 20,
        knn_k = 20,
        dims.use = NULL,
        do.center = FALSE,
        max_sample = 1000,
        eps = 0.9,
        refine.knn = TRUE,
        clusterName = "H.norm_cluster",
        rand.seed = 1,
        verbose = getOption("ligerVerbose")
) {
    lifecycle::deprecate_warn("1.99.0", "quantile_norm()", "quantileNorm()")
    quantileNorm(
        object = object,
        quantiles = quantiles,
        reference = ref_dataset,
        minCells = min_cells,
        nNeighbors = knn_k,
        useDims = dims.use,
        center = do.center,
        maxSample = max_sample,
        eps = eps,
        refineKNN = refine.knn,
        seed = rand.seed,
        verbose = verbose
    )
} # nocov end
