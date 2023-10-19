selectGeneGlobalRank <- function(
        object,
        n = 4000,
        alpha = 0.99,
        useDatasets = NULL,
        unsharedDatasets = NULL,
        chunk = 1000,
        verbose = getOption("ligerVerbose")
) {
    .checkObjVersion(object)
    # A bunch of input checks at first ####
    useDatasets <- .checkUseDatasets(object, useDatasets)
    object <- recordCommand(object, dependencies = "hdf5r")
    if (!is.null(unsharedDatasets))
        unsharedDatasets <- .checkUseDatasets(object, unsharedDatasets)
    involved <- unique(c(useDatasets, unsharedDatasets))

    shared.features <- Reduce(intersect, lapply(datasets(object)[involved],
                                                rownames))
    perDatasetSelect <- list()
    for (d in involved) {
        ld <- dataset(object, d)
        if (is.null(normData(ld))) {
            warning("Dataset \"", d, "\" is not normalized, skipped")
            next
        }
        ## Make sure that all required feature meta values exist ####
        if (isH5Liger(ld)) {
            ld <- calcGeneVars.H5(ld, chunkSize = chunk,
                                  verbose = verbose)
        } else {
            featureMeta(ld, check = FALSE)$geneMeans <-
                Matrix::rowMeans(normData(ld))
            featureMeta(ld, check = FALSE)$geneVars <-
                rowVars_sparse_rcpp(normData(ld), featureMeta(ld)$geneMeans)
        }
        datasets(object, check = FALSE)[[d]] <- ld
        ## The real calculation starts here ####
        geneMeans <- featureMeta(ld)$geneMeans
        geneVars <- featureMeta(ld)$geneVars
        trx_per_cell <- cellMeta(object, "nUMI", cellIdx = object$dataset == d)
        nolan_constant <- mean((1 / trx_per_cell))
        alphathresh.corrected <- alpha / nrow(ld)
        geneMeanUpper <- geneMeans +
            stats::qnorm(1 - alphathresh.corrected / 2) *
            sqrt(geneMeans * nolan_constant / ncol(ld))
        basegenelower <- log10(geneMeans * nolan_constant)
        pass.upper <- geneVars / nolan_constant > geneMeanUpper
        pass.lower <- log10(geneVars) > basegenelower
        preselected <- data.frame(
            gene = rownames(ld),
            dataset = d,
            shared = d %in% useDatasets,
            unshared = d %in% unsharedDatasets,
            varianceDiff = geneVars - basegenelower
        )
        perDatasetSelect[[d]] <- preselected[pass.upper & pass.lower,]
    }
    perDatasetSelect <- Reduce(rbind, perDatasetSelect)
    # For shared
    shareTable <- perDatasetSelect[perDatasetSelect$gene %in% shared.features,]
    rank <- order(shareTable$varianceDiff, decreasing = TRUE)
    shareTable <- shareTable[rank,]
    line <- 1
    while (length(unique(shareTable$gene[seq(line)])) < n) line <- line + 1
    return(unique(shareTable$gene[seq(line)]))
}


.scaleH5Matrix <- function(ld, featureIdx, resultH5Path, chunk, verbose) {
    features <- rownames(ld)[featureIdx]
    geneSumSq <- featureMeta(ld)$geneSumSq[featureIdx]
    nCells <- ncol(ld)
    geneRootMeanSumSq = sqrt(geneSumSq / (nCells - 1))
    h5file <- getH5File(ld)
    safeH5Create(
        ld,
        dataPath = resultH5Path,
        dims = c(length(features), nCells),
        dtype = "double",
        chunkSize = c(length(features), chunk)
    )
    H5Apply(
        ld,
        useData = "normData",
        chunkSize = chunk,
        verbose = verbose,
        FUN = function(chunk, sparseXIdx, cellIdx, values) {
            chunk <- chunk[featureIdx, , drop = FALSE]
            chunk = as.matrix(chunk)
            chunk = sweep(chunk, 1, geneRootMeanSumSq, "/")
            rownames(chunk) <- features
            chunk[is.na(chunk)] = 0
            chunk[chunk == Inf] = 0
            h5file[[resultH5Path]][seq_along(features),
                                   cellIdx] <- chunk
        }
    )
    h5fileInfo(ld, "scaleData", check = FALSE) <- resultH5Path
    safeH5Create(
        ld,
        dataPath = paste0(resultH5Path, ".featureIdx"),
        dims = length(features),
        dtype = "int"
    )
    h5file[[paste0(resultH5Path, ".featureIdx")]][1:length(featureIdx)] <-
        featureIdx
    return(ld)
}


#' #' Perform iNMF on scaled datasets
#' #' @description
#' #' Performs integrative non-negative matrix (iNMF) factorization to return
#' #' factorized \eqn{H}, \eqn{W}, and \eqn{V} matrices. It optimizes the iNMF
#' #' objective function using block coordinate descent (alternating non-negative
#' #' least squares), where the number of factors is set by \code{k}. TODO: include
#' #' objective function equation here in documentation (using deqn)
#' #'
#' #' For each dataset, this factorization produces an \eqn{H} matrix (cells by k),
#' #' a \eqn{V} matrix (k by genes), and a shared \eqn{W} matrix (k by genes). The
#' #' \eqn{H} matrices represent the cell factor loadings. \eqn{W} is held
#' #' consistent among all datasets, as it represents the shared components of the
#' #' metagenes across datasets. The \eqn{V} matrices represent the
#' #' dataset-specific components of the metagenes.
#' #' @param object A \linkS4class{liger} object or a named list of matrix object,
#' #' where the names represents dataset names and matrices are scaled on the same
#' #' set of variable features, with rows as features and columns as cells.
#' #' @param k Inner dimension of factorization (number of factors). Run
#' #' \code{\link{suggestK}} to determine appropriate value; a general rule of
#' #' thumb is that a higher \code{k} will be needed for datasets with more
#' #' sub-structure.
#' #' @param lambda Regularization parameter. Larger values penalize
#' #' dataset-specific effects more strongly (i.e. alignment should increase as
#' #' \code{lambda} increases). Default \code{5}.
#' #' @param thresh Convergence threshold. Convergence occurs when
#' #' \eqn{|obj_0-obj|/(mean(obj_0,obj)) < thresh}. Default \code{1e-6}.
#' #' @param maxIter Maximum number of block coordinate descent iterations to
#' #' perform. Default \code{30}.
#' #' @param nrep Number of restarts to perform (iNMF objective function is
#' #' non-convex, so taking the best objective from multiple successive
#' #' initialization is recommended). For easier reproducibility, this increments
#' #' the random seed by 1 for each consecutive restart, so future factorization
#' #' of the same dataset can be run with one rep if necessary. Default \code{1}.
#' #' @param H.init Initial values to use for \eqn{H} matrices. A list object where
#' #' each element is the initial \eqn{H} matrix of each dataset. Default
#' #' \code{NULL}.
#' #' @param W.init Initial values to use for \eqn{W} matrix. A matrix object.
#' #' Default \code{NULL}.
#' #' @param V.init Initial values to use for \eqn{V} matrices. A list object where
#' #' each element is the initial \eqn{V} matrix of each dataset. Default
#' #' \code{NULL}.
#' #' @param method NNLS subproblem solver. Choose from \code{"liger"} (default
#' #' original implementation), \code{"planc"} or \code{"rcppml"}.
#' #' @param useUnshared Logical, whether to include unshared variable features and
#' #' run optimizeUANLS algorithm. Defaul \code{FALSE}. Running
#' #' \code{\link{selectGenes}} with \code{unshared = TRUE} and then running
#' #' \code{\link{scaleNotCenter}} is required.
#' #' @param seed Random seed to allow reproducible results. Default \code{1}.
#' #' @param readH5 \code{TRUE} to force reading H5 based data into memory and
#' #' conduct factorization. \code{"auto"} reads H5 dataset with less than 8000
#' #' cells. \code{FALSE} will stop users from running if H5 data presents.
#' #' @param verbose Logical. Whether to show information of the progress. Default
#' #' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' #' @param max.iters,use.unshared,rand.seed \bold{Deprecated}. See Usage section
#' #' for replacement.
#' #' @param print.obj \bold{Defunct}. Whether to print objective function values
#' #' after convergence when \code{verbose = TRUE}. Now always print when verbose.
#' #' @return \code{object} with \code{W} slot updated with the result \eqn{W}
#' #' matrix, and the \code{H} and \code{V} slots of each
#' #' \linkS4class{ligerDataset} object in the \code{datasets} slot updated with
#' #' the dataset specific \eqn{H} and \eqn{V} matrix, respectively.
#' #' @rdname runINMF_R
#' #' @examples
#' #' pbmc <- normalize(pbmc)
#' #' pbmc <- selectGenes(pbmc)
#' #' pbmc <- scaleNotCenter(pbmc)
#' #' # Only running a few iterations for fast examples
#' #' pbmc <- runINMF(pbmc, k = 20, maxIter = 2)
#' setGeneric(
#'     "runINMF_R",
#'     function(
#'         object,
#'         k,
#'         lambda = 5.0,
#'         thresh = 1e-6,
#'         maxIter = 30,
#'         nrep = 1,
#'         H.init = NULL,
#'         W.init = NULL,
#'         V.init = NULL,
#'         method = c("planc", "liger", "rcppml"),
#'         useUnshared = FALSE,
#'         seed = 1,
#'         readH5 = "auto",
#'         verbose = getOption("ligerVerbose")
#'     ) standardGeneric("runINMF_R")
#' )
#'
#' #' @rdname runINMF_R
#' setMethod(
#'     "runINMF_R",
#'     signature(object = "liger"),
#'     function(
#'         object,
#'         k,
#'         lambda = 5.0,
#'         thresh = 1e-6,
#'         maxIter = 30,
#'         nrep = 1,
#'         H.init = NULL,
#'         W.init = NULL,
#'         V.init = NULL,
#'         useUnshared = FALSE,
#'         seed = 1,
#'         readH5 = "auto",
#'         verbose = getOption("ligerVerbose")
#'     ) {
#'         .checkObjVersion(object)
#'         object <- recordCommand(object)
#'         if (isFALSE(useUnshared)) {
#'             object <- removeMissing(object, orient = "cell",
#'                                     verbose = verbose)
#'             data <- lapply(datasets(object), function(ld) {
#'                 if (is.null(scaleData(ld)))
#'                     stop("Scaled data not available. ",
#'                          "Run `scaleNotCenter(object)` first")
#'                 if (isH5Liger(ld)) {
#'                     if (!isFALSE(readH5)) {
#'                         h5d <- scaleData(ld)
#'                         if (readH5 == "auto") {
#'                             if (h5d$dims[2] <= 8000) {
#'                                 warning("Automatically reading H5 based ",
#'                                         "scaled dense matrix into memory. ",
#'                                         "Dim: ", h5d$dims[1], "x", h5d$dims[2],
#'                                         immediate. = verbose)
#'                                 return(h5d[,])
#'                             } else {
#'                                 stop("Scaled data in H5 based dataset with ",
#'                                      "more than 8000 cells will not be ",
#'                                      "automatically read into memory. Use ",
#'                                      "`readH5 = TRUE` to force reading, or ",
#'                                      "try `online_iNMF()` instead.")
#'                             }
#'                         } else if (isTRUE(readH5)) {
#'                             return(h5d[,])
#'                         } else {
#'                             stop("Can only set `readH5` to TRUE, FALSE, ",
#'                                  "or 'auto'.")
#'                         }
#'                     } else {
#'                         stop("H5 based dataset detected while `readH5` is ",
#'                              "set to FALSE.")
#'                     }
#'                 } else {
#'                     return(scaleData(ld))
#'                 }
#'             })
#'             out <- runINMF_R(
#'                 object = data,
#'                 k = k,
#'                 lambda = lambda,
#'                 thresh = thresh,
#'                 maxIter = maxIter,
#'                 nrep = nrep,
#'                 H.init = H.init,
#'                 W.init = W.init,
#'                 V.init = V.init,
#'                 useUnshared = FALSE,
#'                 seed = seed,
#'                 verbose = verbose
#'             )
#'             object@W <- out$W
#'             for (d in names(object)) {
#'                 ld <- dataset(object, d)
#'                 ld@H <- out$H[[d]]
#'                 ld@V <- out$V[[d]]
#'                 datasets(object, check = FALSE)[[d]] <- ld
#'             }
#'             object@uns$factorization$k <- k
#'             object@uns$factorization$lambda <- lambda
#'         } else {
#'             object <- runUINMF(
#'                 object = object,
#'                 k = k,
#'                 lambda = lambda,
#'                 thresh = thresh,
#'                 maxIter = maxIter,
#'                 nrep = nrep,
#'                 seed = seed,
#'                 verbose = verbose
#'             )
#'         }
#'         return(object)
#'     }
#' )
#'
#' #' @rdname runINMF_R
#' setMethod(
#'     "runINMF_R",
#'     signature(object = "list"),
#'     function(
#'         object,
#'         k,
#'         lambda = 5.0,
#'         maxIter = 30,
#'         nrep = 1,
#'         H.init = NULL,
#'         W.init = NULL,
#'         V.init = NULL,
#'         method = c("planc", "liger", "rcppml"),
#'         useUnshared = FALSE,
#'         seed = 1,
#'         readH5 = "auto",
#'         verbose = getOption("ligerVerbose")
#'     ) {
#'         # E ==> cell x gene scaled matrices
#'         E <- object
#'         nDatasets <- length(E)
#'         nCells <- sapply(E, ncol)
#'         nGenes <- nrow(E[[1]])
#'         if (k >= nGenes) {
#'             stop("Select k lower than the number of variable genes: ", nGenes)
#'         }
#'         Wm <- matrix(0, nGenes, k)
#'         Vm <- rep(list(matrix(0, nGenes, k)), nDatasets)
#'         Hm <- lapply(nCells, function(n) matrix(0, n, k))
#'
#'         bestObj <- Inf
#'         bestSeed <- seed
#'         for (i in seq(nrep)) {
#'             set.seed(seed = seed + i - 1)
#'             startTime <- Sys.time()
#'             if (!is.null(W.init))
#'                 W <- .checkInit(W.init, nCells, nGenes, k, "W")
#'             else W <- matrix(stats::runif(nGenes * k, 0, 2), nGenes, k)
#'
#'             if (!is.null(V.init)) {
#'                 V <- .checkInit(V.init, nCells, nGenes, k, "V")
#'             } else
#'                 V <- lapply(seq(nDatasets), function(i) {
#'                     matrix(stats::runif(nGenes * k, 0, 2), nGenes, k)})
#'
#'             if (!is.null(H.init)) {
#'                 H <- .checkInit(H.init, nCells, nGenes, k, "H")
#'                 H <- lapply(H, t)
#'             } else
#'                 H <- lapply(nCells, function(n) {
#'                     matrix(stats::runif(n * k, 0, 2), n, k)
#'                 })
#'
#'             if (isTRUE(verbose)) {
#'                 .log("Start iNMF with seed: ", seed + i - 1, "...")
#'                 if (maxIter > 0)
#'                     pb <- utils::txtProgressBar(0, maxIter, style = 3)
#'             }
#'             iter <- 0
#'             while (iter < maxIter) {
#'                 H <- inmfSolveH(W = W, V = V, E = E, lambda = lambda)
#'                 V <- inmfSolveV(W = W, H = H, E = E, lambda = lambda)
#'                 W <- inmfSolveW(H = H, V = V, E = E, lambda = lambda)
#'                 iter <- iter + 1
#'                 if (isTRUE(verbose) && maxIter > 0)
#'                     utils::setTxtProgressBar(pb, value = iter)
#'
#'             }
#'             if (isTRUE(verbose) && maxIter > 0) {
#'                 utils::setTxtProgressBar(pb, value = maxIter)
#'                 cat("\n")
#'             }
#'             obj <- inmf_calcObj(E, H, W, V, lambda)
#'             if (obj < bestObj) {
#'                 Wm <- W
#'                 Hm <- H
#'                 Vm <- V
#'                 bestObj <- obj
#'                 bestSeed <- seed + i - 1
#'             }
#'             endTime <- difftime(time1 = Sys.time(), time2 = startTime,
#'                                 units = "auto")
#'             if (isTRUE(verbose)) {
#'                 .log("Finished in ", endTime, " ", units(endTime),
#'                      "\nObjective error: ", bestObj)
#'                 .log("Objective: ", obj)
#'                 .log("Best results with seed ", bestSeed)
#'             }
#'         }
#'         out <- list(H = lapply(Hm, t), V = Vm, W = Wm)
#'         factorNames <- paste0("Factor_", seq(k))
#'         for (i in seq(nDatasets)) {
#'             dimnames(out$H[[i]]) <- list(factorNames, colnames(object[[i]]))
#'             dimnames(out$V[[i]]) <- list(rownames(object[[i]]), factorNames)
#'         }
#'         names(out$V) <- names(out$H) <- names(object)
#'         dimnames(out$W) <- list(rownames(object[[1]]), factorNames)
#'         return(out)
#'     }
#' )
#'
#' inmf_calcObj <- function(E, H, W, V, lambda) {
#'     # E - dgCMatrix
#'     # H, W, V - matrix
#'     obj <- 0
#'     for (i in seq_along(E)) {
#'         obj <- obj +
#'             Matrix::norm(E[[i]] - (W + V[[i]]) %*% t(H[[i]]), "F") ^ 2 +
#'             lambda*norm(V[[i]] %*% t(H[[i]]), "F") ^ 2
#'     }
#'     return(obj)
#' }
#'
#' inmfSolveH <- function(W, V, E, lambda) {
#'     H <- list()
#'     for (i in seq_along(E)) {
#'         CtC <- t(W + V[[i]]) %*% (W + V[[i]]) + lambda*(t(V[[i]]) %*% V[[i]])
#'         CtB <- as.matrix(t(W + V[[i]]) %*% E[[i]])
#'         H[[i]] <- t(RcppPlanc::bppnnls_prod(CtC, CtB))
#'     }
#'     return(H)
#' }
#'
#' inmfSolveV <- function(W, H, E, lambda) {
#'     V <- list()
#'     for (i in seq_along(E)) {
#'         CtC <- (1 + lambda)*(t(H[[i]]) %*% H[[i]])
#'         CtB <- as.matrix(t(H[[i]]) %*% t(E[[i]]))
#'         CtB <- CtB - t(H[[i]]) %*% H[[i]] %*% t(W)
#'         V[[i]] <- t(RcppPlanc::bppnnls_prod(CtC, CtB))
#'     }
#'     return(V)
#' }
#'
#' inmfSolveW <- function(H, V, E, lambda) {
#'     m <- nrow(E[[1]])
#'     k <- ncol(H[[1]])
#'     CtC <- matrix(0, k, k)
#'     CtB <- matrix(0, k, m)
#'     for (i in seq_along(E)) {
#'         CtC <- CtC + t(H[[i]]) %*% H[[i]]
#'         CtB <- CtB + as.matrix(t(H[[i]]) %*% t(E[[i]])) -
#'             t(H[[i]]) %*% H[[i]] %*% t(V[[i]])
#'     }
#'     return(t(RcppPlanc::bppnnls_prod(CtC, CtB)))
#' }


#' #' @export
#' #' @rdname runUINMF
#' #' @method runUINMF Seurat
#' #' @param datasetVar Metadata variable name that stores the dataset source
#' #' annotation. Default \code{"orig.ident"}.
#' #' @param useLayer For Seurat>=4.9.9, the name of layer to retrieve input
#' #' non-negative scaled data. Default \code{"ligerScaleData"}. For older Seurat,
#' #' always retrieve from \code{scale.data} slot.
#' #' @param assay Name of assay to use. Default \code{NULL} uses current active
#' #' assay.
#' runUINMF.Seurat <- function(
#'         object,
#'         unsharedList,
#'         k = 20,
#'         lambda = 5,
#'         datasetVar = "orig.ident",
#'         useLayer = "ligerScaleData",
#'         assay = NULL,
#'         nIteration = 30,
#'         nRandomStarts = 1,
#'         seed = 1,
#'         verbose = getOption("ligerVerbose"),
#'         ...
#' ) {
#'     mat <- .getSeuratData(object, layer = useLayer, slot = "scale.data",
#'                           assay = assay)
#'     if (any(mat < 0)) {
#'         stop("Negative data encountered for integrative Non-negative Matrix ",
#'              "Factorization. Please run `scaleNotCenter()` first.")
#'     }
#'     # the last [,1] converts data.frame to the vector/factor
#'     datasetVar <- object[[datasetVar]][,1]
#'     if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
#'     datasetVar <- droplevels(datasetVar)
#'
#'     Es <- lapply(levels(datasetVar), function(d) {
#'         as(mat[, datasetVar == d], "CsparseMatrix")
#'     })
#'     names(Es) <- levels(datasetVar)
#'     .runUINMF.list(
#'         object = Es,
#'         unsharedList = unsharedList,
#'         k = k,
#'         lambda = lambda,
#'         nIteration = nIteration,
#'         nRandomStarts = nRandomStarts,
#'         seed = seed,
#'         verbose = verbose
#'     )
#' }








#' Perform online iNMF on scaled datasets
#' @noRd
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
#' a \eqn{V} matrix (genes by k), and a shared \eqn{W} matrix (genes by k). The
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
#' @param X_new New datasets for scenario 2 or scenario 3. See detail for usage.
#' @param projection Perform data integration with scenario 3. See description.
#' Default \code{FALSE}.
#' @param W.init Optional initialization for W. See detail. Default \code{NULL}.
#' @param V.init Optional initialization for V. See detail. Default \code{NULL}.
#' @param A.init Optional initialization for A. See detail. Default \code{NULL}.
#' @param B.init Optional initialization for B. See detail. Default \code{NULL}.
#' @param k Inner dimension of factorization--number of metagenes. A value in
#' the range 20-50 works well for most analyses. Default \code{20}.
#' @param lambda Regularization parameter. Larger values penalize
#' dataset-specific effects more strongly (i.e. alignment should increase as
#' lambda increases). We recommend always using the default value except
#' possibly for analyses with relatively small differences (biological
#' replicates, male/female comparisons, etc.) in which case a lower value such
#' as 1.0 may improve reconstruction quality. Default \code{5.0}.
#' @param max.epochs The number of epochs to iterate through. See detail.
#' Default \code{5}.
#' @param miniBatch_max_iters Maximum number of block coordinate descent (HALS
#' algorithm) iterations to perform for each update of \eqn{W} and \eqn{V}.
#' Default \code{1}. Changing this parameter is not recommended.
#' @param miniBatch_size Total number of cells in each minibatch. See detail.
#' Default \code{5000}.
#' @param h5_chunk_size Chunk size of input hdf5 files. Default \code{1000}. The
#' chunk size should be no larger than the batch size.
#' @param seed Random seed to allow reproducible results. Default \code{123}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return \code{object} with \code{W} slot updated with resulting \eqn{W}
#' matrix; the \code{H}, \code{V}, \code{A} and \code{B} slots of each
#' \linkS4class{ligerDataset} object in \code{datasets} slot is updated with the
#' corresponding result matrices.
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
online_iNMFOld <- function(
        object,
        X_new = NULL,
        projection = FALSE,
        W.init = NULL,
        V.init = NULL,
        A.init = NULL,
        B.init = NULL,
        k = 20,
        lambda = 5,
        max.epochs = 5,
        miniBatch_max_iters = 1,
        miniBatch_size = 5000,
        h5_chunk_size = 1000,
        seed = 123,
        verbose = getOption("ligerVerbose")
) {
    if (!inherits(object, "liger"))
        stop("Please use a liger object.")
    .checkObjVersion(object)
    object <- recordCommand(object, dependencies = "hdf5r")
    nPrevDataset <- 0
    nNewDataset <- 0
    if (!is.null(X_new)) {
        nPrevDataset <- length(object)
        if (inherits(X_new, "liger")) {
            # concatenate `object` and `X_new`
            for (d in names(X_new)) {
                dataset(object, d, qc = FALSE) <- dataset(X_new, d)
            }
            nNewDataset <- length(X_new)
        } else if (is.list(X_new)) {
            # Add new dataset from the list to `object`
            for (i in seq_along(X_new)) {
                ld <- as.ligerDataset(X_new[[i]])
                dataset(object, names(X_new)[i], qc = FALSE) <- ld
                nNewDataset <- nNewDataset + 1
            }
            object <- normalize(object, useDatasets = names(X_new))
            object <- scaleNotCenter(object, useDatasets = names(X_new))
        } else {
            stop("`X_new` of class \"", class(X_new), "\" is not supported. ",
                 "Please use either a liger object or a named list.")
        }
    } else {
        nNewDataset <- length(object)
    }
    if (isTRUE(verbose)) .log(nNewDataset, " new datasets detected.")

    W <- getMatrix(object, "W")
    V <- getMatrix(object, "V")
    H <- getMatrix(object, "H")
    A <- lapply(datasets(object), function(x) x@A)
    B <- lapply(datasets(object), function(x) x@B)

    # Initialize minibatch sizes ####
    dataIdx <- seq_along(object)
    dataIdxNew <- seq(nPrevDataset + 1, length(object))
    dataIdxPrev <- setdiff(dataIdx, dataIdxNew)
    barcodes <- lapply(datasets(object), colnames)
    nCells <- unlist(lapply(datasets(object), ncol)) # vector of ncell per data
    nCellsNew <- nCells[dataIdxNew]
    minibatchSizes <- rep(0, length(object))
    nGenes <- length(varFeatures(object))
    for (i in dataIdxNew) {
        minibatchSizes[i] <- round(nCells[i] / sum(nCellsNew) *
                                       miniBatch_size)
        if (minibatchSizes[i] > nCells[i]) {
            stop("Number of cells to be sampled (n = ",  minibatchSizes[i],
                 ") is greater than the size of input dataset ", i, " (n = ",
                 nCells[i], "). Please use a smaller `miniBatch_size`.")
        }
    }
    minibatchSizesOrig <- minibatchSizes

    AOld <- rep(list(NULL), length(object))
    BOld <- rep(list(NULL), length(object))
    if (isFALSE(projection)) {
        if (!is.null(seed)) set.seed(seed)
        # Initialize all kinds of matrices ####
        if (is.null(X_new)) {
            if (isTRUE(verbose)) .log("Scenario 1, init arguments ignored.")

            W <- matrix(stats::runif(nGenes * k, 0, 2), nGenes, k)
            # TODO: W <- W / sqrt(colSums(W ^ 2))
            for (j in seq(k)) W[, j] <- W[, j] / sqrt(sum(W[, j] ^ 2))
            V <- list()
            for (i in dataIdx) {
                # VInitIdx <- sample(nCells[i], k)
                VInitIdx <- sample_cpp(nCells[i], k)#as.numeric(RcppPlanc::sample_cpp(nCells[i]))[1:k] + 1
                # pick k sample from datasets as initial H matrix
                V[[i]] = as.matrix(scaleData(object, i)[1:nGenes, VInitIdx])
                for (j in seq(k)) {
                    # normalize columns of dictionaries
                    V[[i]][, j] = V[[i]][, j] / sqrt(sum(V[[i]][, j]^2))
                }
            }
            # "Old"s for saving information older than 2 epochs
            A <- AOld <- rep(list(matrix(0, k, k)), nNewDataset)
            B <- BOld <- rep(list(matrix(0, nGenes, k)), nNewDataset)
        } else {
            if (isTRUE(verbose)) .log("Scenario 2, initiating parameters")

            if (!is.null(W.init)) W <- .checkInit(W.init, NULL, nGenes, k, "W")
            W <- .checkMatrixValid(W, k, name = "W")

            if (!is.null(V.init)) {
                V <- .checkInit(V.init, nCells[dataIdxPrev], nGenes, k, "V")
            }
            V[dataIdxPrev] <- lapply(V[dataIdxPrev],
                                     function(v)
                                         .checkMatrixValid(v, k, name = "V"))
            for (i in dataIdxNew) {
                VInitIdx <- sample_cpp(nCells[i], k)
                # initialize the Vi for new dataset
                V[[i]] <- as.matrix(scaleData(object, i)[1:nGenes, VInitIdx])
                for (j in seq(k))
                    V[[i]][, j] <- V[[i]][, j] / sqrt(sum(V[[i]][, j]^2))
            }

            if (!is.null(A.init)) A[dataIdxPrev] <- A.init
            if (!is.null(B.init)) B[dataIdxPrev] <- B.init
            A[dataIdxNew] <- rep(list(matrix(0, k, k)), nNewDataset)
            B[dataIdxNew] <- rep(list(matrix(0, nGenes, k)), nNewDataset)
            AOld[dataIdxNew] <- rep(list(matrix(0, k, k)), nNewDataset)
            # save information older than 2 epochs
            BOld[dataIdxNew] <- rep(list(matrix(0, nGenes, k)), nNewDataset)
            # save information older than 2 epochs
        }

        iter <- 1
        # intialize the number of epoch for each dataset
        epoch <- rep(0, length(object))
        # intialize the previous number of epoch for each dataset
        epochPrev <- rep(0, length(object))
        epochNext <- rep(FALSE, length(object))
        sqrtLambda <- sqrt(lambda)

        if (isTRUE(verbose)) .log("Initialization done.")

        # chunk permutation: shuffle cell index by H5 chunks
        allIdx <- rep(list(NULL), length(object))
        for (i in dataIdxNew) {
            allIdx[[i]] <- .permuteChunkIdx(object, i, h5_chunk_size)
        }

        totalIters <- floor(sum(nCellsNew) * max.epochs / miniBatch_size)
        if (isTRUE(verbose)) {
            .log("Starting Online iNMF...")
            pb <- utils::txtProgressBar(1, totalIters + 1, style = 3)
        }
        while (epoch[dataIdxNew[1]] < max.epochs) {
            # track epochs
            # indices of samples in each dataest used for this iteration
            minibatchIdx <- rep(list(NULL), length(object))
            # check if the size of the last mini-batch equals to pre-specified
            # mini-batch size
            if (max.epochs*nCellsNew[1] >= iter*minibatchSizes[dataIdxNew[1]]) {
                for (i in dataIdxNew) {
                    # calculate the current epoch
                    epoch[i] <- (iter*minibatchSizes[i]) %/% nCells[i]
                    # if current iteration cycles through the data and start
                    # a new cycle
                    batchStartIdx <- ((iter - 1)*minibatchSizes[i]) %% nCells[i] + 1
                    if (epochPrev[i] != epoch[i]) {
                        # If entering a new epoch
                        epochNext[i] <- TRUE
                        epochPrev[i] <- epoch[i]
                        minibatchIdx[[i]] <- allIdx[[i]][batchStartIdx:nCells[i]]
                        # print(paste0(
                        #     "start: ", batchStartIdx-1, ", end: ", nCells[i]-1, " + "
                        # ))
                        allIdx[[i]] <- .permuteChunkIdx(object, i, h5_chunk_size)
                        if ((iter * minibatchSizes[i]) %% nCells[i] != 0) {
                            # print(paste0("start: 0, end: ", (iter * minibatchSizes[i]) %% nCells[i] - 1))
                            minibatchIdx[[i]] <- c(minibatchIdx[[i]],
                                                   allIdx[[i]][seq((iter * minibatchSizes[i]) %% nCells[i])])
                        }
                    } else {
                        # if current iter stays within a epoch
                        # print(
                        #     paste0("start: ", batchStartIdx-1, ", end: ", ((iter * minibatchSizes[i]) %% nCells[i]) - 1)
                        # )
                        minibatchIdx[[i]] <- allIdx[[i]][batchStartIdx:((iter * minibatchSizes[i]) %% nCells[i])]
                    }
                }
            } else {
                # last iteration
                for (i in dataIdxNew) {
                    minibatchSizes[i] <- max.epochs*nCells[i] - (iter - 1)*minibatchSizes[i]
                    minibatchIdx[[i]] <- (((iter - 1)*minibatchSizesOrig[i]) %% nCells[i] + 1):nCells[i]
                }
                epoch[dataIdxNew[1]] <- max.epochs
            }
            if (length(minibatchIdx[[dataIdxNew[1]]]) != minibatchSizesOrig[dataIdxNew[1]])
                next

            X_minibatch = rep(list(NULL), length(object))
            for (i in dataIdxNew) {
                X_minibatch[[i]] = as.matrix(scaleData(object, i)[1:nGenes, minibatchIdx[[i]]])
            }

            # update H_i by ANLS Hi_minibatch[[i]]
            H_minibatch = rep(list(NULL), length(object))
            for (i in dataIdxNew) {
                H_minibatch[[i]] <- solveNNLS(
                    rbind(W + V[[i]], sqrtLambda * V[[i]]),
                    rbind(X_minibatch[[i]],
                          matrix(0, nGenes, minibatchSizes[i]))
                )
            }

            # updata A and B matrices
            if (iter == 1) {
                scale_param <- rep(0, length(object))
            } else if (iter == 2) {
                scale_param <- c(rep(0, nPrevDataset),
                                 rep(1, nNewDataset) /
                                     minibatchSizes[dataIdxNew])
            } else {
                scale_param <- c(rep(0, nPrevDataset),
                                 rep((iter - 2) / (iter - 1), nNewDataset))
            }

            for (i in dataIdxNew) {
                # print(scale_param[i])
                if (epoch[dataIdxNew[1]] > 0 & epochNext[dataIdxNew[1]]) {
                    # remove information older than 2 epochs
                    A[[i]] <- A[[i]] - AOld[[i]]
                    AOld[[i]] <- scale_param[i] * A[[i]]
                    B[[i]] <- B[[i]] - BOld[[i]]
                    BOld[[i]] <- scale_param[i] * B[[i]]
                } else {
                    AOld[[i]] <- scale_param[i] * AOld[[i]]
                    BOld[[i]] <- scale_param[i] * BOld[[i]]
                }

                # HiHit
                A[[i]] <- scale_param[i] * A[[i]] +
                    H_minibatch[[i]] %*% t(H_minibatch[[i]]) / minibatchSizes[i]
                diag(A[[i]])[diag(A[[i]]) == 0] <- 1e-15
                # XiHit
                B[[i]] <- scale_param[i] * B[[i]] +
                    X_minibatch[[i]] %*% t(H_minibatch[[i]]) / minibatchSizes[i]
                # print(A[[i]][1:4,1:4])
                # print(B[[i]][1:4, 1:4])
            }

            # update W, V_i by HALS
            iter_miniBatch <- 1
            max_iters_miniBatch <- miniBatch_max_iters


            while (iter_miniBatch <= max_iters_miniBatch) {
                # update W
                for (j in seq(k)) {
                    W_update_numerator <- rep(0, nGenes)
                    W_update_denominator <- 0
                    for (i in dataIdx) {
                        W_update_numerator <- W_update_numerator +
                            B[[i]][, j] - (W + V[[i]]) %*% A[[i]][, j]
                        W_update_denominator <- W_update_denominator +
                            A[[i]][j,j]
                    }

                    W[, j] <- nonneg(
                        W[, j] + W_update_numerator / W_update_denominator
                    )
                }

                # update V_i
                for (j in seq(k)) {
                    for (i in dataIdxNew) {
                        V[[i]][, j] <- nonneg(
                            V[[i]][, j] +
                                (B[[i]][, j] - (W + (1 + lambda)*V[[i]]) %*% A[[i]][, j]) /
                                ((1 + lambda) * A[[i]][j, j])
                        )
                    }
                }
                iter_miniBatch <- iter_miniBatch + 1
            }
            # reset epoch change indicator
            epochNext <- rep(FALSE, length(object))
            iter <- iter + 1
            if (isTRUE(verbose)) {
                utils::setTxtProgressBar(pb = pb, value = iter)
            }
        }

        if (isTRUE(verbose)) {
            cat("\n")
            .log("Calculating metagene loadings...")
        }
        H <- rep(list(NULL), length(object))
        for (i in dataIdx) {
            batchIdxs <- .batchCellIdx(nCells[i], miniBatch_size)
            for (j in seq_along(batchIdxs)) {
                cellIdx <- batchIdxs[[j]]
                batch <- as.matrix(scaleData(object, i)[1:nGenes, cellIdx])
                H[[i]] <- cbind(
                    H[[i]],
                    solveNNLS(rbind(W + V[[i]], sqrtLambda * V[[i]]),
                              rbind(batch, matrix(0, nGenes,length(cellIdx))))
                )
            }
            colnames(H[[i]]) <- barcodes[[i]]
        }
    } else {
        if (isTRUE(verbose))
            .log("Scenario 3, metagene projection")
        if (!is.null(W.init)) W <- W.init
        H[dataIdxNew] <- rep(list(NULL), nNewDataset)
        V[dataIdxNew] <- rep(list(NULL), nNewDataset)
        for (i in dataIdxNew) {
            batchIdxs <- .batchCellIdx(nCells[i], miniBatch_size)
            for (j in seq_along(batchIdxs)) {
                batch <- as.matrix(scaleData(object, i)[1:nGenes, batchIdxs[[j]]])
                H[[i]] <- cbind(H[[i]], solveNNLS(W, batch))
            }
            colnames(H[[i]]) <- barcodes[[i]]
            V[[i]] <- matrix(0, nGenes, k)
        }
    }
    factorNames <- paste0("Factor_", seq(k))
    colnames(W) <- factorNames
    rownames(W) <- varFeatures(object)
    object@W <- W
    for (i in dataIdx) {
        ld <- dataset(object, i)
        rownames(H[[i]]) <- factorNames
        ld@H <- H[[i]]
        colnames(V[[i]]) <- factorNames
        rownames(V[[i]]) <- varFeatures(object)
        ld@V <- V[[i]]
        #dimnames(A[[i]]) <- list(factorNames, factorNames)
        ld@A <- A[[i]]
        #colnames(B[[i]]) <- factorNames
        ld@B <- B[[i]]
        datasets(object, check = FALSE)[[i]] <- ld
    }
    object@uns$factorization$k <- k
    object@uns$factorization$lambda <- lambda
    object
}


.permuteChunkIdx <- function(object, dataset, chunkSize = NULL) {
    ld <- dataset(object, dataset)
    if (isH5Liger(ld)) chunkSize <- scaleData(ld)$chunk_dims[2]
    nChunks <- ceiling(ncol(ld) / chunkSize)
    chunkIdx <- sample_cpp(nChunks, nChunks)#as.numeric(RcppPlanc::sample_cpp(nChunks)) + 1
    # chunkIdx <- sample(nChunks, nChunks)
    unlist(lapply(chunkIdx, function(i) {
        if (i != nChunks) seq(1 + chunkSize * (i - 1), i * chunkSize)
        else seq((1 + chunkSize * (i - 1)), ncol(ld))
    }), use.names = FALSE)
}

.batchCellIdx <- function(nCell, size) {
    nBatch <- ceiling(nCell / size)
    result <- list()
    for (i in seq(nBatch)) {
        result[[i]] <- if (i != nBatch) seq((i - 1) * size + 1, i * size)
        else seq((i - 1) * size + 1, nCell)
    }
    result
}

.checkMatrixValid <- function(m, k, name) {
    if (is.null(m)) {
        stop("Matrix ", name, " has to be provided, either from ",
             "pre-calculation or supply by `", name, ".init`.")
    }
    if (ncol(m) != k) {
        stop("Factors from matrix ", name,
             " do not match specification of `k` (", k, ").")
    }
    m
}

nonneg <- function(x, eps = 1e-16) {
    x[x < eps] <- eps
    return(x)
}
