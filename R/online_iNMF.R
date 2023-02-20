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
#' @param verbose Logical. Whether to show information of the progress.
#' Default \code{TRUE}.
#' @return \code{object} with \code{W} slot updated with resulting \eqn{W}
#' matrix; the \code{H}, \code{V}, \code{A} and \code{B} slots of each
#' \linkS4class{ligerDataset} object in \code{datasets} slot is updated with the
#' corresponding result matrices.
#' @export
online_iNMF <- function(
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
        seed = 123,
        verbose = TRUE
) {
    if (!inherits(object, "liger"))
        stop("Please use a liger object.")
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
            for (i in length(X_new)) {
                if (inherits(X_new[[i]], "liger")) {
                    for (d in names(X_new[[i]])) {
                        dataset(object, d, qc = FALSE) <- dataset(X_new[[i]], d)
                    }
                    nNewDataset <- nNewDataset + length(X_new[[i]])
                } else {
                    ld <- as.ligerDataset(X_new)
                    dataset(object, names(X_new)[i], qc = FALSE) <- ld
                    nNewDataset <- nNewDataset + 1
                }
            }
            object <- normalize(object)
            object <- scaleNotCenter(object)
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
    nGenes <- length(var.features(object))
    for (i in dataIdxNew) {
        minibatchSizes[i] <- round(nCells[i] / sum(nCells[dataIdxNew]) *
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

            W <- matrix(runif(nGenes * k, 0, 2), nGenes, k)
            # TODO: W <- W / sqrt(colSums(W ^ 2))
            for (j in seq(k)) W[, j] <- W[, j] / sqrt(sum(W[, j] ^ 2))

            V <- list()
            for (i in dataIdx) {
                VInitIdx <- sample(nCells[i], k)
                # pick k sample from datasets as initial H matrix
                V[[i]] = scale.data(object, i)[1:nGenes, VInitIdx]
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
            W <- .checkMatrixValid(W, k, nGenes, name = "W")

            if (!is.null(V.init)) {
                V <- .checkInit(V.init, nCells[dataIdxPrev], nGenes, k, "V")
            }
            V[dataIdxPrev] <- lapply(V[dataIdxPrev],
                                     function(v)
                                         .checkMatrixValid(v, k, name = "V"))
            for (i in dataIdxNew) {
                VInitIdx <- sample(nCells[i], k)
                # initialize the Vi for new dataset
                V[[i]] <- scale.data(object, i)[1:nGenes, VInitIdx]
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

        .log("Initialization done.")

        # chunk permutation: shuffle cell index by H5 chunks
        allIdx <- rep(list(NULL), length(object))
        for (i in dataIdxNew) {
            allIdx[[i]] <- .permuteChunkIdx(object, i, 1000)
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
                        allIdx[[i]] <- .permuteChunkIdx(object, i, 1000)
                        if ((iter * minibatchSizes[i]) %% nCells[i] != 0) {
                            minibatchIdx[[i]] <- c(minibatchIdx[[i]],
                                                  allIdx[[i]][seq((iter * minibatchSizes[i]) %% nCells[i])])
                        }
                    } else {
                        # if current iter stays within a epoch
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
                X_minibatch[[i]] = scale.data(object, i)[1:nGenes, minibatchIdx[[i]]]
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
                setTxtProgressBar(pb = pb, value = iter)
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
                batch <- scale.data(object, i)[1:nGenes, cellIdx]
                H[[i]] = cbind(H[[i]],
                               solveNNLS(rbind(W + V[[i]], sqrtLambda * V[[i]]),
                                         rbind(batch, matrix(0, nGenes,length(cellIdx))))
                               )
            }
            colnames(H[[i]]) = barcodes[[i]]
        }
        rownames(W) <- var.features(object)
        colnames(W) <- NULL

        for (i in dataIdx) {
            rownames(V[[i]]) <- var.features(object)
            colnames(V[[i]]) <- NULL
        }
    } else {
        if (isTRUE(verbose))
            .log("Scenario 3, metagene projection")
        if (!is.null(W.init)) W <- W.init
        H[dataIdxNew] = rep(list(NULL), nNewDataset)
        V[dataIdxNew] = rep(list(NULL), nNewDataset)
        for (i in dataIdxNew) {
            batchIdxs <- .batchCellIdx(nCells[i], miniBatch_size)
            for (j in seq_along(batchIdxs)) {
                batch <- scale.data(object, i)[1:nGenes, batchIdxs[[j]]]
                H[[i]] <- cbind(H[[i]], solveNNLS(W, batch))
            }
            colnames(H[[i]]) = barcodes[[i]]
            V[[i]] = matrix(0, nGenes, k)
        }
    }
    object@W <- W
    for (i in dataIdx) {
        ld <- dataset(object, i)
        ld@H <- H[[i]]
        ld@V <- V[[i]]
        ld@A <- A[[i]]
        ld@B <- B[[i]]
        datasets(object, check = FALSE)[[i]] <- ld
    }
    object@k <- ncol(W)
    object
}


.permuteChunkIdx <- function(object, dataset, chunkSize = NULL) {
    ld <- dataset(object, dataset)
    if (isH5Liger(ld)) chunkSize <- scale.data(ld)$chunk_dims[2]
    nChunks <- ceiling(ncol(ld)/chunkSize)
    chunkIdx <- sample(nChunks, nChunks)
    unlist(lapply(chunkIdx, function(i) {
        if (i != nChunks) seq(1 + chunkSize*(i - 1), i*chunkSize)
        else seq((1 + chunkSize*(i - 1)), ncol(ld))
    }), use.names = FALSE)
}

.batchCellIdx <- function(nCell, size) {
    nBatch <- ceiling(nCell / size)
    result <- list()
    for (i in seq(nBatch)) {
        result[[i]] <- if (i != nBatch) seq((i - 1)*size + 1, i*size)
        else seq((i - 1)*size + 1, nCell)
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
