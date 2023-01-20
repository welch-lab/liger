online_iNMF <- function(
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
        seed = 123,
        verbose = TRUE
) {
    if (!inherits(object, "liger"))
        stop("Please use a liger object.")
    W <- object@W
    V <- lapply(datasets(object), function(x) x@V)
    H <- lapply(datasets(object), function(x) x@H)
    A <- lapply(datasets(object), function(x) x@A)
    B <- lapply(datasets(object), function(x) x@B)
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
        }
    } else {
        nNewDataset <- length(object)
    }
    if (isTRUE(verbose))
        message(date(), " ... ", nNewDataset, " new datasets detected.")

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
            if (isTRUE(verbose))
                message(date(), " ... Scenario 1, init arguments ignored.")
            W <- matrix(runif(nGenes * k, 0, 2), nGenes, k)
            for (j in seq(k)) W[, j] <- W[, j] / sqrt(sum(W[, j] ^ 2))
            V <- list()
            for (i in dataIdx) {
                VInitIdx <- sample(nCells[i], k)
                # pick k sample from datasets as initial H matrix
                V[[i]] = scale.data(object, i)[1:nGenes, VInitIdx]
            }

            # normalize the columns of H_i, H_s matrices
            for (j in seq(k)) {
                for (i in dataIdx) {
                    # normalize columns of dictionaries
                    V[[i]][, j] = V[[i]][, j] / sqrt(sum(V[[i]][, j]^2))
                }
            }
            H <- rep(list(NULL), nNewDataset)
            H_minibatch <- list()
            A <- rep(list(matrix(0, k, k)), nNewDataset)
            B <- rep(list(matrix(0, nGenes, k)), nNewDataset)
            AOld <- rep(list(matrix(0, k, k)), nNewDataset)
            # save information older than 2 epochs
            BOld <- rep(list(matrix(0, nGenes, k)), nNewDataset)
        } else {
            if (isTRUE(verbose))
                message(date(), " ... Scenario 2, initiating parameters")
            if (!is.null(W.init)) W <- W.init
            if (!is.null(V.init)) V <- V.init
            for (i in dataIdxNew) {
                VInitIdx <- sample(nCells[i], k)
                # initialize the Vi for new dataset
                V[[i]] <- scale.data(object, i)[1:nGenes, VInitIdx]
                for (j in seq(k))
                    V[[i]][, j] <- V[[i]][, j] / sqrt(sum(V[[i]][, j]^2))
            }
            if (!is.null(H.init)) H[dataIdxPrev] <- H.init
            H[dataIdxNew] <- rep(list(NULL), nNewDataset)
            HMinibatch <- list()
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

        message(date(), " ... Initialization done.")

        # chunk permutation: shuffle cell index by H5 chunks
        allIdx <- list()
        for (i in seq_along(dataIdxNew)) {
            allIdx[[i]] <- .permuteChunkIdx(object, i, 1000)
        }

        totalIters <- floor(sum(nCellsNew) * max.epochs / miniBatch_size)
        if (isTRUE(verbose)) {
            message(date(), " ... Starting Online iNMF...")
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
                    B[[i]] <- B[[i]] - BOld[[i]]
                }
                AOld[[i]] <- scale_param[i] * AOld[[i]]
                BOld[[i]] <- scale_param[i] * BOld[[i]]
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
            message(date(), " ... Calculating metagene loadings...")
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
            message(date(), " ... Scenario 3, metagene projection")
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
