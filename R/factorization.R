#' Perform online iNMF on scaled datasets
#'
#' @description
#' Perform online integrative non-negative matrix factorization to represent multiple single-cell datasets
#' in terms of H, W, and V matrices. It optimizes the iNMF objective function using online learning (non-negative
#' least squares for H matrix, hierarchical alternating least squares for W and V matrices), where the
#' number of factors is set by k. The function allows online learning in 3 scenarios: (1) fully observed datasets;
#' (2) iterative refinement using continually arriving datasets; and (3) projection of new datasets without updating
#' the existing factorization. All three scenarios require fixed memory independent of the number of cells.
#'
#' For each dataset, this factorization produces an H matrix (cells by k), a V matrix (k by genes),
#' and a shared W matrix (k by genes). The H matrices represent the cell factor loadings.
#' W is identical among all datasets, as it represents the shared components of the metagenes
#' across datasets. The V matrices represent the dataset-specific components of the metagenes.
#'
#' @param object \code{liger} object with data stored in HDF5 files. Should normalize, select genes, and scale before calling.
#' @param X_new List of new datasets for scenario 2 or scenario 3. Each list element should be the name of an HDF5 file.
#' @param projection Perform data integration by shared metagene (W) projection (scenario 3). (default FALSE)
#' @param W.init Optional initialization for W. (default NULL)
#' @param V.init Optional initialization for V (default NULL)
#' @param H.init Optional initialization for H (default NULL)
#' @param A.init Optional initialization for A (default NULL)
#' @param B.init Optional initialization for B (default NULL)
#' @param k Inner dimension of factorization--number of metagenes (default 20). A value in the range 20-50 works well for most analyses.
#' @param lambda Regularization parameter. Larger values penalize dataset-specific effects more
#'   strongly (ie. alignment should increase as lambda increases). We recommend always using the default value except
#'   possibly for analyses with relatively small differences (biological replicates, male/female comparisons, etc.)
#'   in which case a lower value such as 1.0 may improve reconstruction quality. (default 5.0).
#' @param max.epochs Maximum number of epochs (complete passes through the data). (default 5)
#' @param miniBatch_max_iters Maximum number of block coordinate descent (HALS algorithm) iterations to perform for
#' each update of W and V (default 1). Changing this parameter is not recommended.
#' @param miniBatch_size Total number of cells in each minibatch (default 5000). This is a reasonable default, but a smaller value
#' such as 1000 may be necessary for analyzing very small datasets. In general, minibatch size should be no larger than the number
#' of cells in the smallest dataset.
#' @param h5_chunk_size Chunk size of input hdf5 files (default 1000). The chunk size should be no larger than the batch size.
#' @param seed Random seed to allow reproducible results (default 123).
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return \code{liger} object with H, W, V, A and B slots set.
#'
#' @import hdf5r
#'
#' @export
#' @examples
#' \dontrun{
#' # Requires preprocessed liger object
#' # Get factorization using 20 factors and mini-batch of 5000 cells
#' # (default setting, can be adjusted for ideal results)
#' ligerex <- online_iNMF(ligerex, k = 20, lambda = 5, miniBatch_size = 5000)
#' }

online_iNMF <- function(object,
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
                        verbose = TRUE) {
    if (!is.null(X_new)) {
        # if there is new dataset
        raw.data_prev = object@raw.data
        norm.data_prev = object@norm.data
        h5file.info_prev = object@h5file.info
        scale.data_prev = object@scale.data
        cell.data_prev = object@cell.data
        names(raw.data_prev) = names(object@raw.data)

        # assuming only one new dataset arrives at a time
        raw.data = c()
        norm.data = c()
        h5file.info = c()
        scale.data = c()
        cell.data = c()
        for (i in 1:length(X_new)) {
            raw.data = c(raw.data, X_new[[i]]@raw.data)
            norm.data = c(norm.data, X_new[[i]]@norm.data)
            h5file.info = c(h5file.info, X_new[[i]]@h5file.info)
            scale.data = c(scale.data, X_new[[i]]@scale.data)
            cell.data = rbind(cell.data, X_new[[i]]@cell.data)
        }
        object@raw.data = raw.data
        object@norm.data = norm.data
        object@h5file.info = h5file.info
        object@scale.data = scale.data
        object@cell.data = cell.data

        # check whether X_new needs to be processed
        for (i in 1:length(object@raw.data)) {
            if (class(object@raw.data[[i]])[1] == "H5File") {
                processed = object@raw.data[[i]]$exists("scale.data")
            } else {
                processed = !is.null(X_new[[i]]@scale.data)
            }

            if (processed) {
                if (verbose) {
                    cat("New dataset",
                        i,
                        "already preprocessed.",
                        "\n")
                }
            } else {
                if (verbose) {
                    cat("New dataset",
                        i,
                        "not preprocessed. Preprocessing...",
                        "\n")
                }
                object = normalize(object, chunk = h5_chunk_size)
                object = scaleNotCenter(object,
                                        remove.missing = TRUE,
                                        chunk = h5_chunk_size)
                if (verbose) {
                    cat("New dataset", i, "Processed.", "\n")
                }
            }
        }


        object@raw.data = c(raw.data_prev, object@raw.data)
        object@norm.data = c(norm.data_prev, object@norm.data)
        object@h5file.info = c(h5file.info_prev, object@h5file.info)
        object@scale.data = c(scale.data_prev, object@scale.data)
        object@cell.data = rbind(cell.data_prev, object@cell.data)
        # k x gene -> gene x k & cell x k-> k x cell
        object@W = t(object@W)
        object@V = lapply(object@V, t)
        object@H = lapply(object@H, t)
    }

    for (i in 1:length(object@raw.data)) {
        if (class(object@raw.data[[i]])[1] != "H5File")
            object@scale.data[[i]] = t(object@scale.data[[i]])
    }

    ## extract required information and initialize algorithm
    num_files = length(object@raw.data) # number of total input hdf5 files
    num_prev_files = 0 # number of input hdf5 files processed in last step
    num_new_files = 0 # number of new input hdf5 files since last step
    if (is.null(X_new)) {
        num_prev_files = 0 # start from scratch
        num_new_files = num_files
    } else {
        num_new_files = length(X_new)
        num_prev_files = num_files - num_new_files
        if (verbose) {
            cat(num_new_files, "new datasets detected.", "\n")
        }
    }

    file_idx = 1:num_files # indices for all input files
    file_idx_new = (num_prev_files + 1):num_files # indices only for new input files
    file_idx_prev = setdiff(file_idx, file_idx_new)

    vargenes = object@var.genes
    file_names = names(object@raw.data)
    gene_names = vargenes # genes selected for analysis
    num_genes = length(vargenes) # number of the selected genes

    cell_barcodes = list() # cell barcodes for each dataset
    for (i in file_idx) {
        cell_barcodes[[i]] = rownames(object@cell.data)[object@cell.data$dataset == file_names[i]]
    }
    num_cells = unlist(lapply(cell_barcodes, length)) # number of cells in each dataset
    num_cells_new = num_cells[(num_prev_files + 1):num_files]
    minibatch_sizes = rep(0, num_files)

    for (i in file_idx_new) {
        minibatch_sizes[i] = round((num_cells[i] / sum(num_cells[file_idx_new])) * miniBatch_size)
        if (minibatch_sizes[i] > num_cells[i]) {
            stop(
                "\nNumber of cells to be sampled (n=",
                minibatch_sizes[i],
                ") is larger than the size of input dataset ",
                i,
                " (n=",
                num_cells[i],
                ").",
                "\nPlease use a smaller mini-batch size."
            )
        }
    }
    minibatch_sizes_orig = minibatch_sizes

    if (!projection) {
        if (!is.null(seed)) {
            set.seed(seed)
        }

        # W matrix initialization
        if (is.null(X_new)) {
            object@W = matrix(abs(runif(num_genes * k, 0, 2)), num_genes, k)
            for (j in 1:k) {
                object@W[, j] = object@W[, j] / sqrt(sum(object@W[, j] ^ 2))
            }
        } else {
            object@W = if (!is.null(W.init))
                W.init
            else
                object@W
        }
        # V_i matrix initialization
        if (is.null(X_new)) {
            object@V = list()
            for (i in file_idx) {
                V_init_idx = sample(1:num_cells_new[i], k)
                # pick k sample from datasets as initial H matrix
                object@V[[i]] = object@scale.data[[i]][1:num_genes, V_init_idx]
                #object@V[[i]] = matrix(data = abs(x = runif(n = num_genes * k, min = 0, max = 2)),
                #                       nrow = num_genes,
                #                       ncol = k)
            }

            # normalize the columns of H_i, H_s matrices
            for (j in 1:k) {
                for (i in file_idx) {
                    # normalize columns of dictionaries
                    object@V[[i]][, j] = object@V[[i]][, j] / sqrt(sum(object@V[[i]][, j] ^
                                                                           2))
                }
            }
        } else {
            # if previous Vs are provided
            object@V[file_idx_prev] = if (!is.null(V.init))
                V.init
            else
                object@V
            V_init_idx = list()
            for (i in file_idx_new) {
                V_init_idx = sample(1:num_cells[i], k)
                object@V[[i]] = object@scale.data[[i]][1:num_genes, V_init_idx] # initialize the Vi for new dataset
                for (j in 1:k) {
                    object@V[[i]][, j] = object@V[[i]][, j] / sqrt(sum(object@V[[i]][, j] ^
                                                                           2))
                }
            }
        }
        # H_i matrices initialization
        if (is.null(X_new)) {
            object@H = rep(list(NULL), num_files)
            H_minibatch = list()
        } else {
            # if previous Hs are provided
            object@H[file_idx_prev] = if (!is.null(H.init))
                H.init
            else
                object@H
            object@H[file_idx_new] = rep(list(NULL), num_new_files)
            H_minibatch = list()
        }
        # A = HiHi^t, B = XiHit
        A_old = list()
        B_old = list()

        if (is.null(X_new)) {
            object@A = rep(list(matrix(0, k, k)), num_new_files)
            object@B = rep(list(matrix(0, num_genes, k)), num_new_files)
            A_old = rep(list(matrix(0, k, k)), num_new_files)
            # save information older than 2 epochs
            B_old = rep(list(matrix(0, num_genes, k)), num_new_files)
            # save information older than 2 epochs

        } else {
            object@A[file_idx_prev] = if (!is.null(A.init))
                A.init
            else
                object@A
            object@B[file_idx_prev] = if (!is.null(B.init))
                B.init
            else
                object@B
            A_old[file_idx_prev] = rep(list(NULL), num_prev_files)
            B_old[file_idx_prev] = rep(list(NULL), num_prev_files)
            object@A[(num_prev_files + 1):num_files] = rep(list(matrix(0, k, k)),
                                                           num_new_files)
            object@B[(num_prev_files + 1):num_files] = rep(list(matrix(0,
                                                                       num_genes,
                                                                       k)),
                                                           num_new_files)
            A_old[(num_prev_files + 1):num_files] = rep(list(matrix(0, k, k)),
                                                        num_new_files)
            # save information older than 2 epochs
            B_old[(num_prev_files + 1):num_files] = rep(list(matrix(0, k, k)),
                                                        num_new_files)
            # save information older than 2 epochs
        }

        iter = 1
        epoch = rep(0, num_files)
        # intialize the number of epoch for each dataset
        epoch_prev = rep(0, num_files)
        # intialize the previous number of epoch for each dataset
        epoch_next = rep(FALSE, num_files)
        sqrt_lambda = sqrt(lambda)
        total_time = 0
        # track the total amount of time used for the online learning


        num_chunks = rep(NULL, num_files)
        chunk_idx = rep(list(NULL), num_files)
        all_idx = rep(list(NULL), num_files)

        # chunk permutation
        for (i in file_idx_new) {
            num_chunks[i] = ceiling(num_cells[i] / h5_chunk_size)
            chunk_idx[[i]] = sample(1:num_chunks[i], num_chunks[i])
            # idx in the first chunk
            if (chunk_idx[[i]][1] != num_chunks[i]) {
                all_idx[[i]] = (1 + h5_chunk_size * (chunk_idx[[i]][1] - 1)):(chunk_idx[[i]][1] *
                                                                                  h5_chunk_size)
            } else {
                all_idx[[i]] = (1 + h5_chunk_size * (chunk_idx[[i]][1] - 1)):(num_cells[i])
            }

            for (j in chunk_idx[[i]][-1]) {
                if (j != num_chunks[i]) {
                    all_idx[[i]] = c(all_idx[[i]],
                                     (1 + h5_chunk_size * (j - 1)):(j * h5_chunk_size))
                } else {
                    all_idx[[i]] = c(all_idx[[i]],
                                     (1 + h5_chunk_size * (j - 1)):num_cells[i])
                }
            }
        }

        total.iters = floor(sum(num_cells_new) * max.epochs / miniBatch_size)
        if (verbose) {
            cat("Starting Online iNMF...", "\n")
            pb <-
                txtProgressBar(min = 1,
                               max = total.iters + 1,
                               style = 3)
        }

        while (epoch[file_idx_new[1]] < max.epochs) {
            # track epochs
            minibatch_idx = rep(list(NULL), num_files) # indices of samples in each dataest used for this iteration
            if ((max.epochs * num_cells_new[1] - (iter - 1) * minibatch_sizes[file_idx_new[1]]) >= minibatch_sizes[file_idx_new[1]]) {
                # check if the size of the last mini-batch == pre-specified mini-batch size
                for (i in file_idx_new) {
                    epoch[i] = (iter * minibatch_sizes[i]) %/% num_cells[i] # caculate the current epoch
                    if ((epoch_prev[i] != epoch[i]) &
                        ((iter * minibatch_sizes[i]) %% num_cells[i] != 0)) {
                        # if current iter cycles through the data and start a new cycle
                        epoch_next[i] = TRUE
                        epoch_prev[i] = epoch[i]
                        # shuffle dataset before the next epoch
                        minibatch_idx[[i]] = all_idx[[i]][c(((((iter - 1) * minibatch_sizes[i]) %% num_cells[i]
                        ) + 1):num_cells[i])]
                        chunk_idx[[i]] = sample(1:num_chunks[i], num_chunks[i])
                        all_idx[[i]] = 0
                        for (j in chunk_idx[[i]]) {
                            if (j != num_chunks[i]) {
                                all_idx[[i]] = c(all_idx[[i]],
                                                 (1 + h5_chunk_size * (j - 1)):(j * h5_chunk_size))
                            } else{
                                all_idx[[i]] = c(all_idx[[i]],
                                                 (1 + h5_chunk_size * (j - 1)):num_cells[i])
                            }
                        }
                        all_idx[[i]] = all_idx[[i]][-1] # remove the first element 0
                        minibatch_idx[[i]] = c(minibatch_idx[[i]], all_idx[[i]][1:((iter * minibatch_sizes[i]) %% num_cells[i])])

                    } else if ((epoch_prev[i] != epoch[i]) &
                               ((iter * minibatch_sizes[i]) %% num_cells[i] == 0)) {
                        # if current iter finishes this cycle without start a a new cycle
                        epoch_next[i] = TRUE
                        epoch_prev[i] = epoch[i]

                        minibatch_idx[[i]] = all_idx[[i]][((((iter - 1) * minibatch_sizes[i]
                        ) %% num_cells[i]) + 1):num_cells[i]]
                        chunk_idx[[i]] = sample(1:num_chunks[i], num_chunks[i])
                        all_idx[[i]] = 0
                        for (j in chunk_idx[[i]]) {
                            if (j != num_chunks[i]) {
                                all_idx[[i]] = c(all_idx[[i]],
                                                 (1 + h5_chunk_size * (j - 1)):(j * h5_chunk_size))
                            } else{
                                all_idx[[i]] = c(all_idx[[i]],
                                                 (1 + h5_chunk_size * (j - 1)):num_cells[i])
                            }
                        }
                        all_idx[[i]] = all_idx[[i]][-1] # remove the first element 0
                    } else {
                        # if current iter stays within a single cycle
                        minibatch_idx[[i]] = all_idx[[i]][(((iter - 1) * minibatch_sizes[i]) %% num_cells[i] + 1):((iter * minibatch_sizes[i]) %% num_cells[i])]
                    }
                }
            } else {
                for (i in file_idx_new) {
                    minibatch_sizes[i] = max.epochs * num_cells[i] - (iter - 1) * minibatch_sizes[i]
                    minibatch_idx[[i]] = (((iter - 1) * minibatch_sizes_orig[i] + 1) %% num_cells[i]):num_cells[i]
                }
                epoch[file_idx_new[1]] = max.epochs # last epoch
            }


            if (length(minibatch_idx[[file_idx_new[1]]]) == minibatch_sizes_orig[file_idx_new[1]]) {
                X_minibatch = rep(list(NULL), num_files)
                for (i in file_idx_new) {
                    X_minibatch[[i]] = object@scale.data[[i]][1:num_genes , minibatch_idx[[i]]]
                }

                # update H_i by ANLS Hi_minibatch[[i]]
                H_minibatch = rep(list(NULL), num_files)
                for (i in file_idx_new) {
                    H_minibatch[[i]] = solveNNLS(
                        rbind(
                            object@W + object@V[[i]],
                            sqrt_lambda * object@V[[i]]
                        ),
                        rbind(
                            X_minibatch[[i]],
                            matrix(0, num_genes, minibatch_sizes[i])
                        )
                    )
                }

                # updata A and B matrices
                if (iter == 1) {
                    scale_param = c(rep(0, num_prev_files),
                                    rep(0, num_new_files))
                } else if (iter == 2) {
                    scale_param = c(
                        rep(0, num_prev_files),
                        rep(1, num_new_files) / minibatch_sizes[file_idx_new]
                    )
                } else {
                    scale_param = c(rep(0, num_prev_files),
                                    rep((iter - 2) / (iter - 1), num_new_files))
                }


                if (epoch[file_idx_new[1]] > 0 &
                    epoch_next[file_idx_new[1]] == TRUE) {
                    # remove information older than 2 epochs
                    for (i in file_idx_new) {
                        object@A[[i]] = object@A[[i]] - A_old[[i]]
                        A_old[[i]] = scale_param[i] * object@A[[i]]
                        object@B[[i]] = object@B[[i]] - B_old[[i]]
                        B_old[[i]] = scale_param[i] * object@B[[i]]
                    }
                } else{
                    # otherwise scale the old information
                    for (i in file_idx_new) {
                        A_old[[i]] = scale_param[i] * A_old[[i]]
                        B_old[[i]] = scale_param[i] * B_old[[i]]
                    }
                }

                for (i in file_idx_new) {
                    object@A[[i]] = scale_param[i] * object@A[[i]] + H_minibatch[[i]] %*% t(H_minibatch[[i]]) / minibatch_sizes[i]   # HiHit
                    diag(object@A[[i]])[diag(object@A[[i]]) == 0] = 1e-15
                    object@B[[i]] = scale_param[i] * object@B[[i]] + X_minibatch[[i]] %*% t(H_minibatch[[i]]) / minibatch_sizes[i]   # XiHit
                }


                # update W, V_i by HALS
                iter_miniBatch = 1
                delta_miniBatch = Inf
                max_iters_miniBatch = miniBatch_max_iters

                while (iter_miniBatch <= max_iters_miniBatch) {
                    # update W
                    for (j in 1:k) {
                        W_update_numerator = rep(0, num_genes)
                        W_update_denominator = 0
                        for (i in file_idx) {
                            W_update_numerator = W_update_numerator + object@B[[i]][, j] - (object@W + object@V[[i]]) %*% object@A[[i]][, j]
                            W_update_denominator = W_update_denominator +  object@A[[i]][j, j]
                        }

                        object@W[, j] = nonneg(object@W[, j] + W_update_numerator / W_update_denominator)
                    }

                    # update V_i
                    for (j in 1:k) {
                        for (i in file_idx_new) {
                            object@V[[i]][, j] = nonneg(
                                object@V[[i]][, j] + (
                                    object@B[[i]][, j] - (object@W + (1 + lambda) * object@V[[i]]) %*% object@A[[i]][, j]
                                ) /
                                    ((1 + lambda) * object@A[[i]][j, j])
                            )
                        }
                    }

                    iter_miniBatch = iter_miniBatch + 1
                }
                epoch_next = rep(FALSE, num_files) # reset epoch change indicator
                iter = iter + 1
                if (verbose) {
                    setTxtProgressBar(pb = pb, value = iter)
                }
            }
        }
        if (verbose) {
            cat("\nCalculate metagene loadings...", "\n")
        }
        object@H = rep(list(NULL), num_files)
        for (i in file_idx) {
            if (num_cells[i] %% miniBatch_size == 0)
                num_batch = num_cells[i] %/% miniBatch_size
            else
                num_batch = num_cells[i] %/% miniBatch_size + 1
            if (num_batch == 1) {
                X_i = object@scale.data[[i]][1:num_genes,]
                object@H[[i]] = solveNNLS(
                    rbind(object@W + object@V[[i]], sqrt_lambda * object@V[[i]]),
                    rbind(X_i, matrix(0, num_genes , num_cells[i]))
                )
            } else {
                for (batch_idx in 1:num_batch) {
                    if (batch_idx != num_batch) {
                        cell_idx = ((batch_idx - 1) * miniBatch_size + 1):(batch_idx * miniBatch_size)
                    } else {
                        cell_idx = ((batch_idx - 1) * miniBatch_size + 1):num_cells[i]
                    }
                    X_i_batch = object@scale.data[[i]][1:num_genes, cell_idx]
                    object@H[[i]] = cbind(object@H[[i]], solveNNLS(
                        rbind(
                            object@W + object@V[[i]],
                            sqrt_lambda * object@V[[i]]
                        ),
                        rbind(X_i_batch, matrix(
                            0, num_genes , length(cell_idx)
                        ))
                    ))
                }
            }
            colnames(object@H[[i]]) = cell_barcodes[[i]]
        }

        rownames(object@W) = gene_names
        colnames(object@W) = NULL

        for (i in file_idx) {
            rownames(object@V[[i]]) = gene_names
            colnames(object@V[[i]]) = NULL
        }

    } else {
        if (verbose) {
            cat("Metagene projection", "\n")
        }
        object@W = if (!is.null(W.init))
            W.init
        else
            object@W
        object@H[file_idx_new] = rep(list(NULL), num_new_files)
        object@V[file_idx_new] = rep(list(NULL), num_new_files)
        for (i in file_idx_new) {
            if (num_cells[i] %% miniBatch_size == 0)
                num_batch = num_cells[i] %/% miniBatch_size
            else
                num_batch = num_cells[i] %/% miniBatch_size + 1
            if (num_cells[i] <= miniBatch_size) {
                object@H[[i]] = solveNNLS(object@W, object@scale.data[[i]][1:num_genes,])
            } else {
                for (batch_idx in 1:num_batch) {
                    if (batch_idx != num_batch) {
                        cell_idx = ((batch_idx - 1) * miniBatch_size + 1):(batch_idx * miniBatch_size)
                    } else {
                        cell_idx = ((batch_idx - 1) * miniBatch_size + 1):num_cells[i]
                    }
                    object@H[[i]] = cbind(object@H[[i]],
                                          solveNNLS(object@W, object@scale.data[[i]][1:num_genes, cell_idx]))
                }
            }
            colnames(object@H[[i]]) = cell_barcodes[[i]]
            object@V[[i]] = matrix(0, num_genes, k)
        }
    }

    # gene x k -> k x gene & k x cell -> cell x k
    object@W = t(object@W)
    object@V = lapply(object@V, t)
    object@H = lapply(object@H, t)
    for (i in 1:length(object@raw.data)) {
        if (class(object@raw.data[[i]])[1] != "H5File")
            object@scale.data[[i]] = t(object@scale.data[[i]])
    }

    if (!is.null(X_new)) {
        names(object@scale.data) <-
            names(object@raw.data) <-
            c(names(raw.data_prev), names(X_new))
    }
    names(object@H) <- names(object@V) <- names(object@raw.data)
    return(object)
}


#' Perform thresholding on dense matrix
#'
#' @description
#' Perform thresholding on the input dense matrix. Remove any values samller than eps by eps.
#' Helper function for online_iNMF
#'
#' @param x Dense matrix.
#' @param eps Threshold. Should be a small positive value. (default 1e-16)
#' @return Dense matrix with smallest values equal to eps.

nonneg <- function(x, eps = 1e-16) {
    x[x < eps] = eps
    return(x)
}


#' Perform iNMF on scaled datasets
#'
#' @description
#' Perform integrative non-negative matrix factorization to return factorized H, W, and V matrices.
#' It optimizes the iNMF objective function using block coordinate descent (alternating non-negative
#' least squares), where the number of factors is set by k. TODO: include objective function
#' equation here in documentation (using deqn)
#'
#' For each dataset, this factorization produces an H matrix (cells by k), a V matrix (k by genes),
#' and a shared W matrix (k by genes). The H matrices represent the cell factor loadings.
#' W is held consistent among all datasets, as it represents the shared components of the metagenes
#' across datasets. The V matrices represent the dataset-specific components of the metagenes.
#'
#' @param object \code{liger} object. Should normalize, select genes, and scale before calling.
#' @param k Inner dimension of factorization (number of factors). Run suggestK to determine
#'   appropriate value; a general rule of thumb is that a higher k will be needed for datasets with
#'   more sub-structure.
#' @param lambda Regularization parameter. Larger values penalize dataset-specific effects more
#'   strongly (ie. alignment should increase as lambda increases). Run suggestLambda to determine
#'   most appropriate value for balancing dataset alignment and agreement (default 5.0).
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh.
#'   (default 1e-6)
#' @param max.iters Maximum number of block coordinate descent iterations to perform (default 30).
#' @param nrep Number of restarts to perform (iNMF objective function is non-convex, so taking the
#'   best objective from multiple successive initializations is recommended). For easier
#'   reproducibility, this increments the random seed by 1 for each consecutive restart, so future
#'   factorizations of the same dataset can be run with one rep if necessary. (default 1)
#' @param H.init Initial values to use for H matrices. (default NULL)
#' @param W.init Initial values to use for W matrix (default NULL)
#' @param V.init Initial values to use for V matrices (default NULL)
#' @param rand.seed Random seed to allow reproducible results (default 1).
#' @param print.obj Print objective function values after convergence (default FALSE).
#' @param verbose Print progress bar/messages (TRUE by default)
#' @param ... Arguments passed to other methods
#'
#' @return \code{liger} object with H, W, and V slots set.
#'
#' @export
#' @examples
#' \dontrun{
#' # Requires preprocessed liger object (only for objected not based on HDF5 files)
#' # Get factorization using 20 factors and mini-batch of 5000 cells
#' # (default setting, can be adjusted for ideal results)
#' ligerex <- optimizeALS(ligerex, k = 20, lambda = 5, nrep = 1)
#' }

optimizeALS <- function(object,
                        ...) {
    UseMethod(generic = 'optimizeALS', object = object)
}

#' @rdname optimizeALS
#' @importFrom stats runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
#' @method optimizeALS list
#'
optimizeALS.list <- function(object,
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
                             verbose = TRUE,
                             ...) {
    if (!all(sapply(X = object, FUN = is.matrix))) {
        stop("All values in 'object' must be a matrix")
    }
    E <- object
    N <- length(x = E)
    ns <- sapply(X = E, FUN = nrow)
    #if (k >= min(ns)) {
    #  stop('Select k lower than the number of cells in smallest dataset: ', min(ns))
    #}
    tmp <- gc()
    g <- ncol(x = E[[1]])
    if (k >= g) {
        stop('Select k lower than the number of variable genes: ', g)
    }
    W_m <- matrix(data = 0,
                  nrow = k,
                  ncol = g)
    V_m <- lapply(
        X = 1:N,
        FUN = function(i) {
            return(matrix(
                data = 0,
                nrow = k,
                ncol = g
            ))
        }
    )
    H_m <- lapply(
        X = ns,
        FUN = function(n) {
            return(matrix(
                data = 0,
                nrow = n,
                ncol = k
            ))
        }
    )
    tmp <- gc()
    best_obj <- Inf
    run_stats <- matrix(data = 0,
                        nrow = nrep,
                        ncol = 2)
    for (i in 1:nrep) {
        set.seed(seed = rand.seed + i - 1)
        start_time <- Sys.time()
        W <- matrix(data = abs(x = runif(
            n = g * k,
            min = 0,
            max = 2
        )),
        nrow = k,
        ncol = g)

        V <- lapply(
            X = 1:N,
            FUN = function(i) {
                return(matrix(
                    data = abs(x = runif(
                        n = g * k,
                        min = 0,
                        max = 2
                    )),
                    nrow = k,
                    ncol = g
                ))
            }
        )

        H <- lapply(
            X = ns,
            FUN = function(n) {
                return(matrix(
                    data = abs(x = runif(
                        n = n * k,
                        min = 0,
                        max = 2
                    )),
                    nrow = n,
                    ncol = k
                ))
            }
        )
        tmp <- gc()
        if (!is.null(x = W.init)) {
            W <- W.init
        }
        if (!is.null(x = V.init)) {
            V <- V.init
        }
        if (!is.null(x = H.init)) {
            H <- H.init
        }
        delta <- 1
        iters <- 0
        pb <- txtProgressBar(min = 0,
                             max = max.iters,
                             style = 3)
        sqrt_lambda <- sqrt(x = lambda)
        obj0 <- sum(sapply(
            X = 1:N,
            FUN = function(i) {
                return(norm(x = E[[i]] - H[[i]] %*% (W + V[[i]]), type = "F") ^ 2)
            }
        )) +
            sum(sapply(
                X = 1:N,
                FUN = function(i) {
                    return(lambda * norm(x = H[[i]] %*% V[[i]], type = "F") ^ 2)
                }
            ))
        tmp <- gc()
        while (delta > thresh & iters < max.iters) {
            H <- lapply(
                X = 1:N,
                FUN = function(i) {
                    return(t(x = solveNNLS(
                        C = rbind(t(x = W) + t(x = V[[i]]), sqrt_lambda * t(x = V[[i]])),
                        B = rbind(t(x = E[[i]]), matrix(
                            data = 0,
                            nrow = g,
                            ncol = ns[i]
                        ))
                    )))
                }
            )
            tmp <- gc()
            V <- lapply(
                X = 1:N,
                FUN = function(i) {
                    return(solveNNLS(
                        C = rbind(H[[i]], sqrt_lambda * H[[i]]),
                        B = rbind(E[[i]] - H[[i]] %*% W, matrix(
                            data = 0,
                            nrow = ns[[i]],
                            ncol = g
                        ))
                    ))
                }
            )
            tmp <- gc()
            W <- solveNNLS(C = rbindlist(mat_list = H),
                           B = rbindlist(mat_list = lapply(
                               X = 1:N,
                               FUN = function(i) {
                                   return(E[[i]] - H[[i]] %*% V[[i]])
                               }
                           )))
            tmp <- gc()
            obj <- sum(sapply(
                X = 1:N,
                FUN = function(i) {
                    return(norm(x = E[[i]] - H[[i]] %*% (W + V[[i]]), type = "F") ^ 2)
                }
            )) +
                sum(sapply(
                    X = 1:N,
                    FUN = function(i) {
                        return(lambda * norm(x = H[[i]] %*% V[[i]], type = "F") ^ 2)
                    }
                ))
            tmp <- gc()
            delta <- abs(x = obj0 - obj) / (mean(obj0, obj))
            obj0 <- obj
            iters <- iters + 1
            setTxtProgressBar(pb = pb, value = iters)
        }
        setTxtProgressBar(pb = pb, value = max.iters)
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
        end_time <-
            difftime(time1 = Sys.time(),
                     time2 = start_time,
                     units = "auto")
        run_stats[i, 1] <- as.double(x = end_time)
        run_stats[i, 2] <- iters
        if (verbose) {
            cat(
                "\nFinished in ",
                run_stats[i, 1],
                " ",
                units(x = end_time),
                ", ",
                iters,
                " iterations.\n",
                "Max iterations set: ",
                max.iters,
                ".\n",
                "Final objective delta: ",
                delta,
                '.\n',
                sep = ""
            )
        }

        if (verbose) {
            if (print.obj) {
                cat("Objective:", obj, "\n")
            }
            cat("Best results with seed ", best_seed, ".\n", sep = "")
        }
    }
    out <- list()
    out$H <- H_m
    for (i in 1:length(x = object)) {
        rownames(x = out$H[[i]]) <- rownames(x = object[[i]])
    }
    out$V <- V_m
    names(x = out$V) <- names(x = out$H) <- names(x = object)
    out$W <- W_m
    return(out)
}

#' @importFrom methods slot<-
#'
#' @rdname optimizeALS
#' @export
#' @method optimizeALS liger
#'
optimizeALS.liger <- function(object,
                              k,
                              lambda = 5.0,
                              thresh = 1e-6,
                              max.iters = 30,
                              nrep = 1,
                              H.init = NULL,
                              W.init = NULL,
                              V.init = NULL,
                              use.unshared = FALSE,
                              rand.seed = 1,
                              print.obj = FALSE,
                              verbose = TRUE,
                              ...) {
    if (use.unshared == FALSE) {
        object <- removeMissingObs(
            object = object,
            slot.use = 'scale.data',
            use.cols = FALSE,
            verbose = TRUE
        )
        out <- optimizeALS(
            object = object@scale.data,
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
        names(x = out$H) <-
            names(x = out$V) <- names(x = object@raw.data)
        for (i in 1:length(x = object@scale.data)) {
            rownames(x = out$H[[i]]) <- rownames(x = object@scale.data[[i]])
        }
        colnames(x = out$W) <- object@var.genes
        for (i in names(x = out)) {
            slot(object = object, name = i) <- out[[i]]
        }
        object@parameters$lambda <- lambda
        return(object)
    }
    if (use.unshared == TRUE) {
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

#' Perform factorization for new value of k
#'
#' This uses an efficient strategy for updating that takes advantage of the information in the
#' existing factorization. It is most recommended for values of k smaller than current value,
#' where it is more likely to speed up the factorization.
#'
#' @param object \code{liger} object. Should call optimizeALS before calling.
#' @param k.new Inner dimension of factorization (number of factors)
#' @param lambda Regularization parameter. By default, this will use the lambda last used with
#'   optimizeALS.
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#'   (default 1e-4).
#' @param max.iters Maximum number of block coordinate descent iterations to perform (default 100).
#' @param rand.seed Random seed to set. Only relevant if k.new > k. (default 1)
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return \code{liger} object with H, W, and V slots reset.
#'
#' @importFrom plyr rbind.fill.matrix
#'
#' @export
#' @examples
#' \dontrun{
#' # decide to run with k = 15 instead (keeping old lambda the same)
#' ligerex <- optimizeNewK(ligerex, k.new = 15)
#' }

optimizeNewK <-
    function(object,
             k.new,
             lambda = NULL,
             thresh = 1e-4,
             max.iters = 100,
             rand.seed = 1,
             verbose = TRUE) {
        if (is.null(lambda)) {
            lambda <- object@parameters$lambda
        }
        k <- ncol(object@H[[1]])
        if (k.new == k) {
            return(object)
        }
        H <- object@H
        W <- object@W
        V <- object@V

        if (k.new > k) {
            set.seed(rand.seed)
            sqrt_lambda <- sqrt(lambda)
            g <- ncol(W)
            N <- length(H)
            ns <- sapply(H, nrow)
            W_new <- matrix(abs(runif(g * k, 0, 2)), k.new - k, g)
            V_new <- lapply(1:N, function(i) {
                matrix(abs(runif(g * (
                    k.new - k
                ), 0, 2)), k.new - k, g)
            })
            H_new <- lapply(ns, function(n) {
                matrix(abs(runif(n * (
                    k.new - k
                ), 0, 2)), n, k.new - k)
            })
            H_new <- lapply(1:N, function(i) {
                t(solveNNLS(rbind(
                    t(W_new) + t(V_new[[i]]), sqrt_lambda * t(V_new[[i]])
                ),
                rbind(
                    t(object@scale.data[[i]] - H[[i]] %*% (W + V[[i]])),
                    matrix(0, nrow = g, ncol = ns[i])
                )))
            })
            V_new <- lapply(1:N, function(i) {
                solveNNLS(
                    rbind(H_new[[i]], sqrt_lambda * H_new[[i]]),
                    rbind(
                        object@scale.data[[i]] - H[[i]] %*% (W + V[[i]]) - H_new[[i]] %*% W_new,
                        matrix(0, nrow = ns[[i]], ncol = g)
                    )
                )
            })
            W_new <- solveNNLS(rbind.fill.matrix(H_new),
                               rbind.fill.matrix(lapply(1:N, function(i) {
                                   object@scale.data[[i]] - H[[i]] %*% (W + V[[i]]) - H_new[[i]] %*% V_new[[i]]
                               })))
            H <- lapply(1:N, function(i) {
                cbind(H[[i]], H_new[[i]])
            })
            V <- lapply(1:N, function(i) {
                rbind(V[[i]], V_new[[i]])
            })
            W <- rbind(W, W_new)
        }
        else {
            deltas <- rep(0, k)
            for (i in 1:length(object@H))
            {
                deltas <- deltas + sapply(1:k, function(x) {
                    norm(H[[i]][, k] %*% t(W[k, ] + V[[i]][k, ]), "F")
                })
            }
            k.use <- order(deltas, decreasing = TRUE)[1:k.new]
            W <- W[k.use, ]
            H <- lapply(H, function(x) {
                x[, k.use]
            })
            V <- lapply(V, function(x) {
                x[k.use, ]
            })
        }
        object <- optimizeALS(
            object,
            k.new,
            lambda = lambda,
            thresh = thresh,
            max.iters = max.iters,
            H.init = H,
            W.init = W,
            V.init = V,
            rand.seed = rand.seed,
            verbose = verbose
        )
        return(object)
    }

#' Perform factorization for new data
#'
#' Uses an efficient strategy for updating that takes advantage of the information in the existing
#' factorization. Assumes that selected genes (var.genes) are represented in the new datasets.
#'
#' @param object \code{liger} object. Should call optimizeALS before calling.
#' @param new.data List of raw data matrices (one or more). Each list entry should be named.
#' @param which.datasets List of datasets to append new.data to if add.to.existing is true.
#'   Otherwise, the most similar existing datasets for each entry in new.data.
#' @param add.to.existing Add the new data to existing datasets or treat as totally new datasets
#'   (calculate new Vs?) (default TRUE)
#' @param lambda Regularization parameter. By default, this will use the lambda last used with
#'   optimizeALS.
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#'   (default 1e-4).
#' @param max.iters Maximum number of block coordinate descent iterations to perform (default 100).
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return \code{liger} object with H, W, and V slots reset. Raw.data, norm.data, and scale.data will
#'   also be updated to include the new data.
#'
#' @export
#' @examples
#' \dontrun{
#' # Given preprocessed liger object: ligerex (contains two datasets Y and Z)
#' # get factorization using three restarts and 20 factors
#' ligerex <- optimizeALS(ligerex, k = 20, lambda = 5, nrep = 3)
#' # acquire new data (Y_new, Z_new) from the same cell type, let's add it to existing datasets
#' new_data <- list(Y_set = Y_new, Z_set = Z_new)
#' ligerex2 <- optimizeNewData(ligerex, new.data = new_data, which.datasets = list('y_set', 'z_set'))
#' # acquire new data from different cell type (X), we'll just add another dataset
#' # it's probably most similar to y_set
#' ligerex <- optimizeNewData(ligerex, new.data = list(x_set = X), which.datasets = list('y_set'),
#'                            add.to.existing = FALSE)
#' }

optimizeNewData <-
    function(object,
             new.data,
             which.datasets,
             add.to.existing = TRUE,
             lambda = NULL,
             thresh = 1e-4,
             max.iters = 100,
             verbose = TRUE) {
        if (is.null(lambda)) {
            lambda <- object@parameters$lambda
        }
        if (add.to.existing) {
            for (i in 1:length(new.data)) {
                if (verbose) {
                    message(dim(object@raw.data[[which.datasets[[i]]]]))
                }
                object@raw.data[[which.datasets[[i]]]] <-
                    cbind(object@raw.data[[which.datasets[[i]]]],
                          new.data[[i]])
                if (verbose) {
                    message(dim(object@raw.data[[which.datasets[[i]]]]))
                }
            }
            object <- normalize(object)
            object <- scaleNotCenter(object)
            sqrt_lambda <- sqrt(lambda)
            g <- ncol(object@W)
            H_new <- lapply(1:length(new.data), function(i) {
                t(solveNNLS(rbind(
                    t(object@W) + t(object@V[[which.datasets[[i]]]]),
                    sqrt_lambda * t(object@V[[which.datasets[[i]]]])
                ),
                rbind(
                    t(object@scale.data[[which.datasets[[i]]]][colnames(new.data[[i]]), ]),
                    matrix(0, nrow = g, ncol = ncol(new.data[[i]]))
                )))
            })
            for (i in 1:length(new.data)) {
                object@H[[which.datasets[[i]]]] <-
                    rbind(object@H[[which.datasets[[i]]]], H_new[[i]])
            }
        } else {
            old.names <- names(object@raw.data)
            new.names <- names(new.data)
            combined.names <- c(old.names, new.names)
            for (i in 1:length(which.datasets)) {
                object@V[[names(new.data)[i]]] <- object@V[[which.datasets[[i]]]]
            }
            object@raw.data <- c(object@raw.data, new.data)
            names(object@raw.data) <-
                names(object@V) <- combined.names
            object <- normalize(object)
            object <- scaleNotCenter(object)
            ns <- lapply(object@raw.data, ncol)
            N <- length(ns)
            g <- ncol(object@W)
            sqrt_lambda <- sqrt(lambda)
            for (i in 1:N) {
                print(ns[[i]])
                print(dim(object@raw.data[[i]]))
                print(dim(object@norm.data[[i]]))
                print(dim(object@scale.data[[i]]))
                print(dim(object@V[[i]]))
            }
            H_new <- lapply(1:length(new.data), function(i) {
                t(solveNNLS(rbind(
                    t(object@W) + t(object@V[[new.names[i]]]),
                    sqrt_lambda * t(object@V[[new.names[i]]])
                ),
                rbind(
                    t(object@scale.data[[new.names[i]]]),
                    matrix(0, nrow = g, ncol = ncol(new.data[[i]]))
                )))
            })
            object@H <- c(object@H, H_new)
            names(object@H) <- combined.names
        }
        k <- ncol(object@H[[1]])
        object <- optimizeALS(
            object,
            k,
            lambda,
            thresh,
            max.iters,
            H.init = object@H,
            W.init = object@W,
            V.init = object@V,
            verbose = verbose
        )
        return(object)
    }

#' Perform factorization for subset of data
#'
#' Uses an efficient strategy for updating that takes advantage of the information in the existing
#' factorization. Can use either cell names or cluster names to subset. For more basic subsetting
#' functionality (without automatic optimization), see subsetLiger.
#'
#' @param object \code{liger} object. Should call optimizeALS before calling.
#' @param cell.subset List of cell names to retain from each dataset (same length as number of
#'   datasets).
#' @param cluster.subset Clusters for which to keep cells (ie. c(1, 5, 6)). Should pass in either
#'   cell.subset or cluster.subset but not both.
#' @param lambda Regularization parameter. By default, uses last used lambda.
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#'   (default 1e-4).
#' @param max.iters Maximum number of block coordinate descent iterations to perform (default 100).
#' @param datasets.scale Names of datasets to rescale after subsetting (default NULL).
#'
#' @return \code{liger} object with H, W, and V slots reset. Scale.data
#'   (if desired) will also be updated to reflect the subset.
#'
#' @export
#' @examples
#' \dontrun{
#' # now want to look at only subset of data
#' # Requires a vector of cell names from data 1 and a vector of cell names from data 2
#' ligerex2 <- optimizeSubset(ligerex, cell.subset = list(cell_names_1, cell_names_2))
#' }

optimizeSubset <-
    function(object,
             cell.subset = NULL,
             cluster.subset = NULL,
             lambda = NULL,
             thresh = 1e-4,
             max.iters = 100,
             datasets.scale = NULL) {
        if (is.null(lambda)) {
            lambda <- object@parameters$lambda
        }
        if (is.null(cell.subset) & is.null(cluster.subset)) {
            stop("Please specify a cell subset or cluster subset.")
        }
        else if (is.null(cell.subset) & !is.null(cluster.subset)) {
            cell.subset <- lapply(1:length(object@scale.data), function(i) {
                which(object@clusters[rownames(object@scale.data[[i]])] %in% cluster.subset)
            })
        }
        old_names <- names(object@raw.data)
        H <- object@H
        H <- lapply(1:length(object@H), function(i) {
            object@H[[i]][cell.subset[[i]], ]
        })
        object@raw.data <-
            lapply(1:length(object@raw.data), function(i) {
                object@raw.data[[i]][, cell.subset[[i]]]
            })
        object@cell.data <-
            droplevels(object@cell.data[cell.subset, ])
        for (i in 1:length(object@norm.data)) {
            object@norm.data[[i]] <- object@norm.data[[i]][, cell.subset[[i]]]
            if (names(object@norm.data)[i] %in% datasets.scale) {
                object@scale.data[[i]] <-
                    scale(t(object@norm.data[[i]][object@var.genes, ]),
                          scale = TRUE,
                          center = FALSE)
                object@scale.data[[i]][is.na(object@scale.data[[i]])] <-
                    0
            } else {
                object@scale.data[[i]] <-
                    t(object@norm.data[[i]][object@var.genes, ])
            }
        }

        names(object@raw.data) <-
            names(object@norm.data) <- names(object@H) <- old_names
        k <- ncol(H[[1]])
        object <-
            optimizeALS(
                object,
                k = k,
                lambda = lambda,
                thresh = thresh,
                max.iters = max.iters,
                H.init = H,
                W.init = object@W,
                V.init = object@V
            )
        return(object)
    }

#' Perform factorization for new lambda value
#'
#' Uses an efficient strategy for updating that takes advantage of the information in the existing
#' factorization; uses previous k. Recommended mainly when re-optimizing for higher lambda and when
#' new lambda value is significantly different; otherwise may not return optimal results.
#'
#' @param object \code{liger} object. Should call optimizeALS before calling.
#' @param new.lambda Regularization parameter. Larger values penalize dataset-specific effects more
#' strongly.
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max.iters Maximum number of block coordinate descent iterations to perform (default 100).
#' @param rand.seed Random seed for reproducibility (default 1).
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return \code{liger} object with optimized factorization values
#'
#' @export
#' @examples
#' \dontrun{
#' # decide to run with lambda = 15 instead (keeping k the same)
#' ligerex <- optimizeNewLambda(ligerex, new.lambda = 15)
#' }

optimizeNewLambda <-
    function(object,
             new.lambda,
             thresh = 1e-4,
             max.iters = 100,
             rand.seed = 1,
             verbose = TRUE) {
        k <- ncol(object@H[[1]])
        H <- object@H
        W <- object@W
        if (new.lambda < object@parameters$lambda && verbose) {
            message(
                "New lambda less than current lambda; new factorization may not be optimal. ",
                "Re-optimization with optimizeAlS recommended instead."
            )
        }
        object <-
            optimizeALS(
                object,
                k,
                lambda = new.lambda,
                thresh = thresh,
                max.iters = max.iters,
                H.init = H,
                W.init = W,
                rand.seed = rand.seed,
                verbose = verbose
            )
        return(object)
    }

#' Visually suggest appropriate lambda value
#'
#' Can be used to select appropriate value of lambda for factorization of particular dataset. Plot
#' alignment and agreement for various test values of lambda. Most appropriate lambda
#' is likely around the "elbow" of the alignment plot (when alignment stops increasing). This will
#' likely also correspond to slower decrease in agreement. Depending on number of cores used,
#' this process can take 10-20 minutes.
#'
#' @param object \code{liger} object. Should normalize, select genes, and scale before calling.
#' @param k Number of factors to use in test factorizations. See optimizeALS documentation.
#' @param lambda.test Vector of lambda values to test. If not given, use default set spanning
#'   0.25 to 60
#' @param rand.seed Random seed for reproducibility (default 1).
#' @param num.cores Number of cores to use for optimizing factorizations in parallel (default 1).
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max.iters Maximum number of block coordinate descent iterations to perform
#' @param knn_k Number of nearest neighbors for within-dataset knn in quantileAlignSNF (default 20).
#' @param k2 Horizon parameter for quantileAlignSNF (default 500).
#' @param ref_dataset Reference dataset for quantileAlignSNF (defaults to larger dataset).
#' @param resolution Resolution for quantileAlignSNF (default 1).
#' @param gen.new Do not use optimizeNewLambda in factorizations. Recommended to set TRUE
#'   when looking at only a small range of lambdas (ie. 1:7) (default FALSE)
#' @param nrep Number restarts to perform at each lambda value tested (increase to produce
#'   smoother curve if results unclear) (default 1).
#' @param return.data Whether to return list of data matrices (raw) or dataframe (processed)
#'   instead of ggplot object (default FALSE).
#' @param return.raw If return.results TRUE, whether to return raw data (in format described below),
#'   or dataframe used to produce ggplot object. Raw data is matrix of alignment values for each
#'   lambda value tested (each column represents a different rep for nrep).(default FALSE)
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return Matrix of results if indicated or ggplot object. Plots alignment vs. lambda to console.
#'
#' @import doParallel
#' @import parallel
#' @importFrom foreach foreach
#' @importFrom foreach "%dopar%"
#' @importFrom ggplot2 ggplot aes geom_point geom_line guides guide_legend labs theme theme_classic
#'
#' @export
#' @examples
#' \dontrun{
#' # Requires preprocessed liger object
#' # examine plot for most appropriate lambda, use multiple cores for faster results
#' suggestLambda(ligerex, k = 20, num.cores = 4)
#' }

suggestLambda <-
    function(object,
             k,
             lambda.test = NULL,
             rand.seed = 1,
             num.cores = 1,
             thresh = 1e-4,
             max.iters = 100,
             knn_k = 20,
             k2 = 500,
             ref_dataset = NULL,
             resolution = 1,
             gen.new = FALSE,
             nrep = 1,
             return.data = FALSE,
             return.raw = FALSE,
             verbose = TRUE) {
        if (is.null(lambda.test)) {
            lambda.test <- c(seq(0.25, 1, 0.25), seq(2, 10, 1), seq(15, 60, 5))
        }
        time_start <- Sys.time()
        # optimize smallest lambda value first to take advantage of efficient updating
        if (verbose) {
            message(
                "This operation may take several minutes depending on number of values being tested"
            )
        }
        rep_data <- list()
        for (r in 1:nrep) {
            if (verbose) {
                message(
                    "Preprocessing for rep ",
                    r,
                    ": optimizing initial factorization with smallest test lambda=",
                    lambda.test[1]
                )
            }
            object <-
                optimizeALS(
                    object,
                    k = k,
                    thresh = thresh,
                    lambda = lambda.test[1],
                    max.iters = max.iters,
                    nrep = 1,
                    rand.seed = (rand.seed + r - 1),
                    verbose = verbose
                )
            if (verbose) {
                message('Testing different choices of lambda values')
            }
            #cl <- makeCluster(num.cores)
            cl <- parallel::makeCluster(num.cores)
            #registerDoSNOW(cl)
            doParallel::registerDoParallel(cl)
            #pb <- txtProgressBar(min = 0, max = length(lambda.test), style = 3, initial = 1, file = "")
            # define progress bar function
            #progress <- function(n) setTxtProgressBar(pb, n)
            #opts <- list(progress = progress)
            i <- 0
            data_matrix <-
                foreach(
                    i = 1:length(lambda.test),
                    .combine = "rbind",
                    .packages = 'rliger'
                ) %dopar% {
                    if (i != 1) {
                        if (gen.new) {
                            ob.test <- optimizeALS(
                                object,
                                k = k,
                                lambda = lambda.test[i],
                                thresh = thresh,
                                max.iters = max.iters,
                                rand.seed = (rand.seed + r - 1)
                            )
                        } else {
                            ob.test <- optimizeNewLambda(
                                object,
                                new.lambda = lambda.test[i],
                                thresh = thresh,
                                max.iters = max.iters,
                                rand.seed = (rand.seed + r - 1)
                            )
                        }
                    } else {
                        ob.test <- object
                    }
                    ob.test <-
                        quantileAlignSNF(
                            ob.test,
                            knn_k = knn_k,
                            k2 = k2,
                            resolution = resolution,
                            ref_dataset = ref_dataset,
                            id.number = i
                        )
                    calcAlignment(ob.test)
                }
            #close(pb)
            parallel::stopCluster(cl)
            rep_data[[r]] <- data_matrix
        }

        aligns <- Reduce(cbind, rep_data)
        if (is.null(dim(aligns))) {
            aligns <- matrix(aligns, ncol = 1)
        }
        mean_aligns <- apply(aligns, 1, mean)

        time_elapsed <-
            difftime(Sys.time(), time_start, units = "auto")
        if (verbose) {
            cat(paste(
                "\nCompleted in:",
                as.double(time_elapsed),
                units(time_elapsed)
            ))
        }
        # make dataframe
        df_al <-
            data.frame(align = mean_aligns, lambda = lambda.test)

        p1 <-
            ggplot(df_al, aes_string(x = 'lambda', y = 'mean_aligns')) + geom_line(size =
                                                                                       1) +
            geom_point() +
            theme_classic() + labs(y = 'Alignment', x = 'Lambda') +
            guides(col = guide_legend(title = "", override.aes = list(size = 2))) +
            theme(legend.position = 'top')

        if (return.data) {
            print(p1)
            if (return.raw) {
                rownames(aligns) <- lambda.test
                return(aligns)
            }
            return(df_al)
        }
        return(p1)
    }

#' Visually suggest appropiate k value
#'
#' @description
#' This can be used to select appropriate value of k for factorization of particular dataset.
#' Plots median (across cells in all datasets) K-L divergence from uniform for cell factor loadings
#' as a function of k. This should increase as k increases but is expected to level off above
#' sufficiently high number of factors (k). This is because cells should have factor loadings which
#' are not uniformly distributed when an appropriate number of factors is reached.
#'
#' Depending on number of cores used, this process can take 10-20 minutes.
#'
#' @param object \code{liger} object. Should normalize, select genes, and scale before calling.
#' @param k.test Set of factor numbers to test (default seq(5, 50, 5)).
#' @param lambda Lambda to use for all foctorizations (default 5).
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max.iters Maximum number of block coordinate descent iterations to perform
#' @param num.cores Number of cores to use for optimizing factorizations in parallel (default 1)
#' @param rand.seed Random seed for reproducibility (default 1).
#' @param gen.new Do not use optimizeNewK in factorizations. Results in slower factorizations.
#'   (default FALSE).
#' @param nrep Number restarts to perform at each k value tested (increase to produce
#'   smoother curve if results unclear) (default 1).
#' @param plot.log2 Plot log2 curve for reference on K-L plot (log2 is upper bound and con
#'   sometimes help in identifying "elbow" of plot). (default TRUE)
#' @param return.data Whether to return list of data matrices (raw) or dataframe (processed)
#'   instead of ggplot object (default FALSE).
#' @param return.raw If return.results TRUE, whether to return raw data (in format described below),
#'   or dataframe used to produce ggplot object. Raw data is list of matrices of K-L divergences
#'   (length(k.test) by n_cells). Length of list corresponds to nrep. (default FALSE)
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return Matrix of results if indicated or ggplot object. Plots K-L divergence vs. k to console.
#'
#' @import doParallel
#' @import parallel
#' @importFrom foreach foreach
#' @importFrom foreach "%dopar%"
#' @importFrom ggplot2 ggplot aes geom_point geom_line guides guide_legend labs theme theme_classic
#'
#' @export
#' @examples
#' \dontrun{
#' # Requires preprocessed liger object
#' # examine plot for most appropriate k, use multiple cores for faster results
#' suggestK(ligerex, num.cores = 4)
#' }

suggestK <-
    function(object,
             k.test = seq(5, 50, 5),
             lambda = 5,
             thresh = 1e-4,
             max.iters = 100,
             num.cores = 1,
             rand.seed = 1,
             gen.new = FALSE,
             nrep = 1,
             plot.log2 = TRUE,
             return.data = FALSE,
             return.raw = FALSE,
             verbose = TRUE) {
        if (length(object@scale.data) == 0) {
            stop("scaleNotCenter should be run on the object before running suggestK.")
        }
        time_start <- Sys.time()
        # optimize largest k value first to take advantage of efficient updating
        if (verbose) {
            message(
                "This operation may take several minutes depending on number of values being tested"
            )
        }
        rep_data <- list()
        for (r in 1:nrep) {
            if (verbose) {
                message(
                    "Preprocessing for rep ",
                    r,
                    ": optimizing initial factorization with largest test k=",
                    k.test[length(k.test)]
                )
            }
            object <-
                optimizeALS(
                    object,
                    k = k.test[length(k.test)],
                    lambda = lambda,
                    thresh = thresh,
                    max.iters = max.iters,
                    nrep = 1,
                    rand.seed = (rand.seed + r - 1)
                )
            if (verbose) {
                message('Testing different choices of k')
            }
            cl <- parallel::makeCluster(num.cores)
            doParallel::registerDoParallel(cl)
            #pb <- txtProgressBar(min = 0, max = length(k.test), style = 3, initial = 1, file = "")
            # define progress bar function
            #progress <- function(n) setTxtProgressBar(pb, n)
            #opts <- list(progress = progress)
            i <- 0
            data_matrix <-
                foreach(
                    i = length(k.test):1,
                    .combine = "rbind",
                    .packages = 'rliger'
                ) %dopar% {
                    if (i != length(k.test)) {
                        if (gen.new) {
                            ob.test <- optimizeALS(
                                object,
                                k = k.test[i],
                                lambda = lambda,
                                thresh = thresh,
                                max.iters = max.iters,
                                rand.seed = (rand.seed + r - 1)
                            )
                        } else {
                            ob.test <- optimizeNewK(
                                object,
                                k.new = k.test[i],
                                lambda = lambda,
                                thresh = thresh,
                                max.iters = max.iters,
                                rand.seed = (rand.seed + r - 1)
                            )
                        }
                    } else {
                        ob.test <- object
                    }
                    dataset_split <-
                        kl_divergence_uniform(ob.test)
                    unlist(dataset_split)
                }
            #close(pb)
            parallel::stopCluster(cl)
            data_matrix <- data_matrix[nrow(data_matrix):1, ]
            rep_data[[r]] <- data_matrix
        }

        medians <-
            Reduce(cbind, lapply(rep_data, function(x) {
                apply(x, 1, median)
            }))
        if (is.null(dim(medians))) {
            medians <- matrix(medians, ncol = 1)
        }
        mean_kls <- apply(medians, 1, mean)

        time_elapsed <-
            difftime(Sys.time(), time_start, units = "auto")
        if (verbose) {
            cat(paste(
                "\nCompleted in:",
                as.double(time_elapsed),
                units(time_elapsed)
            ))
        }
        # make dataframe
        df_kl <-
            data.frame(
                median_kl = c(mean_kls, log2(k.test)),
                k = c(k.test, k.test),
                calc = c(rep('KL_div', length(k.test)), rep('log2(k)', length(k.test)))
            )
        if (!plot.log2) {
            df_kl <- df_kl[df_kl$calc == 'KL_div', ]
        }

        p1 <-
            ggplot(df_kl, aes_string(x = 'k', y = 'median_kl', col = 'calc')) + geom_line(size =
                                                                                              1) +
            geom_point() +
            theme_classic() + labs(y = 'Median KL divergence (across all cells)', x = 'K') +
            guides(col = guide_legend(title = "", override.aes = list(size = 2))) +
            theme(legend.position = 'top')

        if (return.data) {
            print(p1)
            if (return.raw) {
                rep_data <- lapply(rep_data, function(x) {
                    rownames(x) <- k.test
                    return(x)
                })
                return(rep_data)
            }
            return(df_kl)
        }
        return(p1)
    }

# Conversion/Transformation ####

#' Perform iNMF on scaled datasets, and include unshared, scaled and normalized, features
#' @param object \code{liger} object. Should normalize, select genes, and scale before calling.
#' @param k Inner dimension of factorization (number of factors).
#' @param lambda The lambda penalty. Default 5
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh.
#'   (default 1e-6)
#' @param max.iters Maximum number of block coordinate descent iterations to perform (default 30).
#' @param nrep Number of restarts to perform (iNMF objective function is non-convex, so taking the
#'   best objective from multiple successive initializations is recommended). For easier
#'   reproducibility, this increments the random seed by 1 for each consecutive restart, so future
#'   factorizations of the same dataset can be run with one rep if necessary. (default 1)
#' @param rand.seed Random seed to allow reproducible results (default 1).
#' @param print.obj  Print objective function values after convergence (default FALSE).
#' @param vectorized.lamba Whether or not to expect a vectorized lambda parameter
#' ##########################################################################
optimize_UANLS = function(object,
                          k = 30,
                          lambda = 5,
                          max.iters = 30,
                          nrep = 1,
                          thresh = 1e-10,
                          rand.seed = 1,
                          print.obj = FALSE,
                          vectorized.lambda = FALSE) {
    set.seed(seed = rand.seed)
    #Account for vectorized lambda
    print('Performing Factorization using UINMF and unshared features')
    if (vectorized.lambda == FALSE) {
        lambda = rep(lambda, length(names(object@raw.data)))
    }

    # Get a list of all the matrices
    mlist = list()
    xdim =  list()
    for (i in 1:length(object@scale.data)) {
        mlist[[i]] = t(object@scale.data[[i]])
        xdim[[i]] = dim(mlist[[i]])
    }

    #return what datasets have unshared features, and the dimensions of those unshared features
    u_dim <- c()
    max_feats = 0
    unshared <- c()
    ulist <- c()
    for (i in 1:length(object@var.unshared.features)) {
        if (length(object@var.unshared.features[[i]])) {
            u_dim[[i]] <- dim(object@scale.unshared.data[[i]])
            names(u_dim[i]) <- i
            unshared = c(unshared, i)
            if (u_dim[[i]][2] > max_feats) {
                max_feats = u_dim[[i]][1]
            }
            ulist[[i]] = t(object@scale.unshared.data[[i]])
        }
    }
    ############## For every set of additional features less than the maximum, append an additional zero matrix s.t. it matches the maximum
    for (i in 1:length(object@scale.data)) {
        if (i %in% unshared) {
            mlist[[i]] <-  rbind(mlist[[i]], object@scale.unshared.data[[i]])
        }
        #For the U matrix with the maximum amount of features, append the whole thing
        else {
            mlist[[i]] <- rbind(mlist[[i]])
        }
    }

    X <- mlist
    ################# Create an 0 matrix the size of U for all U's, s.t. it can be stacked to W
    zero_matrix_u_full <- c()
    zero_matrix_u_partial <- c()
    for (i in 1:length(object@raw.data)) {
        if (i %in% unshared) {
            zero_matrix_u_full[[i]] <-
                matrix(0, nrow = u_dim[[i]][1], ncol = u_dim[[i]][2])
            zero_matrix_u_partial[[i]] <-
                matrix(0, nrow = u_dim[[i]][1], ncol = k)
        }
    }

    num_cells = c()
    for (i in 1:length(X)) {
        num_cells = c(num_cells, ncol(X[[i]]))
    }

    num_genes = length(object@var.genes)

    best_obj <- Inf
    for (i in 1:nrep) {
        print("Processing")
        current <- rand.seed + i - 1
        # initialization
        idX = list()
        for (i in 1:length(X)) {
            idX[[i]] = sample(1:num_cells[i], k)
        }
        V = list()

        #Establish V from only the RNA dimensions

        for (i in 1:length(X)) {
            V[[i]] = t(object@scale.data[[i]])[, idX[[i]]]
        }
        #Establish W from the shared gene dimensions

        W = matrix(abs(runif(num_genes * k, 0, 2)), num_genes, k)

        H = list()

        #Initialize U
        U = list()
        for (i in 1:length(X)) {
            if (i %in% unshared) {
                U[[i]] = t(ulist[[i]])[, idX[[i]]]
            }
        }

        iter = 0
        total_time = 0
        pb <- txtProgressBar(min = 0,
                             max = max.iters,
                             style = 3)
        sqrt_lambda = list()
        for (i in 1:length(X)) {
            sqrt_lambda[[i]] = sqrt(lambda[[i]])
        }
        ############################ Initial Training Objects

        obj_train_approximation = 0
        obj_train_penalty = 0

        for (i in 1:length(X)) {
            H[[i]] = matrix(abs(runif(k * num_cells[i], 0, 2)), k, num_cells[i])
            if (i %in% unshared) {
                obj_train_approximation = obj_train_approximation + norm(X[[i]] - (
                    rbind(W, zero_matrix_u_partial[[i]]) + rbind(V[[i]], U[[i]])
                ) %*% H[[i]], "F") ^ 2
                obj_train_penalty = obj_train_penalty + lambda[[i]] * norm(rbind(V[[i]], U[[i]]) %*% H[[i]], "F") ^
                    2
            }
            else {
                obj_train_approximation = obj_train_approximation + norm(X[[i]] - (W + V[[i]]) %*% H[[i]], "F") ^
                    2
                obj_train_penalty = obj_train_penalty + lambda[[i]] * norm(V[[i]] %*% H[[i]], "F") ^
                    2

            }
        }
        obj_train = obj_train_approximation + obj_train_penalty

        ######################### Initialize Object Complete ###########################
        ########################## Begin Updates########################################
        delta = Inf
        objective_value_list = list()

        iter = 1
        while (delta > thresh & iter <= max.iters) {
            iter_start_time = Sys.time()


            #H- Updates
            for (i in 1:length(X)) {
                if (!(i %in% unshared)) {
                    H[[i]] = solveNNLS(rbind((W + V[[i]]), sqrt_lambda[[i]] * V[[i]]),
                                       rbind(X[[i]], matrix(0, num_genes, xdim[[i]][2])))
                }
                else{
                    H[[i]] = solveNNLS(rbind(
                        rbind(W, zero_matrix_u_partial[[i]]) + rbind((V[[i]]), U[[i]]),
                        sqrt_lambda[[i]] * rbind(V[[i]], U[[i]])
                    ),
                    rbind((X[[i]]),
                          matrix(0, num_genes + u_dim[[i]][1], xdim[[i]][2])
                    ))
                }
            }

            #V - updates
            for (i in 1:length(X)) {
                V[[i]] = t(solveNNLS(
                    rbind(t(H[[i]]), sqrt_lambda[[i]] * t(H[[i]])),
                    rbind(
                        t(X[[i]][0:num_genes, ] - W %*% H[[i]]),
                        matrix(0, num_cells[i], num_genes)
                    )
                ))
            }
            ################################################# Updating U##################################

            for (i in 1:length(X)) {
                if (i %in% unshared) {
                    U[[i]] = t(solveNNLS(rbind(
                        t(H[[i]]), sqrt_lambda[[i]] * t(H[[i]])
                    ), rbind(
                        t(X[[i]][(num_genes + 1):(u_dim[[i]][1] + num_genes),]), t(zero_matrix_u_full[[i]])
                    )))
                }
            }


            ##############################################################################################
            ################################################# Updating W #################################
            H_t_stack = c()
            for (i in 1:length(X)) {
                H_t_stack = rbind(H_t_stack, t(H[[i]]))
            }
            diff_stack_w = c()
            for (i in 1:length(X)) {
                diff_stack_w = rbind(diff_stack_w, t(X[[i]][0:num_genes, ] - V[[i]] %*% H[[i]]))
            }
            W = t(solveNNLS(H_t_stack, diff_stack_w))

            ############################################################################################
            iter_end_time = Sys.time()
            iter_time = as.numeric(difftime(iter_end_time, iter_start_time, units = "secs"))
            total_time = total_time + iter_time

            #Updating training object
            obj_train_prev = obj_train
            obj_train_approximation = 0
            obj_train_penalty = 0



            for (i in 1:length(X)) {
                if (i %in% unshared) {
                    obj_train_approximation = obj_train_approximation + norm(X[[i]] - (
                        rbind(W, zero_matrix_u_partial[[i]]) + rbind(V[[i]], U[[i]])
                    ) %*% H[[i]], "F") ^ 2
                    obj_train_penalty = obj_train_penalty + lambda[[i]] *
                        norm(rbind(V[[i]], U[[i]]) %*% H[[i]], "F") ^ 2
                }
                else {
                    obj_train_approximation = obj_train_approximation + norm(X[[i]] - (W + V[[i]]) %*% H[[i]], "F") ^
                        2
                    obj_train_penalty = obj_train_penalty + lambda[[i]] *
                        norm(V[[i]] %*% H[[i]], "F") ^ 2

                }
            }

            obj_train = obj_train_approximation + obj_train_penalty
            delta = abs(obj_train_prev - obj_train) / mean(c(obj_train_prev, obj_train))
            iter = iter + 1
            setTxtProgressBar(pb = pb, value = iter)
        }
        setTxtProgressBar(pb = pb, value = max.iters)
        cat("\nCurrent seed ",
            current ,
            " current objective ",
            obj_train)
        if (obj_train < best_obj) {
            W_m <- W
            H_m <- H
            V_m <- V
            U_m <- U
            best_obj <- obj_train
            best_seed <- current
        }
    }

    rownames(W_m) = rownames(X[[1]][0:xdim[[i]][1], ])
    colnames(W_m) = NULL

    for (i in 1:length(X)) {
        if (i %in% unshared) {
            rownames(U_m[[i]]) = rownames(X[[i]][(num_genes + 1):(u_dim[[i]][1] + num_genes),])
            colnames(U_m[[i]]) = NULL
        }
        rownames(V_m[[i]]) = rownames(X[[i]][0:xdim[[i]][1], ])
        colnames(V_m[[i]]) = NULL
        colnames(H_m[[i]]) = colnames(X[[i]])
    }

    ################################## Returns Results Section #########################################################
    object@W <- t(W_m)
    for (i in 1:length(X)) {
        object@V[[i]] <- t(V_m[[i]])
        object@H[[i]] <- t(H_m[[i]])
        if (i %in% unshared) {
            object@U[[i]] <- t(U_m[[i]])
        }
    }
    titles <- names(object@raw.data)
    names(object@H) <- titles
    names(object@V) <- titles
    if (i %in% unshared) {
        names(object@U) <- titles
    }
    if (print.obj) {
        cat("\n", "Objective:", best_obj, "\n")
    }

    rel_cells = list()
    for (i in 1:length(X)) {
        rel_cells <- c(rel_cells, rownames(object@scale.data[[i]]))
    }
    rel_cells <- unlist(rel_cells)

    object@cell.data <- object@cell.data[rel_cells, ]
    cat("\n", "Best results with seed ", best_seed, ".\n", sep = "")
    return(object)
}


#' Calculate loadings for each factor
#'
#' Calculates the contribution of each factor of W,V, and U to the reconstruction.
#'
#' @param object \code{liger} object. Should call quantileNorm before calling.
#' @return A dataframe, such that each column represents the contribution of a specific matrix (W, V_1, V_2, etc. )
#' @export
calcNormLoadings = function(object) {
    H_norm = object@H.norm
    W_norm = object@W
    V_norm = object@V
    U_norm = object@U
    ##### Calculation of Contribution #########################
    w_loadings = list()
    u_loadings = list()
    for (i in 1:length(object@raw.data)) {
        u_loadings[[i]] = list()
    }
    v_loadings = list()
    for (i in 1:length(object@raw.data)) {
        v_loadings[[i]] = list()
    }
    for (i in 1:dim(object@H.norm)[[2]]) {
        hi = as.matrix(H_norm[, i])
        ####### Calculate W
        wi = t(as.matrix(W_norm[i, ]))
        hw = hi %*% wi
        forb_hw = norm(hw, type = "F") / dim(W_norm)[[2]]
        w_loadings = append(w_loadings, forb_hw)

        ###### Calculate V
        for (j in 1:length(object@raw.data)) {
            temp_v = t(as.matrix(V_norm[[j]][i, ]))
            hv_temp = hi %*% temp_v
            forb_hv = norm(hv_temp, type = "F") / dim(V_norm[[j]])[[2]]
            v_loadings[[j]] = append(v_loadings[[j]], forb_hv)
        }
        if (length(object@U) != 0) {
            ###### Calculate U
            for (j in 1:length(object@raw.data)) {
                if (length(object@U[[j]]) != 0) {
                    temp_u = t(as.matrix(U_norm[[j]][i, ]))
                    hu_temp = hi %*% temp_u
                    forb_hu = norm(hu_temp, type = "F") / dim(U_norm[[j]])[[2]]
                    u_loadings[[j]] = append(u_loadings[[j]], forb_hu)
                }
            }
        }
    }

    ################# Format the return object
    w_loadings = unlist(w_loadings)
    factors = 1:dim(object@H.norm)[[2]]
    results = data.frame(factors, w_loadings)

    # For all V
    for (j in 1:length(object@raw.data)) {
        results = cbind(results, unlist(v_loadings[[j]]))
        colnames(results)[[2 + j]] = paste0("V_", j, "_loadings")
    }
    if (length(object@U) != 0) {
        # For all U
        for (j in 1:length(object@raw.data)) {
            name_di = dim(results)[[2]]
            if (length(object@U[[j]]) != 0) {
                results = cbind(results, unlist(u_loadings[[j]]))
                colnames(results)[[name_di + 1]] = paste0("U_", j, "_loadings")
            }
        }
    }

    return(results)
}
