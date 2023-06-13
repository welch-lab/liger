#' Perform iNMF on scaled datasets, and include unshared features
#' @param object \linkS4class{liger} object. Should run
#' \code{\link{selectGenes}} with \code{unshared = TRUE} and then run
#' \code{\link{scaleNotCenter}} in advance.
#' @param k Integer, inner dimension of factorization (number of factors).
#' Default \code{30}.
#' @param lambda Numeric, the lambda penalty. Default \code{5}.
#' @param thresh Numeric, convergence threshold. Convergence occurs when
#' \eqn{|obj_0-obj|/(mean(obj_0,obj)) < thresh}. Default \code{1e-10}.
#' @param maxIter Maximum number of block coordinate descent iterations to
#' perform. Default \code{30}.
#' @param nrep Number of restarts to perform (iNMF objective function is
#' non-convex, so taking the best objective from multiple successive
#' initializations is recommended). For easier reproducibility, this increments
#' the random seed by 1 for each consecutive restart, so future factorizations
#' of the same dataset can be run with one rep if necessary. Default \code{1}.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @noRd
optimizeUANLS <- function(
        object,
        k = 30,
        lambda = 5,
        maxIter = 30,
        nrep = 1,
        thresh = 1e-10,
        seed = 1,
        verbose = getOption("ligerVerbose")
) {
    set.seed(seed)
    if (isTRUE(verbose))
        .log('Performing Factorization using UINMF and unshared features')
    if (length(lambda) == 1)
        lambda <- rep(lambda, length(object))
    else if (length(lambda) != length(object))
        stop("Vectorized lambda should have one value for each dataset.")

    # Get a list of all the matrices
    # scaleData in mlist: gene x cell
    mlist <- getMatrix(object, "scaleData", returnList = TRUE)
    mlist <- lapply(mlist, as.matrix)
    xdim <- lapply(mlist, dim)

    # Return what datasets have unshared features,
    # and the dimensions of those unshared features
    ulist <- getMatrix(object, "scaleUnsharedData", returnList = TRUE)
    ulist <- lapply(ulist, as.matrix)
    udim <- lapply(ulist, dim)
    unshared <- which(sapply(ulist, function(x) !is.null(x)))
    max_feats <- max(unlist(lapply(ulist, nrow)))

    # For every set of additional features less than the maximum,
    # append an additional zero matrix s.t. it matches the maximum
    for (i in seq_along(object)) {
        if (i %in% unshared)
            # X: (g+u) x c
            mlist[[i]] <- rbind(mlist[[i]], ulist[[i]])
    }

    X <- mlist

    # Create an 0 matrix the size of U for all U's, s.t. it can be stacked to W
    zero_matrix_u_full <- c()
    zero_matrix_u_partial <- c()
    for (i in seq_along(object)) {
        if (i %in% unshared) {
            # u x c
            zero_matrix_u_full[[i]] <-
                matrix(0, nrow = udim[[i]][1], ncol = udim[[i]][2])
            # u x k
            zero_matrix_u_partial[[i]] <-
                matrix(0, nrow = udim[[i]][1], ncol = k)
        }
    }

    nCells <- sapply(X, ncol)
    nGenes <- length(varFeatures(object))
    bestObj <- Inf

    for (i in seq(nrep)) {
        current <- seed + i - 1
        # initialization
        idX <- list()

        for (i in seq_along(X)) idX[[i]] <- sample(nCells[i], k)

        # Establish V from only the RNA dimensions
        # V matrices: g x k
        V <- list()
        for (i in seq_along(X))
            V[[i]] <- as.matrix(scaleData(object, i)[, idX[[i]]])
        # Establish W from the shared gene dimensions
        # W matrices: g x k
        W <- matrix(stats::runif(nGenes*k, 0, 2), nGenes, k)
        H <- list()
        # Initialize U
        U <- list()
        for (i in seq_along(X)) {
            if (i %in% unshared) {
                # u x k
                U[[i]] <- ulist[[i]][, idX[[i]]]
            }
        }

        iter <- 0
        total_time <- 0
        if (isTRUE(verbose))
            pb <- utils::txtProgressBar(min = 0, max = maxIter, style = 3)
        sqrtLambda <- lapply(lambda, sqrt)

        # Initial Training Objects

        obj_train_approximation <- 0
        obj_train_penalty <- 0

        for (i in seq_along(X)) {
            # H: k x c
            H[[i]] = matrix(stats::runif(k * nCells[i], 0, 2), k, nCells[i])
            if (i %in% unshared) {
                obj_train_approximation <- obj_train_approximation +
                    norm(X[[i]] - (rbind(W, zero_matrix_u_partial[[i]]) +
                                       rbind(V[[i]], U[[i]])) %*% H[[i]],
                         "F") ^ 2
                obj_train_penalty <- obj_train_penalty +
                    lambda[[i]] * norm(rbind(V[[i]], U[[i]]) %*% H[[i]],
                                       "F") ^ 2
            } else {
                obj_train_approximation <- obj_train_approximation +
                    norm(X[[i]] - (W + V[[i]]) %*% H[[i]], "F") ^ 2
                obj_train_penalty <- obj_train_penalty +
                    lambda[[i]] * norm(V[[i]] %*% H[[i]], "F") ^ 2

            }
        }
        obj_train <- obj_train_approximation + obj_train_penalty

        #### Initialize Object Complete
        #### Begin Updates
        delta <- Inf
        iter <- 1
        while (delta > thresh & iter <= maxIter) {
            iter_start_time <- Sys.time()

            # H - Updates
            # H: k x c
            for (i in seq_along(X)) {
                if (!(i %in% unshared)) {
                    H[[i]] <- solveNNLS(
                        rbind(W + V[[i]], sqrtLambda[[i]] * V[[i]]),
                        rbind(X[[i]], matrix(0, nGenes, xdim[[i]][2]))
                    )
                } else {
                    H[[i]] = solveNNLS(
                        rbind(
                            rbind(W, zero_matrix_u_partial[[i]]) +
                                rbind(V[[i]], U[[i]]),
                            sqrtLambda[[i]]*rbind(V[[i]], U[[i]])
                        ),
                        rbind(X[[i]],
                              matrix(0, nGenes + udim[[i]][1], xdim[[i]][2]))
                    )
                }
            }

            # V - updates
            for (i in seq_along(X)) {
                V[[i]] = t(solveNNLS(
                    rbind(t(H[[i]]), sqrtLambda[[i]] * t(H[[i]])),
                    rbind(
                        t(X[[i]][seq(nGenes), ] - W %*% H[[i]]),
                        matrix(0, nCells[i], nGenes)
                    )
                ))
            }

            # U - updates
            # U: u x k
            for (i in seq_along(X)) {
                if (i %in% unshared) {
                    U[[i]] = t(solveNNLS(
                        rbind(t(H[[i]]), sqrtLambda[[i]] * t(H[[i]])),
                        rbind(
                            t(X[[i]][seq(nGenes + 1, udim[[i]][1] + nGenes),]),
                            t(zero_matrix_u_full[[i]])
                        )
                    ))
                }
            }


            # W - updates
            # H_t_stack: C x k
            # ("C" for all cells, "c" for number of cells in each dataset)
            H_t_stack <- t(Reduce(cbind, H))
            # diff_stack_w: C x g
            diff_stack_w = c()
            for (i in seq_along(X)) {
                diff_stack_w = cbind(
                    diff_stack_w,
                    X[[i]][seq(nGenes), ] - V[[i]] %*% H[[i]]
                )
            }
            diff_stack_w <- t(diff_stack_w)
            # W: g x k
            W = t(solveNNLS(H_t_stack, diff_stack_w))

            ## End of updates
            iter_end_time <- Sys.time()
            iter_time <- as.numeric(difftime(iter_end_time, iter_start_time,
                                             units = "secs"))
            total_time <- total_time + iter_time

            #Updating training object
            obj_train_prev <- obj_train
            obj_train_approximation <- 0
            obj_train_penalty <- 0

            for (i in seq_along(X)) {
                if (i %in% unshared) {
                    obj_train_approximation <- obj_train_approximation +
                        norm(
                            X[[i]] - (rbind(W, zero_matrix_u_partial[[i]]) +
                                          rbind(V[[i]], U[[i]])) %*% H[[i]],
                            "F") ^ 2
                    obj_train_penalty <- obj_train_penalty + lambda[[i]] *
                        norm(rbind(V[[i]], U[[i]]) %*% H[[i]], "F") ^ 2
                }
                else {
                    obj_train_approximation <- obj_train_approximation +
                        norm(X[[i]] - (W + V[[i]]) %*% H[[i]], "F") ^ 2
                    obj_train_penalty = obj_train_penalty + lambda[[i]] *
                        norm(V[[i]] %*% H[[i]], "F") ^ 2
                }
            }

            obj_train <- obj_train_approximation + obj_train_penalty
            delta <- abs(obj_train_prev - obj_train) /
                mean(c(obj_train_prev, obj_train))
            iter <- iter + 1
            if (isTRUE(verbose)) utils::setTxtProgressBar(pb = pb, value = iter)
        }
        if (isTRUE(verbose)) {
            utils::setTxtProgressBar(pb = pb, value = maxIter)
            cat("\n")
            .log("Current seed: ", current, "; Current objective: ", obj_train)
        }
        if (obj_train < bestObj) {
            W_m <- W
            H_m <- H
            V_m <- V
            U_m <- U
            bestObj <- obj_train
            best_seed <- current
        }
    }

    factorNames <- paste0("factor_", seq(k))

    rownames(W_m) <- varFeatures(object)
    colnames(W_m) <- factorNames
    object@W <- W_m

    for (i in seq_along(object)) {
        ld <- dataset(object, i)

        rownames(V_m[[i]]) <- varFeatures(object)
        colnames(V_m[[i]]) <- factorNames
        ld@V <- V_m[[i]]

        rownames(H_m[[i]]) <- factorNames
        colnames(H_m[[i]]) <- colnames(X[[i]])
        ld@H <- H_m[[i]]

        if (i %in% unshared) {
            rownames(U_m[[i]]) <- ld@varUnsharedFeatures
            colnames(U_m[[i]]) <- factorNames
            ld@U <- U_m[[i]]
        }
        datasets(object, check = FALSE)[[i]] <- ld
    }
    param <- list(k = k, lambda = lambda)
    object@uns$factorization <- param

    if (isTRUE(verbose))
        .log("Objective: ", bestObj, "\nBest results with seed: ", best_seed)

    return(object)
}
