optimizeNewK <- function(
        object,
        k.new,
        lambda = 5,
        thresh = 1e-4,
        max.iters = 100,
        rand.seed = 1,
        verbose = TRUE
) {
    .checkValidFactorResult(object)
    if (is.null(lambda)) {
        #lambda <- object@parameters$lambda
    }
    k <- ncol(object@W)
    if (is.null(k)) {
        stop("Unable to detect valid existing factorization result. ",
             "Please run `optimizeALS()` or `online_iNMF()` first.")
    }
    if (k.new == k) {
        return(object)
    }

    H <- lapply(datasets(object), function(ld) t(ld@H))
    W <- t(object@W)
    V <- lapply(datasets(object), function(ld) t(ld@V))
    E <- lapply(datasets(object), function(ld) t(scale.data(ld)))
    nGenes <- ncol(W)
    nDatasets <- length(object)
    nCells <- sapply(datasets(object), ncol)
    if (isTRUE(verbose)) .log("Updating initial input matrices")
    if (k.new > k) {
        set.seed(rand.seed)
        sqrtLambda <- sqrt(lambda)
        W_new <- matrix(runif(nGenes*(k.new - k), 0, 2), k.new - k, nGenes)
        V_new <- lapply(seq(nDatasets), function(i)
            matrix(runif(nGenes*(k.new - k), 0, 2), k.new - k, nGenes))
        H_new <- lapply(nCells, function(n)
            matrix(runif(n*(k.new - k), 0, 2), n, k.new - k))
        H_new <- lapply(seq(nDatasets), function(i) {
            t(solveNNLS(rbind(t(W_new + V_new[[i]]),
                              sqrtLambda * t(V_new[[i]])),
                        rbind(t(E[[i]] - H[[i]] %*% (W + V[[i]])),
                              matrix(0, nGenes, nCells[i]))
                        ))
        })
        V_new <- lapply(seq(nDatasets), function(i) {
            solveNNLS(
                rbind(H_new[[i]], sqrtLambda * H_new[[i]]),
                rbind(
                    E[[i]] - H[[i]] %*% (W + V[[i]]) - H_new[[i]] %*% W_new,
                    matrix(0, nCells[[i]], nGenes)
                )
            )
        })
        W_new <- solveNNLS(rbind.fill.matrix(H_new),
                           rbind.fill.matrix(
                               lapply(seq(nDatasets), function(i) {
                                   E[[i]] - H[[i]] %*% (W + V[[i]]) -
                                       H_new[[i]] %*% V_new[[i]]
                               })))
        H <- lapply(seq(nDatasets), function(i) t(cbind(H[[i]], H_new[[i]])))
        V <- lapply(seq(nDatasets), function(i) t(rbind(V[[i]], V_new[[i]])))
        W <- t(rbind(W, W_new))
    } else {
        deltas <- rep(0, k)
        for (i in seq(nDatasets))
            deltas <- deltas + sapply(seq(k), function(x)
                norm(H[[i]][, k] %*% t(W[k, ] + V[[i]][k, ]), "F")
            )
        k.use <- order(deltas, decreasing = TRUE)[seq(k.new)]
        W <- t(W[k.use, ])
        H <- lapply(H, function(x) t(x[, k.use]))
        V <- lapply(V, function(x) t(x[k.use, ]))
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
