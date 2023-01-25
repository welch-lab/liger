#' Perform factorization for new value of k
#' @description This uses an efficient strategy for updating that takes
#' advantage of the information in the existing factorization. It is most
#' recommended for values of \code{k.new} smaller than current value (\code{k},
#' which is set when running \code{\link{optimizeALS}}), where it
#' is more likely to speed up the factorization.
#' @param object \linkS4class{liger} object. Should call
#' \code{\link{optimizeALS}} in advance.
#' @param k.new Number of factors of factorization.
#' @param lambda Regularization parameter. By default, this will use the lambda
#' last used with \code{\link{optimizeALS}}.
#' @param thresh Convergence threshold. Convergence occurs when
#' \eqn{|obj_0-obj|/(mean(obj0,obj)) < thresh}. Default \code{1e-4}.
#' @param max.iters Maximum number of block coordinate descent iterations to
#' perform. Default \code{100}.
#' @param rand.seed Random seed to set. Only relevant if \code{k.new} is greater
#' than \code{k}. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{TRUE}.
#' @return \code{object} with \code{W} slot updated with the new \eqn{W}
#' matrix, and the \code{H} and \code{V} slots of each
#' \linkS4class{ligerDataset} object in the \code{datasets} slot updated with
#' the new dataset specific \eqn{H} and \eqn{V} matrix, respectively.
#' @export
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
    if (k.new == k) {
        return(object)
    }
    # g x c
    E <- lapply(datasets(object), function(ld) scale.data(ld))
    # g x k
    W <- object@W
    # g x k
    V <- lapply(datasets(object), function(ld) ld@V)
    # k x c
    H <- lapply(datasets(object), function(ld) ld@H)
    nGenes <- nrow(W)
    nDatasets <- length(object)
    nCells <- sapply(datasets(object), ncol)
    if (isTRUE(verbose)) .log("Updating initial input matrices")
    if (k.new > k) {
        set.seed(rand.seed)
        sqrtLambda <- sqrt(lambda)
        # Initialize W_new g x k_diff
        # TODO: Directly initialize with nGene x k.new-k instead of transposing
        # Doing it now because need to reproduce old result
        W_new <- t(matrix(runif(nGenes*(k.new - k), 0, 2), k.new - k, nGenes))
        # Initialize V_new g x k_diff
        V_new <- lapply(seq(nDatasets), function(i)
            t(matrix(runif(nGenes*(k.new - k), 0, 2), k.new - k, nGenes)))
        # H_new k_diff x c
        H_new <- lapply(seq(nDatasets), function(i) {
            solveNNLS(rbind(W_new + V_new[[i]], sqrtLambda * V_new[[i]]),
                      rbind(E[[i]] - (W + V[[i]]) %*% H[[i]],
                            matrix(0, nGenes, nCells[i])))
        })
        V_new <- lapply(seq(nDatasets), function(i) {
            t(solveNNLS(
                t(cbind(H_new[[i]], sqrtLambda * H_new[[i]])),
                t(cbind(
                    E[[i]] - (W + V[[i]]) %*% H[[i]] - W_new %*% H_new[[i]],
                    matrix(0, nGenes, nCells[[i]])))
            ))
        })
        W_new <- t(solveNNLS(
            t(Reduce(cbind, H_new)),
            t(Reduce(cbind,  lapply(seq(nDatasets), function(i)
                E[[i]] - (W + V[[i]]) %*% H[[i]] - V_new[[i]] %*% H_new[[i]])))
        ))
        # H k.new x c
        H <- lapply(seq(nDatasets), function(i) rbind(H[[i]], H_new[[i]]))
        # V&W g x k.new
        V <- lapply(seq(nDatasets), function(i) cbind(V[[i]], V_new[[i]]))
        W <- cbind(W, W_new)
    } else {
        deltas <- rep(0, k)
        for (i in seq(nDatasets))
            deltas <- deltas + sapply(seq(k), function(ki)
                norm((W[, ki] + V[[i]][, ki]) %*% t(H[[i]][ki,]), "F")
            )
        k.use <- order(deltas, decreasing = TRUE)[seq(k.new)]
        W <- W[, k.use]
        V <- lapply(V, function(x) x[, k.use])
        H <- lapply(H, function(x) x[k.use, ])
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
