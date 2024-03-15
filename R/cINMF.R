#' Perform consensus iNMF on scaled datasets
#' @description
#' Performs consensus integrative non-negative matrix factorization (c-iNMF)
#' to return factorized \eqn{H}, \eqn{W}, and \eqn{V} matrices. We run the
#' regular iNMF multiple times with different random starts, and then take the
#' consensus of frequently appearing factors from gene loading matrices, \eqn{W}
#' and \eqn{V}. The cell factor loading \eqn{H} matrices are eventually solved
#' with the consensus \eqn{W} and \eqn{V} matrices.
#'
#' Please see \code{\link{runINMF}} for detailed introduction to the regular
#' iNMF algorithm which is run multiple times in this function.
#'
#' The consensus iNMF algorithm is developed basing on the consensus NMF (cNMF)
#' method (D. Kotliar et al., 2019).
#' @param object A \linkS4class{liger} object or a Seurat object with
#' non-negative scaled data of variable features (Done with
#' \code{\link{scaleNotCenter}}).
#' @param k Inner dimension of factorization (number of factors). Generally, a
#' higher \code{k} will be needed for datasets with more sub-structure. Default
#' \code{20}.
#' @param lambda Regularization parameter. Larger values penalize
#' dataset-specific effects more strongly (i.e. alignment should increase as
#' \code{lambda} increases). Default \code{5}.
#' @param rho Numeric number between 0 and 1. Fraction for determining the
#' number of nearest neighbors to look at for consensus (by
#' \code{rho * nRandomStarts}). Default \code{0.3}.
#' @param nIteration Total number of block coordinate descent iterations to
#' perform. Default \code{30}.
#' @param nRandomStarts Number of replicate runs for creating the pool of
#' factorization results. Default \code{10}.
#' @param HInit Initial values to use for \eqn{H} matrices. A list object where
#' each element is the initial \eqn{H} matrix of each dataset. Default
#' \code{NULL}.
#' @param WInit Initial values to use for \eqn{W} matrix. A matrix object.
#' Default \code{NULL}.
#' @param VInit Initial values to use for \eqn{V} matrices. A list object where
#' each element is the initial \eqn{V} matrix of each dataset. Default
#' \code{NULL}.
#' @param seed Random seed to allow reproducible results. Default \code{1}.
#' @param nCores The number of parallel tasks to speed up the computation.
#' Default \code{2L}. Only supported for platform with OpenMP support.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param ... Arguments passed to methods.
#' @rdname runCINMF
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
#' @references
#' Joshua D. Welch and et al., Single-Cell Multi-omic Integration Compares and
#' Contrasts Features of Brain Cell Identity, Cell, 2019
#'
#' Dylan Kotliar and et al., Identifying gene expression programs of cell-type
#' identity and cellular activity with single-cell RNA-Seq, eLife, 2019
#' @examples
#' \donttest{
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' if (requireNamespace("RcppPlanc", quietly = TRUE)) {
#'     pbmc <- runCINMF(pbmc)
#' }
#' }
runCINMF <- function(
        object,
        k = 20,
        lambda = 5.0,
        rho = 0.3,
        ...
) {
    UseMethod("runCINMF", object)
}

#' @rdname runCINMF
#' @export
#' @method runCINMF liger
runCINMF.liger <- function(
        object,
        k = 20,
        lambda = 5.0,
        rho = 0.3,
        nIteration = 30,
        nRandomStarts = 10,
        HInit = NULL,
        WInit = NULL,
        VInit = NULL,
        seed = 1,
        nCores = 2L,
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    .checkObjVersion(object)
    object <- recordCommand(object, ..., dependencies = "RcppPlanc")
    object <- removeMissing(object, orient = "cell", verbose = verbose)
    data <- lapply(datasets(object), function(ld) {
        if (is.null(scaleData(ld)))
            cli::cli_abort("Scaled data not available. Run {.fn scaleNotCenter} first.")
        return(scaleData(ld))
    })
    dataClasses <- sapply(data, function(x) class(x)[1])
    if (!all(dataClasses == dataClasses[1])) {
        cli::cli_abort("Currently the {.field scaledData} of all datasets have to be of the same class.")
    }
    out <- .runCINMF.list(
        object = data,
        k = k,
        lambda = lambda,
        rho = rho,
        nIteration = nIteration,
        nRandomStarts = nRandomStarts,
        HInit = HInit,
        WInit = WInit,
        VInit = VInit,
        seed = seed,
        nCores = nCores,
        verbose = verbose
    )
    # return(out)
    object@W <- out$W
    for (d in names(object)) {
        object@datasets[[d]]@H <- out$H[[d]]
        object@datasets[[d]]@V <- out$V[[d]]
        # ld <- dataset(object, d)
        # ld@H <- out$H[[d]]
        # ld@V <- out$V[[d]]
        # datasets(object, check = FALSE)[[d]] <- ld
    }
    object@uns$factorization <- list(k = k, lambda = lambda)
    return(object)
}

#' @rdname runCINMF
#' @export
#' @param datasetVar Metadata variable name that stores the dataset source
#' annotation. Default \code{"orig.ident"}.
#' @param layer For Seurat>=4.9.9, the name of layer to retrieve input
#' non-negative scaled data. Default \code{"ligerScaleData"}. For older Seurat,
#' always retrieve from \code{scale.data} slot.
#' @param assay Name of assay to use. Default \code{NULL} uses current active
#' assay.
#' @param reduction Name of the reduction to store result. Also used as the
#' feature key. Default \code{"cinmf"}.
#' @method runCINMF Seurat
runCINMF.Seurat <- function(
        object,
        k = 20,
        lambda = 5.0,
        rho = 0.3,
        datasetVar = "orig.ident",
        layer = "ligerScaleData",
        assay = NULL,
        reduction = "cinmf",
        nIteration = 30,
        nRandomStarts = 10,
        HInit = NULL,
        WInit = NULL,
        VInit = NULL,
        seed = 1,
        nCores = 2L,
        verbose = getOption("ligerVerbose", TRUE),
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
            cli::cli_abort(
                c("x" = "Negative data encountered for integrative {.emph Non-negative} Matrix Factorization.",
                  "!" = "Please run {.fn scaleNotCenter} first.")
            )
        }
    }

    res <- .runCINMF.list(
        object = Es,
        k = k,
        lambda = lambda,
        rho = rho,
        nIteration = nIteration,
        nRandomStarts = nRandomStarts,
        HInit = HInit,
        WInit = WInit,
        VInit = VInit,
        seed = seed,
        nCores = nCores,
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
.runCINMF.list <- function(
        object,
        k = 20,
        lambda = 5.0,
        rho = 0.3,
        nIteration = 30,
        nRandomStarts = 10,
        HInit = NULL,
        WInit = NULL,
        VInit = NULL,
        seed = 1,
        nCores = 2L,
        verbose = getOption("ligerVerbose", TRUE)
) {
    if (!requireNamespace("RcppPlanc", quietly = TRUE)) # nocov start
        cli::cli_abort(
            "Package {.pkg RcppPlanc} is required for c-iNMF integration.
        Please install it by command:
        {.code devtools::install_github('welch-lab/RcppPlanc')}") # nocov end
    if (nRandomStarts <= 1) {
        cli::cli_abort("{.var nRandomStarts} must be greater than 1 for taking the consensus.")
    }
    if (rho <= 0 || rho > 1) {
        cli::cli_abort(c("Invalid value for {.var rho}.",
                         "i" = "{.var rho} must be in the range (0, 1]."))
    }
    if (round(rho * nRandomStarts) < 1) {
        cli::cli_abort(c("Too few nearest neighbors to take the consensus. ",
                         "i" = "Please use a larger {.var rho} or/and a larger {.var nRandomStarts}."))
    }

    barcodeList <- lapply(object, colnames)
    allFeatures <- lapply(object, rownames)
    features <- Reduce(.same, allFeatures)

    Ws <- list()
    # Each dataset got a list of replicate V matrices
    Vs <- stats::setNames(rep(list(list()), length(object)), names(object))
    # Each dataset got a H matrices
    Hs <- stats::setNames(rep(list(NULL), length(object)), names(object))
    if (isTRUE(verbose))
        cli::cli_progress_bar(name = "Replicating iNMF runs",
                              status = sprintf("[ 0 / %d ]", nRandomStarts),
                              total = nRandomStarts)
    set.seed(seed)
    for (i in seq(nRandomStarts)) {
        repName <- paste0("rep_", i)
        out <- RcppPlanc::inmf(objectList = object, k = k, lambda = lambda,
                               niter = nIteration, Hinit = HInit,
                               Vinit = VInit, Winit = WInit, nCores = nCores,
                               verbose = FALSE)
        Ws[[repName]] <- out$W
        for (n in names(object)) Vs[[n]][[repName]] <- out$V[[n]]
        for (n in names(object)) Hs[[n]][[repName]] <- out$H[[n]]
        if (isTRUE(verbose))
            cli::cli_progress_update(
                status = sprintf("[ %d / %d ]", i, nRandomStarts), set = i
            )
    }
    # return(list(W = Ws, V = Vs, H = Hs))
    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Taking the consensus")
    # Attempt 1 ####
    # W <- takeConsensus(Ws, rho = rho, tao = tao)
    # Vs <- lapply(Vs, takeConsensus, rho = rho, tao = tao)

    # Attempt 2 ####
    # factorSelect <- takeConsensus2(Ws, rho = rho, tao = tao)
    # W <- Reduce(cbind, Ws)[, factorSelect$idx]
    # W <- colNormalize_dense_cpp(W, L = 2)
    # W <- colAggregateMedian_dense_cpp(W, group = factorSelect$cluster - 1, n = k)
    # W <- colNormalize_dense_cpp(W, L = 1)
    # for (i in seq_along(object)) {
    #     V <- Reduce(cbind, Vs[[i]])[, factorSelect$idx]
    #     V <- colNormalize_dense_cpp(V, L = 2)
    #     V <- colAggregateMedian_dense_cpp(V, group = factorSelect$cluster - 1, n = k)
    #     Vs[[i]] <- colNormalize_dense_cpp(V, L = 1)
    # }

    # Attempt 3 ####
    # W <- Reduce(cbind, Ws)
    # selection <- cluster_sel(W, nNeighbor = 3, k = k, resolution = 0.8)
    # W <- W[, selection$idx]
    # W <- colNormalize_dense_cpp(W, L = 2)
    # print(table(selection$cluster - 1))
    # W <- colAggregateMedian_dense_cpp(W, group = selection$cluster, n = k)
    # W <- colNormalize_dense_cpp(W, L = 1)
    # for (i in seq_along(Vs)) {
    #     V <- Reduce(cbind, Vs[[i]])
    #     V <- V[, selection$idx]
    #     V <- colNormalize_dense_cpp(V, L = 2)
    #     V <- colAggregateMedian_dense_cpp(V, group = selection$cluster, n = k)
    #     Vs[[i]] <- colNormalize_dense_cpp(V, L = 1)
    # }

    # Attempt 4 ####
    # nNeighbor <- round(rho * nRandomStarts)
    # geneLoadings <- list(Reduce(cbind, Ws))
    # for (i in seq_along(object)) {
    #     geneLoadings[[i + 1]] <- Reduce(cbind, Vs[[i]])
    # }
    # geneLoadings <- lapply(geneLoadings, colNormalize_dense_cpp, L = 2)
    # # geneLoadings is a list, each element is a ngene by (nFactor * nRandomStart) matrix
    # # The first is for W, the rest are for Vs
    # selection <- factor_cluster_sel(geneLoadings, nNeighbor = nNeighbor,
    #                                 minWeight = 0.6, k = k, resolution = 0.2)
    # W <- geneLoadings[[1]][, selection$idx]
    # W <- colAggregateMedian_dense_cpp(W, group = selection$cluster, n = k)
    # W <- colNormalize_dense_cpp(W, L = 1)
    # for (i in seq_along(Vs)) {
    #     V <- geneLoadings[[i + 1]][, selection$idx]
    #     V <- colAggregateMedian_dense_cpp(V, group = selection$cluster, n = k)
    #     Vs[[i]] <- colNormalize_dense_cpp(V, L = 1)
    # }

    # Attempt 5
    nNeighbor <- round(rho * nRandomStarts)
    geneLoadings <- list(Reduce(cbind, Ws))
    for (i in seq_along(object)) {
        geneLoadings[[i + 1]] <- Reduce(cbind, Vs[[i]])
    }
    geneLoadings <- lapply(geneLoadings, colNormalize_dense_cpp, L = 2)
    # geneLoadings is a list, each element is a ngene by (nFactor * nRandomStart) matrix
    # The first is for W, the rest are for Vs
    selection <- factor_cluster_sel(geneLoadings, nNeighbor = nNeighbor,
                                    minWeight = 0.6, k = k, resolution = 0.2)
    W <- geneLoadings[[1]][, selection$idx]
    W <- colAggregateMedian_dense_cpp(W, group = selection$cluster, n = k)
    W <- colNormalize_dense_cpp(W, L = 1)
    for (i in seq_along(Vs)) {
        V <- geneLoadings[[i + 1]][, selection$idx]
        V <- colAggregateMedian_dense_cpp(V, group = selection$cluster, n = k)
        Vs[[i]] <- colNormalize_dense_cpp(V, L = 1)
    }
    # Vs <- lapply(seq_along(object), function(i) {
    #     matrix(stats::runif(nrow(W) * k, 0, 2), nrow(W), k)
    # })
    Hs <- lapply(seq_along(object), function(i) {
        matrix(stats::runif(ncol(object[[i]]) * k, 0, 2), ncol(object[[i]]), k)
    })

    if (isTRUE(verbose)) cli::cli_process_done(id = cliID)

    msg <- "ANLS optimization with consensus fixed"
    if (isTRUE(verbose)) cliID <- cli::cli_process_start(msg = msg)
    for (iter in seq_len(nIteration*2)) {
        for (i in seq_along(object)) {
            Hs[[i]] <- inmf_solveH(NULL, W, Vs[[i]], object[[i]], lambda, nCores = nCores)
        }
        for (i in seq_along(object)) {
            Vs[[i]] <- inmf_solveV(Hs[[i]], W, NULL, object[[i]], lambda, nCores = nCores)
        }
    }

    # for (i in seq_along(object)) {
    #     Vs[[i]] <- inmf_solveV(Hs[[i]], W, NULL, object[[i]], lambda, nCores = nCores)
    # }
    objErr <- sum(sapply(seq_along(object), function(i)
        inmf_objErr_i(H = Hs[[i]], W = W, V = Vs[[i]], E = object[[i]],
                      lambda = lambda)
    ))
    if (isTRUE(verbose))
        cli::cli_process_done(id = cliID, msg_done = paste(msg, " ... objective error: {objErr}"))

    factors <- paste0("Factor_", seq(k))
    dimnames(W) <- list(features, factors)
    for (i in seq_along(object)) {
        dimnames(Hs[[i]]) <- list(barcodeList[[i]], factors)
        Hs[[i]] <- t(Hs[[i]])
        dimnames(Vs[[i]]) <- list(features, factors)
    }
    names(Hs) <- names(Vs) <- names(object)
    result <- list(W = W, H = Hs, V = Vs, objErr = objErr)
    return(result)

    # out <- .runINMF.list(object, k = k, lambda = lambda, nIteration = 1,
    #                      HInit = Hs, WInit = W, VInit = Vs, verbose = FALSE)
    # if (isTRUE(verbose))
    #     cli::cli_process_done(id = cliID, msg_done = paste(msg, " ... objective error: {out$objErr}"))
    # return(out)
}

# takeConsensus2 <- function(matList, rho = 0.3, tao = 0.1) {
#     # 2nd attempt of methods
#     # Return the selection of factors from the replicate pool as well as the
#     # grouping on the selected ones. So this part of information will be used
#     # to select the V factors that represent the same thing as what we want
#     # from W
#
#     ## According to cNMF method.
#     ## G - The program matrix. what we call W or V
#     ## rho - fraction parameter to set number of nearest neighbors to look at
#     ## tao - distance threshold to filter out factors
#
#     ## matList are matrices of dimensionality gene x factor
#     # Step 1: Concatenate and l2 normalize
#     ## R - a number of replicates
#     G_R <- Reduce(cbind, matList) # ngene by (nFactor * nRandomStart)
#     G_R <- colNormalize_dense_cpp(G_R, L = 2)
#
#     # Step 2: Find kNN matching for each component (factor)
#     # and filter by mean distance against the kNN of each
#     nNeighbor <- round(rho * length(matList))
#     knn <- RANN::nn2(t(G_R), k = nNeighbor + 1)
#     nnDistMeans <- rowMeans(knn$nn.dist[,2:(nNeighbor + 1)])
#     selectedIdx <- nnDistMeans < tao
#     # `select_factor_cpp` is a C++ function that returns the indices of the
#     # factors that are selected by examining whether the mean euclidean distance
#     # to its NNs is less than the threshold `tao`.
#     # selectedIdx <- select_factor_cpp(all_data = G_R,
#     #                                  knn = knn$nn.idx - 1,
#     #                                  threshold = tao)
#
#     G_R <- G_R[, selectedIdx] # genes by selected factors
#     if (ncol(G_R) < ncol(matList[[1]])) {
#         cli::cli_abort(c("Too few factors are selected from the pool to take the consensus.",
#                          "i" = "Please try: ",
#                          "*" = "a larger {.var tao} for loosen threshold",
#                          "*" = "a larger {.var nRandomStarts} for more replicates to be used."))
#     }
#     # Step 3: Kmeans clustering to get the k consensus
#     cluster <- stats::kmeans(t(G_R), centers = ncol(matList[[1]]))$cluster
#     return(list(idx = selectedIdx, cluster = cluster))
#     # G_consensus <- colAggregateMedian_dense_cpp(x = G_R, group = cluster - 1, n = ncol(matList[[1]]))
#     # G_consensus <- colNormalize_dense_cpp(G_consensus, L = 1)
#     # return(G_consensus)
# }

inmf_solveH <- function(H, W, V, E, lambda, nCores = 2L) {
    WV <- W + V
    CtC <- t(WV) %*% WV + lambda * t(V) %*% V
    CtB <- t(WV) %*% E
    H <- RcppPlanc::bppnnls_prod(CtC, as.matrix(CtB), nCores = nCores)
    return(t(H)) # return cell x factor
}

inmf_solveV <- function(H, W, V, E, lambda, nCores = 2L) {
    HtH <- t(H) %*% H
    CtC <- (1 + lambda) * HtH
    CtB <- t(H) %*% t(E) - HtH %*% t(W)
    V <- RcppPlanc::bppnnls_prod(CtC, as.matrix(CtB), nCores = nCores)
    return(t(V))
}

# inmf_solveW <- function(Hs, W, Vs, Es, lambda, nCores = 2L) {
#     CtC <- matrix(0, ncol(Vs[[1]]), ncol(Vs[[1]]))
#     CtB <- matrix(0, ncol(Vs[[1]]), nrow(Vs[[1]]))
#     for (i in seq_along(Es)) {
#         HtH <- t(Hs[[i]]) %*% Hs[[i]]
#         CtC <- CtC + HtH
#         CtB <- CtB + as.matrix(t(Hs[[i]]) %*% t(Es[[i]])) - HtH %*% t(Vs[[i]])
#     }
#     W <- RcppPlanc::bppnnls_prod(CtC, CtB, nCores = nCores)
#     return(t(W))
# }

inmf_objErr_i <- function(H, W, V, E, lambda) {
    # Objective error function was originally stated as:
    # obj_i = ||E_i - (W + V_i)*H_i||_F^2 + lambda * ||V_i*H_i||_F^2
    # (Use caution with the matrix dimensionality, as we might transpose them
    # occasionally in different implementation.)
    #
    # Let L = W + V
    # ||E - LH||_F^2 = ||E||_F^2 - 2*Tr(Ht*(Et*L)) + Tr((Lt*L)*(Ht*H))
    # ||V*H||_F^2 = Tr((Vt*V)*(Ht*H))
    # This improves the performance in both speed and memory usage.
    L <- W + V
    sqnormE <- Matrix::norm(E, "F")^2
    LtL <- t(L) %*% L
    HtH <- t(H) %*% H
    TrLtLHtH <- sum(diag(LtL %*% HtH))
    EtL <- t(E) %*% L
    TrHtEtL <- sum(Matrix::diag(t(H) %*% EtL))
    VtV <- t(V) %*% V
    TrVtVHtH <- sum(diag(VtV %*% HtH))
    obj <- sqnormE - 2 * TrHtEtL + TrLtLHtH + lambda * TrVtVHtH
    return(obj)
}

factor_cluster_sel <- function(geneLoadings, nNeighbor = 3, minWeight = 0.6,
                               k = 20, resolution = 0.2) {
    graphList <- lapply(geneLoadings, function(w) {
        knn <- RANN::nn2(t(w), k = nNeighbor + 1)
        target <- as.integer(t(knn$nn.idx))
        source <- rep(seq_len(nrow(knn$nn.idx)), each = ncol(knn$nn.idx))
        weight <- as.numeric(t(knn$nn.dists))
        weight <- exp(-weight/0.5)
        graphmat <- cbind(source, target, weight)
        graphmat <- graphmat[graphmat[,3] > minWeight,]
        return(graphmat)
    })
    graphmat <- Reduce(rbind, graphList)

    # Leiden cluster. CPP implementation expects 0-based vertex index
    edgevec <- as.integer(t(graphmat[,1:2])) - 1L
    cluster <- leidenAlg::find_partition_with_rep_rcpp(
        edgelist = edgevec, edgelist_length = length(edgevec),
        num_vertices = ncol(geneLoadings[[1]]), direction = FALSE, niter = 5,
        edge_weights = graphmat[,3], resolution = resolution, nrep = 20
    )

    # Select top k clusters
    selIdx <- cluster %in% names(sort(table(cluster), decreasing = TRUE))[seq_len(k)]
    cluster <- cluster[selIdx]
    cluster <- as.integer(factor(cluster)) - 1L
    return(list(idx = selIdx, cluster = cluster))
}
