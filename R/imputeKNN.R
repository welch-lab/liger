#' Impute the query cell expression matrix
#' @description Impute query features from a reference dataset using KNN.
#' @param object \linkS4class{liger} object with aligned factor loading computed
#' in advance.
#' @param nNeighbors The maximum number of nearest neighbors to search. Default
#' \code{20}.
#' @param reference Name of a dataset containing values to impute into query
#' dataset(s).
#' @param queries Names of datasets to be augmented by imputation. Should not
#' include \code{reference}. Default \code{NULL} uses all datasets except the
#' reference.
#' @param weight Logical. Whether to use KNN distances as weight matrix. Default
#' \code{FALSE}.
#' @param norm Logical. Whether to normalize the imputed data. Default
#' \code{TRUE}.
#' @param scale Logical. Whether to scale but not center the imputed data.
#' Default \code{TRUE}.
#' @param verbose Logical. Whether to show the progress and information. Default
#' \code{TRUE}.
#' @param ... Optional arguments to be passed to \code{\link{normalize}} when
#' \code{norm = TRUE}.
#' @param knn_k \bold{Deprecated}. See Usage section for replacement.
#' @return The input \code{object} where queried \linkS4class{ligerDataset}
#' objects in \code{datasets} slot for the queries are replaced. These
#' \linkS4class{ligerDataset} objects contain only the imputed data in the
#' \code{raw.data} slot. They contain newly computed \code{norm.data} and
#' \code{scale.data} optionally depending on arguments \code{norm} and
#' \code{scale}, respectively.
#' @export
imputeKNN <- function(
        object,
        reference,
        queries = NULL,
        nNeighbors = 20,
        weight = TRUE,
        norm = TRUE,
        scale = TRUE,
        verbose = TRUE,
        ...,
        # Deprecated coding style,
        knn_k = nNeighbors
) {
    .deprecateArgs(list(knn_k = "nNeighbors"),
                   call = rlang::call_args(match.call()))
    if (is.null(getMatrix(object, "H.norm")))
        stop("Aligned factor loading has to be available for imputation. ",
             "Please run `quantile_norm()` in advance.")
    if (isTRUE(verbose)) {
        warning("This function will discard the raw data previously stored ",
                "in the liger object and replace the `raw.data` slot with the ",
                "imputed data.", immediate. = TRUE)
    }
    if (length(reference) > 1) {
        stop("Can only have ONE reference dataset")
    }
    reference <- .checkUseDatasets(object, reference)
    showQueryWarn <- ifelse(is.null(queries), FALSE, TRUE)
    queries <- .checkUseDatasets(object, queries)
    if (any(queries %in% reference)) {
        if (showQueryWarn)
            warning("Reference dataset cannot be inclued in the query ",
                    "datasets. Removed from query list.")
        queries <- queries[!queries %in% reference]
    }
    if (isTRUE(verbose)) {
        .log("Imputing all the datasets exept the reference dataset\n",
             "Reference dataset: ", reference, "\n",
             "Query datasets: ", paste(queries, collapse = ", "))
    }

    referenceCells <- colnames(dataset(object, reference))
    for (i in seq_along(queries)) {
        query <- queries[i]
        queryLD <- dataset(object, query)
        queryCells <- colnames(queryLD)

        # creating a (reference cell numbers X query cell numbers) weights
        # matrix for knn weights and unit weights
        nn.k <- FNN::get.knnx(object@H.norm[referenceCells, ],
                              object@H.norm[queryCells, ],
                              k = nNeighbors,
                              algorithm = "CR")
        weights <- Matrix(0, nrow = length(referenceCells),
                          ncol = nrow(nn.k$nn.index), sparse = TRUE)
        if (isTRUE(weight)) {
            # for weighted situation
            # find nearest neighbors for query cell in normed ref datasets
            for (n in seq(nrow(nn.k$nn.index))) {
                # record ref-query cell-cell distances
                weights[nn.k$nn.index[n, ], n] <-
                    exp(-nn.k$nn.dist[n, ]) / sum(exp(-nn.k$nn.dist[n, ]))
            }
        } else{
            # for unweighted situation
            for (n in seq(nrow(nn.k$nn.index))) {
                # simply count the mean
                weights[nn.k$nn.index[n, ], n] <- 1 / nNeighbors
            }
        }

        # (genes by ref cell num) multiply by the weight matrix
        # (ref cell num by query cell num)
        referenceRawData <- getMatrix(object, "raw.data", dataset = reference)
        imputed_vals <- referenceRawData %*% weights
        # assigning dimnames
        colnames(imputed_vals) <- queryCells
        rownames(imputed_vals) <- rownames(referenceRawData)

        if (!inherits(imputed_vals, "dgCMatrix"))
            imputed_vals <- as(imputed_vals, "dgCMatrix")

        imputedLD <- createLigerDataset(raw.data = imputed_vals,
                                        modal = modalOf(object)[query], )
        datasets(object, check = FALSE)[[query]] <- imputedLD
    }

    if (isTRUE(norm))
        object <- normalize(object, ..., useDatasets = query, verbose = verbose)
    if (isTRUE(scale))
        object <- scaleNotCenter(object, useDatasets = query, verbose = verbose)

    return(object)
}
