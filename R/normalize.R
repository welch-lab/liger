#' Normalize raw data in liger object
#' @description Perform normalization on each raw data matrix from each dataset
#' to account for total feature expression across a cell.
#' @param object \linkS4class{liger} object
#' @param log Logical. Whether to do a \code{log(x + 1)} transform on the
#' normalized data. Default \code{TRUE}.
#' @param scaleFactor Numeric. Scale the normalized expression value by this
#' factor before transformation. \code{NULL} for not scaling. Default
#' \code{1e4}.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to be normalized. Should specify ATACseq
#' datasets when using \code{normalizePeak}. Default \code{NULL} normalizes all
#' valid datasets.
#' @param chunk Integer. Number of maximum number of cells in each chunk, when
#' normalization is applied to any HDF5 based dataset. Default \code{1000}.
#' @param verbose Logical. Whether to show information of the progress.
#' Default \code{TRUE}.
#' @return Updated \code{object}. When using \code{normalize()}, the
#' \code{norm.data} slot of each \linkS4class{ligerDataset} object in the
#' \code{datasets} slot is updated with normalized values calculated from
#' \code{raw.data} slot. While \code{normalizePeak()} normalizes raw peak counts
#' in \code{raw.peak} slot of a \linkS4class{ligerATACDataset} object and update
#' result in its \code{norm.peak} slot.
#' @rdname normalize
#' @export
normalize <- function(
        object,
        log = FALSE,
        scaleFactor = NULL,
        useDatasets = NULL,
        chunk = 1000,
        verbose = TRUE
) {
    useDatasets <- .checkUseDatasets(object, useDatasets)
    object <- recordCommand(object, dependencies = "hdf5r")
    if (!is.null(scaleFactor) && (scaleFactor <= 0 | scaleFactor == 1))
        scaleFactor <- NULL
    for (d in useDatasets) {
        # `d` is the name of each dataset
        if (isTRUE(verbose)) .log("Normalizing dataset: ", d)
        ld <- dataset(object, d)
        if (isH5Liger(ld))
            ld <- normalizeDataset.h5(ld, log, scaleFactor, chunk, verbose)
        else
            ld <- normalizeDataset.Matrix(ld, log, scaleFactor, verbose)
        datasets(object, check = FALSE)[[d]] <- ld
    }
    object
}

#' Perform normalization on ligerDataset object with in memory sparse matrix
#' @param object ligerDataset object
#' @param verbose Not used yet
#' @return Updated ligerDataset object
#' @noRd
normalizeDataset.Matrix <- function(object, log = TRUE, scaleFactor = 1e4,
                                    verbose = TRUE) {
    if (inherits(raw.data(object), "dgCMatrix") |
        inherits(raw.data(object), "dgTMatrix")) {
        norm.data(object, check = FALSE) <-
            Matrix.column_norm(raw.data(object))
    } else {
        norm.data(object, check = FALSE) <-
            sweep(x = raw.data(object), MARGIN = 2,
                  STATS = colSums(raw.data(object)), FUN = "/")
    }
    if (!is.null(scaleFactor))
        norm.data(object, check = FALSE) <- norm.data(object) * scaleFactor
    if (isTRUE(log))
        norm.data(object, check = FALSE) <- log1p(norm.data(object))
    object
}

#' Perform normalization on ligerDataset object with HDF5 file link
#' @param object ligerDataset object
#' @param chunkSize Integer for the maximum number of cells in each chunk
#' @param verbose Logical. Whether to show a progress bar.
#' @return Updated ligerDataset object
#' @noRd
normalizeDataset.h5 <- function(object, log = TRUE, scaleFactor = 1e4,
                                chunkSize = 1000, verbose = TRUE) {
    nFeature <- nrow(object)
    nCell <- ncol(object)

    # Initialize result
    results <- list(
        geneSumSq = rep(0, nFeature),
        geneMeans = rep(0, nFeature)
    )
    h5file <- getH5File(object)
    resultH5Path <- "norm.data"
    # Use safe create here in practice
    safeH5Create(
        object = object,
        dataPath = resultH5Path,
        dims = raw.data(object)$dims,
        dtype = "double",
        chunkSize = chunkSize
    )
    safeH5Create(
        object = object,
        dataPath = "gene_means",
        dims = nFeature,
        dtype = "double"
    )
    safeH5Create(
        object = object,
        dataPath = "gene_sum_sq",
        dims = nFeature,
        dtype = "double"
    )
    # Chunk run
    results <- H5Apply(
        object,
        function(chunk, sparseXIdx, cellIdx, values) {
            norm.data <- Matrix.column_norm(chunk)
            if (!is.null(scaleFactor)) norm.data <- norm.data * scaleFactor
            if (isTRUE(log)) norm.data <- log1p(norm.data)
            h5file[[resultH5Path]][sparseXIdx] <- norm.data@x
            row_sums <- rowSums(norm.data)
            values$geneSumSq <- values$geneSumSq +
                rowSums(norm.data * norm.data)
            values$geneMeans <- values$geneMeans + row_sums
            values
        },
        init = results,
        useData = "raw.data",
        chunkSize = chunkSize,
        verbose = verbose)
    results$geneMeans <- results$geneMeans / nCell
    h5file[["gene_means"]][seq_along(results$geneMeans)] <- results$geneMeans
    feature.meta(object, check = FALSE)$geneMeans <- results$geneMeans
    h5file[["gene_sum_sq"]][seq_along(results$geneSumSq)] <- results$geneSumSq
    feature.meta(object, check = FALSE)$geneSumSq <- results$geneSumSq
    norm.data(object, check = FALSE) <- h5file[[resultH5Path]]
    h5file.info(object, "norm.data", check = FALSE) <- resultH5Path
    return(object)
}

#' @rdname normalize
#' @export
normalizePeak <- function(
        object,
        log = FALSE,
        scaleFactor = NULL,
        useDatasets = NULL,
        verbose = TRUE
) {
    useDatasets <- .checkUseDatasets(object, useDatasets, modal = "atac")
    object <- recordCommand(object, dependencies = "hdf5r")
    if (!is.null(scaleFactor) && (scaleFactor <= 0 | scaleFactor == 1))
        scaleFactor <- NULL
    for (d in useDatasets) {
        # `d` is the name of each dataset
        if (isTRUE(verbose)) .log("Normalizing raw peak counts in dataset: ", d)
        ld <- dataset(object, d)
        if (inherits(raw.peak(ld), "dgCMatrix") |
            inherits(raw.peak(ld), "dgTMatrix"))
            norm.peak(ld, check = FALSE) <- Matrix.column_norm(raw.peak(ld))
        else
            norm.peak(ld, check = FALSE) <- sweep(x = raw.peak(ld), MARGIN = 2,
                                                  STATS = colSums(raw.peak(ld)),
                                                  FUN = "/")

        if (!is.null(scaleFactor))
            norm.peak(ld, check = FALSE) <- norm.peak(ld) * scaleFactor
        if (isTRUE(log))
            norm.peak(ld, check = FALSE) <- log1p(norm.peak(ld))
        datasets(object, check = FALSE)[[d]] <- ld
    }
    object
}

# Perform fast and memory-efficient normalization operation on sparse matrix data.
# param A Sparse matrix DGE.
Matrix.column_norm <- function(A) {
    if (class(A)[1] == "dgTMatrix") {
        temp = summary(A)
        A = Matrix::sparseMatrix(i = temp[, 1], j = temp[, 2], x = temp[, 3])
    }
    A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
    return(A)
}
