#' Normalize raw data in liger object
#' @description Perform normalization on each raw data matrix from each dataset
#' to account for total feature expression across a cell.
#' @param object liger object
#' @param log Logical. Whether to do a \code{log(x + 1)} transform on the
#' normalized data. Default \code{TRUE}.
#' @param scaleFactor Numeric. Scale the normalized expression value by this
#' factor before transformation. \code{NULL} for not scaling. Default
#' \code{1e4}.
#' @param chunk Integer. Number of maximum number of cells in each chunk, when
#' normalization is applied to any HDF5 based dataset. Default \code{1000}.
#' @param verbose Logical. Whether to show information of the progress.
#' @return Updated \code{object}
#' @export
normalize <- function(
        object,
        log = TRUE,
        scaleFactor = 1e4,
        chunk = 1000,
        verbose = TRUE
) {
    if (scaleFactor <= 0 | scaleFactor == 1) scaleFactor <- NULL
    for (d in names(object)) {
        # `d` is the name of each dataset
        if (isTRUE(verbose)) message(date(), " ... Normalizing dataset: ", d)
        ld <- dataset(object, d)
        if (isH5Liger(ld)) {
            dataset(object, d, qc = FALSE) <- normalizeDataset.h5(
                ld, log, scaleFactor, chunk, verbose
            )
        } else {
            dataset(object, d, qc = FALSE) <- normalizeDataset.Matrix(
                ld, log, scaleFactor, verbose
            )
        }
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
    resultH5Path <- "norm.data3"
    # Use safe create here in practice
    safeH5Create(
        object = object,
        dataPath = resultH5Path,
        dims = raw.data(object)$dims,
        dtype = "double",
        chunkSize = chunkSize
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
    norm.data(object) <- h5file[[resultH5Path]]
    object@h5file.info$norm.data <- resultH5Path
    return(object)
}
