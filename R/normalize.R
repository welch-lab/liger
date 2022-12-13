#' Normalize raw data in liger object
#' @description Perform normalization on each raw data matrix from each dataset
#' to account for total feature expression across a cell.
#' @param object liger object
#' @param chunk Integer. Number of maximum number of cells in each chunk, when
#' normalization is applied to any HDF5 based dataset. Default \code{1000}.
#' @param verbose Logical. Whether to show information of the progress.
#' @return Updated \code{object}
#' @export
normalize <- function(object, chunk = 1000, verbose = TRUE) {
    for (d in names(object)) {
        # `d` is the name of each dataset
        if (isTRUE(verbose)) message(date(), " ... Normalizing dataset: ", d)
        ld <- dataset(object, d)
        if (isH5Liger(ld)) {
            dataset(object, d) <- normalizeDataset.h5(ld, chunk, verbose)
        } else {
            dataset(object, d) <- normalizeDataset.Matrix(ld, verbose)
        }
    }
    object
}

#' Perform normalization on ligerDataset object with in memory sparse matrix
#' @param object ligerDataset object
#' @param verbose Not used yet
#' @return Updated ligerDataset object
#' @noRd
normalizeDataset.Matrix <- function(object, verbose = TRUE) {
    if (inherits(raw.data(object), "dgCMatrix") |
        inherits(raw.data(object), "dgTMatrix")) {
        norm.data(object) <- Matrix.column_norm(raw.data(object))
    } else {
        norm.data(object) <- sweep(x = raw.data(object), MARGIN = 2,
                                   STATS = colSums(raw.data(object)),
                                   FUN = "/")
    }
    object
}

#' Perform normalization on ligerDataset object with HDF5 file link
#' @param object ligerDataset object
#' @param chunkSize Integer for the maximum number of cells in each chunk
#' @param verbose Logical. Whether to show a progress bar.
#' @return Updated ligerDataset object
#' @noRd
normalizeDataset.h5 <- function(object, chunkSize = 1000, verbose = TRUE) {
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
            #values$nUMI <- c(values$nUMI, colSums(chunk))
            norm.data <- Matrix.column_norm(chunk)
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
