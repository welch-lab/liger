#' Normalize raw data in liger object
#' @description Perform normalization on each raw data matrix from each dataset
#' to account for total feature expression across a cell.
#' @export
normalize <- function(object, chunk = 1000, verbose = TRUE) {
    datasets(object) <- lapply(datasets(object), function(ld) {
        if (isOnlineLiger(ld)) {
            return(normalizeDataset.online(ld, chunk, verbose))
        } else {
            return(normalizeDataset.Matrix(ld, verbose))
        }
    })
    object
}

#' Perform normalization on ligerDataset object with in memory sparse matrix
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
normalizeDataset.online <- function(object, chunkSize = 1000, verbose = TRUE) {
    nFeature <- nrow(object)
    nCell <- ncol(object)

    # Initialize result
    results <- list(
        geneSumSq = rep(0, nFeature),
        geneMeans = rep(0, nFeature)
    )
    h5.meta <- object@h5file.info
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
        h5file, init = results,
        function(chunk, sparseXIdx, values) {
            #values$nUMI <- c(values$nUMI, colSums(chunk))
            norm.data <- Matrix.column_norm(chunk)
            h5file[[resultH5Path]][sparseXIdx] <- norm.data@x
            row_sums <- rowSums(norm.data)
            values$geneSumSq <- values$geneSumSq +
                rowSums(norm.data * norm.data)
            values$geneMeans <- values$geneMeans + row_sums
            values
        },
        dataPath = h5.meta$raw.data,
        indicesPath = h5.meta$indices.name,
        indptrPath = h5.meta$indptr.name,
        cellIDPath = h5.meta$barcodes.name,
        featureIDPath = h5.meta$genes.name,
        chunkSize = chunkSize,
        verbose = verbose)
    results$geneMeans <- results$geneMeans / nCell
    norm.data(object) <- h5file[[resultH5Path]]
    return(object)
}
