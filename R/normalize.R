#' Normalize rawData in liger object
#' @description Perform normalization on each rawData matrix from each dataset
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
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return Updated \code{object}. When using \code{normalize()}, the
#' \code{normData} slot of each \linkS4class{ligerDataset} object in the
#' \code{datasets} slot is updated with normalized values calculated from
#' \code{rawData} slot. While \code{normalizePeak()} normalizes rawPeak counts
#' in \code{rawPeak} slot of a \linkS4class{ligerATACDataset} object and update
#' result in its \code{normPeak} slot.
#' @rdname normalize
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
normalize <- function(
        object,
        log = FALSE,
        scaleFactor = NULL,
        useDatasets = NULL,
        chunk = 1000,
        verbose = getOption("ligerVerbose")
) {
    .checkObjVersion(object)
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
                                    verbose = getOption("ligerVerbose")) {
    # Commented because only allowing dgCMatrix for rawData slot now
    #if (inherits(rawData(object), "dgCMatrix") |
    #    inherits(rawData(object), "dgTMatrix")) {
    normData(object, check = FALSE) <- Matrix.column_norm(rawData(object))
    #} else {
    #    normData(object, check = FALSE) <-
    #        sweep(x = rawData(object), MARGIN = 2,
    #              STATS = colSums(rawData(object)), FUN = "/")
    #}
    if (!is.null(scaleFactor))
        normData(object, check = FALSE) <- normData(object) * scaleFactor
    if (isTRUE(log))
        normData(object, check = FALSE) <- log1p(normData(object))
    object
}

#' Perform normalization on ligerDataset object with HDF5 file link
#' @param object ligerDataset object
#' @param chunkSize Integer for the maximum number of cells in each chunk
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return Updated ligerDataset object
#' @noRd
normalizeDataset.h5 <- function(object, log = TRUE, scaleFactor = 1e4,
                                chunkSize = 1000, verbose = getOption("ligerVerbose")) {
    nFeature <- nrow(object)
    nCell <- ncol(object)

    # Initialize result
    results <- list(
        geneSumSq = rep(0, nFeature),
        geneMeans = rep(0, nFeature)
    )
    h5file <- getH5File(object)
    resultH5Path <- "normData"
    # Use safe create here in practice
    safeH5Create(
        object = object,
        dataPath = resultH5Path,
        dims = rawData(object)$dims,
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
            normData <- Matrix.column_norm(chunk)
            if (!is.null(scaleFactor)) normData <- normData * scaleFactor
            if (isTRUE(log)) normData <- log1p(normData)
            h5file[[resultH5Path]][sparseXIdx] <- normData@x
            row_sums <- rowSums(normData)
            values$geneSumSq <- values$geneSumSq +
                rowSums(normData * normData)
            values$geneMeans <- values$geneMeans + row_sums
            values
        },
        init = results,
        useData = "rawData",
        chunkSize = chunkSize,
        verbose = verbose)
    results$geneMeans <- results$geneMeans / nCell
    h5file[["gene_means"]][seq_along(results$geneMeans)] <- results$geneMeans
    featureMeta(object, check = FALSE)$geneMeans <- results$geneMeans
    h5file[["gene_sum_sq"]][seq_along(results$geneSumSq)] <- results$geneSumSq
    featureMeta(object, check = FALSE)$geneSumSq <- results$geneSumSq
    normData(object, check = FALSE) <- h5file[[resultH5Path]]
    h5fileInfo(object, "normData", check = FALSE) <- resultH5Path
    return(object)
}

#' @rdname normalize
#' @export
normalizePeak <- function(
        object,
        log = FALSE,
        scaleFactor = NULL,
        useDatasets = NULL,
        verbose = getOption("ligerVerbose")
) {
    useDatasets <- .checkUseDatasets(object, useDatasets, modal = "atac")
    object <- recordCommand(object, dependencies = "hdf5r")
    if (!is.null(scaleFactor) && (scaleFactor <= 0 | scaleFactor == 1))
        scaleFactor <- NULL
    for (d in useDatasets) {
        # `d` is the name of each dataset
        if (isTRUE(verbose)) .log("Normalizing rawPeak counts in dataset: ", d)
        ld <- dataset(object, d)
        if (inherits(rawPeak(ld), "dgCMatrix")) {
            # |
            #inherits(rawPeak(ld), "dgTMatrix"))
            normed <- Matrix.column_norm(rawPeak(ld))
            rownames(normed) <- rownames(rawPeak(ld))

            normPeak(ld, check = FALSE) <- normed
        #else
        #    normPeak(ld, check = FALSE) <- sweep(x = rawPeak(ld), MARGIN = 2,
        #                                          STATS = colSums(rawPeak(ld)),
        #                                          FUN = "/")
        }
        if (!is.null(scaleFactor))
            normPeak(ld, check = FALSE) <- normPeak(ld) * scaleFactor
        if (isTRUE(log))
            normPeak(ld, check = FALSE) <- log1p(normPeak(ld))
        datasets(object, check = FALSE)[[d]] <- ld
    }
    object
}

# Perform fast and memory-efficient normalization operation on sparse matrix data.
# param A Sparse matrix DGE.
Matrix.column_norm <- function(A) {
    # Only dgCMatrix is allowed for rawData slot now
    # if (class(A)[1] == "dgTMatrix") {
    #     temp = summary(A)
    #     A = Matrix::sparseMatrix(i = temp[, 1], j = temp[, 2], x = temp[, 3])
    # }
    A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
    return(A)
}
