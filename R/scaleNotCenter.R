#' Scale genes by root-mean-square across cells
#' @description This function scales normalized gene expression data after
#' variable genes have been selected. Note that the data is not mean-centered
#' before scaling because expression values must remain positive (NMF only
#' accepts positive values).
#' @param object \linkS4class{liger} object, \code{link{selectGenes}} should
#' have been run in advance.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to be scaled but not centered. Default
#' \code{NULL} applies to all datasets.
#' @param chunk Integer. Number of maximum number of cells in each chunk, when
#' scaling is applied to any HDF5 based dataset. Default \code{1000}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{TRUE}.
#' @return Updated \code{object}, where the \code{scale.data} slot of each
#' \linkS4class{ligerDataset} object in the \code{datasets} slot is updated.
#' @export
scaleNotCenter <- function(
    object,
    useDatasets = NULL,
    chunk = 1000,
    verbose = TRUE
) {
    if (is.null(var.features(object)) ||
        length(var.features(object)) == 0) {
        stop("No variable feature found. ",
             "Please check the result of `selectGenes()`")
    }
    useDatasets <- .checkUseDatasets(object, useDatasets)
    object <- recordCommand(object, dependencies = "hdf5r")
    for (d in useDatasets) {
        if (isTRUE(verbose)) .log("Scaling dataset: ", d)
        ld <- dataset(object, d)
        if (isH5Liger(ld)) {
            # Scale H5 based data
            featureIdx <- .getOrderedSubsetIdx(rownames(ld),
                                               var.features(object))
            ld <- .scaleH5Matrix(ld, featureIdx = featureIdx,
                                 resultH5Path = "scale.data",
                                 chunk = chunk, verbose = verbose)

            if (!is.null(ld@var.unshared.features)) {
                unsharedIdx <- .getOrderedSubsetIdx(rownames(ld),
                                                    ld@var.unshared.features)
                ld <- .scaleH5Matrix(ld, featureIdx = unsharedIdx,
                                     resultH5Path = "scale.unshared.data",
                                     chunk = chunk, verbose = verbose)
            }
        } else {
            # Scale in memory data
            scaled <- .scaleMemMatrix(norm.data(ld)[var.features(object),])
            scale.data(ld, check = FALSE) <- scaled
            # For unshared var features
            if (!is.null(ld@var.unshared.features)) {
                scaled.unshared <- .scaleMemMatrix(
                    norm.data(ld)[ld@var.unshared.features,]
                )
                ld@scale.unshared.data <- scaled.unshared
            }
        }
        datasets(object, check = FALSE)[[d]] <- ld
    }
    object
}

.scaleH5Matrix <- function(ld, featureIdx, resultH5Path, chunk, verbose) {
    features <- rownames(ld)[featureIdx]
    geneSumSq <- feature.meta(ld)$geneSumSq[featureIdx]
    nCells <- ncol(ld)
    geneRootMeanSumSq = sqrt(geneSumSq / (nCells - 1))
    h5file <- getH5File(ld)
    safeH5Create(
        ld,
        dataPath = resultH5Path,
        dims = c(length(features), nCells),
        dtype = "double",
        chunkSize = c(length(features), chunk)
    )
    H5Apply(
        ld,
        useData = "norm.data",
        chunkSize = chunk,
        verbose = verbose,
        FUN = function(chunk, sparseXIdx, cellIdx, values) {
            chunk <- chunk[featureIdx, , drop = FALSE]
            chunk = as.matrix(chunk)
            chunk = sweep(chunk, 1, geneRootMeanSumSq, "/")
            rownames(chunk) <- features
            chunk[is.na(chunk)] = 0
            chunk[chunk == Inf] = 0
            h5file[[resultH5Path]][seq_along(features),
                                   cellIdx] <- chunk
        }
    )
    h5file.info(ld, resultH5Path, check = FALSE) <- resultH5Path
    safeH5Create(
        ld,
        dataPath = paste0(resultH5Path, ".featureIdx"),
        dims = length(features),
        dtype = "double"
    )
    h5file[[paste0(resultH5Path, ".featureIdx")]][1:length(featureIdx)] <-
        featureIdx
    return(ld)
}

# Input norm.subset should still be feature x cell dimensionality
.scaleMemMatrix <- function(norm.subset) {
    scaled <- scaleNotCenterFast(t(norm.subset))
    # scaled: g x c
    scaled <- as.matrix(t(scaled))
    scaled[is.na(scaled)] <- 0
    scaled[scaled == Inf] = 0
    rownames(scaled) <- rownames(norm.subset)
    colnames(scaled) <- colnames(norm.subset)
    return(scaled)
}
