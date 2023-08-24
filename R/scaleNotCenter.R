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
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return Updated \code{object}, where the \code{scaleData} slot of each
#' \linkS4class{ligerDataset} object in the \code{datasets} slot is updated.
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
scaleNotCenter <- function(
    object,
    useDatasets = NULL,
    chunk = 1000,
    verbose = getOption("ligerVerbose")
) {
    if (inherits(object, "dgCMatrix")) {
        return(.scaleMemMatrix(object))
    }
    .checkObjVersion(object)
    if (is.null(varFeatures(object)) ||
        length(varFeatures(object)) == 0) {
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
                                               varFeatures(object))
            # ld <- .scaleH5Matrix(ld, featureIdx = featureIdx,
            #                      resultH5Path = "scaleData",
            #                      chunk = chunk, verbose = verbose)
            ld <- .scaleH5SpMatrix(ld, featureIdx = featureIdx,
                                   resultH5Path = "scaleDataSparse",
                                   chunk = chunk, verbose = verbose)

            if (!is.null(ld@varUnsharedFeatures) &&
                length(ld@varUnsharedFeatures) > 0) {
                unsharedIdx <- .getOrderedSubsetIdx(rownames(ld),
                                                    ld@varUnsharedFeatures)
                # ld <- .scaleH5Matrix(ld, featureIdx = unsharedIdx,
                #                      resultH5Path = "scaleUnsharedData",
                #                      chunk = chunk, verbose = verbose)
                ld <- .scaleH5SpMatrix(ld, featureIdx = unsharedIdx,
                                       resultH5Path = "scaleUnsharedDataSparse",
                                       chunk = chunk, verbose = verbose)
            }
        } else {
            # Scale in memory data
            scaled <- .scaleMemMatrix(normData(ld)[varFeatures(object), ,
                                                   drop = FALSE])
            scaleData(ld, check = FALSE) <- scaled
            # For unshared var features
            if (length(ld@varUnsharedFeatures) > 0) {
                scaled.unshared <- .scaleMemMatrix(
                    normData(ld)[ld@varUnsharedFeatures, ,drop = FALSE]
                )
                ld@scaleUnsharedData <- scaled.unshared
            }
        }
        datasets(object, check = FALSE)[[d]] <- ld
    }
    object
}

.scaleH5SpMatrix <- function(ld, featureIdx, resultH5Path, chunk, verbose) {
    features <- rownames(ld)[featureIdx]
    geneSumSq <- featureMeta(ld)$geneSumSq[featureIdx]
    nCells <- ncol(ld)
    geneRootMeanSumSq = sqrt(geneSumSq / (nCells - 1))
    h5file <- getH5File(ld)
    # Count the subset nnz first before creating data space
    nnz <- 0
    H5Apply(
        ld, useData = "normData", chunkSize = chunk, verbose = FALSE,
        FUN = function(chunk, sparseXIdx, cellIdx, values) {
            chunk <- chunk[featureIdx, , drop = FALSE]
            nnz <- nnz + length(chunk@x)
        }
    )
    # Create datasets
    dataPath <- paste0(resultH5Path, "/data")
    rowindPath <- paste0(resultH5Path, "/indices")
    colptrPath <- paste0(resultH5Path, "/indptr")
    safeH5Create(ld, dataPath = dataPath, dims = nnz,
                 dtype = "double", chunkSize = 2048)
    safeH5Create(ld, dataPath = rowindPath, dims = nnz,
        dtype = "int", chunkSize = 2048)
    safeH5Create(ld, dataPath = colptrPath, dims = nCells + 1,
        dtype = "int", chunkSize = 1024)
    # Process chunks of sparse normData, and write to sparse scaleData
    h5file[[colptrPath]][1] <- 0
    H5Apply(
        ld,
        useData = "normData",
        init = c(1, 0), # [1] record of nnz idx start [2] record of last colptr
        chunkSize = chunk,
        verbose = verbose,
        FUN = function(chunk, sparseXIdx, cellIdx, values) {
            # Subset variable features
            chunk <- chunk[featureIdx, , drop = FALSE]
            # Calculate scale not center
            chunk <- sweep(chunk, 1, geneRootMeanSumSq, "/")
            chunk[is.na(chunk)] = 0
            chunk[chunk == Inf] = 0
            # Make sure of the sparse format, could be skipped
            chunk <- methods::as(chunk, "CsparseMatrix")
            # Write
            nnzRange <- seq(from = values[1], length.out = length(chunk@i))
            h5file[[rowindPath]][nnzRange] <- chunk@i
            h5file[[dataPath]][nnzRange] <- chunk@x
            values[1] <- values[1] + length(nnzRange)
            increColptr <- chunk@p + values[2]
            h5file[[colptrPath]][cellIdx + 1] =
                increColptr[seq(2, length(increColptr))]
            values[2] <- increColptr[length(increColptr)]
            return(values)
        }
    )
    safeH5Create(ld, dataPath = paste0(resultH5Path, "/featureIdx"),
                 dims = length(features), dtype = "int")
    h5file[[paste0(resultH5Path, "/featureIdx")]][1:length(featureIdx)] <-
        featureIdx
    return(ld)
}

.scaleH5Matrix <- function(ld, featureIdx, resultH5Path, chunk, verbose) {
    features <- rownames(ld)[featureIdx]
    geneSumSq <- featureMeta(ld)$geneSumSq[featureIdx]
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
        useData = "normData",
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
    h5fileInfo(ld, resultH5Path, check = FALSE) <- resultH5Path
    safeH5Create(
        ld,
        dataPath = paste0(resultH5Path, ".featureIdx"),
        dims = length(features),
        dtype = "int"
    )
    h5file[[paste0(resultH5Path, ".featureIdx")]][1:length(featureIdx)] <-
        featureIdx
    return(ld)
}

# Input norm.subset should still be feature x cell dimensionality
.scaleMemMatrix <- function(norm.subset) {
    if (nrow(norm.subset) == 0)
        return(matrix(nrow = 0, ncol = ncol(norm.subset),
                      dimnames = list(NULL, colnames(norm.subset))))
    scaled <- t(scaleNotCenterFast(t(norm.subset)))
    # scaled: g x c
    # scaled <- as.matrix(t(scaled))
    scaled[is.na(scaled)] <- 0
    scaled[scaled == Inf] = 0
    rownames(scaled) <- rownames(norm.subset)
    colnames(scaled) <- colnames(norm.subset)
    return(scaled)
}

#' Create "scaled data" for DNA methylation datasets
#' @description
#' Because gene body mCH proportions are negatively correlated with gene
#' expression level in neurons, we need to reverse the direction of the
#' methylation data. We do this by simply subtracting all values from the
#' maximum methylation value. The resulting values are positively correlated
#' with gene expression. This will only be applied to variable genes detected in
#' prior.
#' @param object A \linkS4class{liger} object, with variable genes identified.
#' @param useDatasets Required. A character vector of the names, a numeric or
#' logical vector of the index of the datasets that should be identified as
#' methylation data where the reversed data will be created.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return The input \linkS4class{liger} object, where the \code{scaleData} slot
#' of the specified datasets will be updated with value as described in
#' Description.
#' @export
#' @examples
#' # Assuming the second dataset in example data "pbmc" is methylation data
#' pbmc <- normalize(pbmc, useDatasets = 1)
#' pbmc <- selectGenes(pbmc, datasets.use = 1)
#' pbmc <- scaleNotCenter(pbmc, useDatasets = 1)
#' pbmc <- reverseMethData(pbmc, useDatasets = 2)
reverseMethData <- function(object, useDatasets,
                            verbose = getOption("ligerVerbose")) {
    useDatasets <- .checkUseDatasets(object, useDatasets)
    if (is.null(varFeatures(object)) || length(varFeatures(object)) == 0) {
        stop("Variable genes have to be identified first. ",
             "Please run `selectGenes(object)`.")
    }
    for (d in useDatasets) {
        ld <- dataset(object, d)
        raw <- rawData(ld)
        if (isTRUE(verbose)) .log("Substracting methylation data: ", d)
        scaleData(ld, check = FALSE) <- max(raw) -
            as.matrix(raw[varFeatures(object),])
        datasets(object, check = FALSE)[[d]] <- ld
    }
    return(object)
}
