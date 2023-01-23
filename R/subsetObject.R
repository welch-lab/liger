#' Subset liger object
#' @description This function subsets a \linkS4class{liger} object with
#' character feature index and any valid cell index. For datasets based on HDF5,
#' the filenames of subset H5 files could only be automatically generated for
#' now. Feature subsetting is based on the intersection of available features
#' from datasets involved by \code{cellIdx}, while \code{featureIdx = NULL} does
#' not take the intersection (i.e. nothing done on the feature axis).
#'
#' a \linkS4class{ligerDataset} object is also allowed for now and meanwhile,
#' setting \code{filename} is supported.
#' @param object A \linkS4class{liger} or \linkS4class{ligerDataset} object.
#' @param i,featureIdx Character vector. Missing or \code{NULL} for all
#' features.
#' @param j,cellIdx Character, logical or numeric index that can subscribe cells.
#' Missing or \code{NULL} for all cells.
#' @param useSlot The slot(s) to only consider. Choose one or more from
#' \code{"raw.data"}, \code{"norm.data"} and \code{"scale.data"}. Default
#' \code{NULL} subsets the whole object including analysis result matrices.
#' @param newH5 Whether to create new H5 files on disk for the subset datasets
#' if involved datasets in the \code{object} is HDF5 based. \code{TRUE} writes a
#' new ones, \code{FALSE} returns in memory data. Default \code{"auto"} writes a
#' new dataset when more than 8000 cells from that dataset is subscribed.
#' @param chunkSize Integer. Number of maximum number of cells in each chunk,
#' Default \code{1000}.
#' @param verbose Logical. Whether to show the progress. Default \code{TRUE}.
#' @param ... Arguments passed to \code{subsetLigerDataset}
#' @return Subset \code{object}
#' @export
#' @rdname subsetLiger
#' @seealso \code{\link{subsetLigerDataset}}
subsetLiger <- function(
        object,
        featureIdx = NULL,
        cellIdx = NULL,
        useSlot = NULL,
        chunkSize = 1000,
        verbose = TRUE,
        newH5 = "auto",
        ...
) {
    if (inherits(object, "ligerDataset")) {
        object <- subsetLigerDataset(
            object = object, featureIdx = featureIdx, cellIdx = cellIdx,
            useSlot = useSlot, chunkSize = chunkSize, verbose = verbose,
            newH5 = newH5, ...)
        return(object)
    }
    if (!inherits(object, "liger")) {
        warning("`object` is not a liger obejct. Nothing to be done.")
        return(object)
    }
    # Check subscription parameters ####
    cellIdx <- .idxCheck(object, cellIdx, "cell")
    barcodes <- colnames(object)[cellIdx]
    useDatasets <- as.vector(unique(object$dataset[cellIdx]))
    # feature idx need different check from ligerDataset's .idxCheck
    if (!is.null(featureIdx)) {
        if (!is.character(featureIdx)) {
            stop("Feature subscription from liger object can only take ",
                 "character vector.")
        }
        genesList <- lapply(datasets(object)[useDatasets], rownames)
        allGenes <- Reduce(intersect, genesList)
        if (!all(featureIdx %in% allGenes)) {
            notFound <- featureIdx[!featureIdx %in% allGenes]

            warning(length(notFound), " out of ", length(featureIdx),
                    " given features were not found in the intersection of the",
                    " features of used datasets")
        }
        featureIdx <- featureIdx[featureIdx %in% allGenes]
    }
    # Subset each involved dataset and create new liger object
    datasetVar <- object$dataset[cellIdx]
    datasets.new <- list()
    for (d in useDatasets) {
        if (isTRUE(verbose)) .log("Subsetting dataset: ",  d)
        ld <- dataset(object, d)
        ld <- subsetLigerDataset(
            object = ld, featureIdx = featureIdx,
            cellIdx = barcodes[datasetVar == d], useSlot = useSlot,
            chunkSize = chunkSize, verbose = verbose, newH5 = newH5)
        datasets.new[[d]] <- ld
    }
    methods::new(
        "liger",
        datasets = datasets.new,
        cell.meta = cell.meta(object)[cellIdx, , drop = FALSE],
        var.features = character(),
        W = object@W[, featureIdx, drop = FALSE],
        H.norm = object@H.norm[cellIdx, , drop = FALSE],
        version = packageVersion("rliger")
    )

}


#' Subset ligerDataset object
#' @description This function subsets a \linkS4class{ligerDataset} object with
#' valid feature and cell indices. For HDF5 based object, options are available
#' for subsetting data into memory or a new on-disk H5 file. Feature and cell
#' subscription is always based on the size of raw data. Therefore, the feature
#' subsetting on scaled data, which usually contains already a subset of
#' features, will select the intersection between the wanted features and the
#' set available from scaled data.
#' @param object \linkS4class{ligerDataset} object. HDF5 based object if using
#' \code{subsetH5LigerDataset}, in-memory data for \code{subsetMemLigerDataset}.
#' @param i,featureIdx Character, logical or numeric index that can subscribe
#' features. Missing or \code{NULL} for all features.
#' @param j,cellIdx Character, logical or numeric index that can subscribe cells.
#' Missing or \code{NULL} for all cells.
#' @param useSlot The slot(s) to only consider. Choose one or more from
#' \code{"raw.data"}, \code{"norm.data"} and \code{"scale.data"}. Default
#' \code{NULL} subsets the whole object including analysis result matrices.
#' @param newH5 Whether to create a new H5 file on disk for the subset dataset
#' if \code{object} is HDF5 based. \code{TRUE} writes a new one, \code{FALSE}
#' returns in memory data. Default \code{"auto"} writes a new one when more than
#' 8000 cells.
#' @param filename Filename of the new H5 file if being created. Default
#' \code{NULL} adds suffix \code{".subset_{yymmdd_HHMMSS}.h5"} to the original
#' name.
#' @param filenameSuffix Instead of specifying the exact filename, set a suffix
#' for the new files so the new filename looks like
#' \code{original.h5.[suffix].h5}. Default \code{NULL}.
#' @param chunkSize Integer. Number of maximum number of cells in each chunk,
#' Default \code{1000}.
#' @param verbose Logical. Whether to show the progress. Default \code{TRUE}.
#' @param ... Arguments passed to \code{subsetH5LigerDataset}
#' @return Subset \code{object}
#' @export
#' @rdname subsetLigerDataset
subsetLigerDataset <- function(
    object,
    featureIdx = NULL,
    cellIdx = NULL,
    useSlot = NULL,
    newH5 = "auto",
    filename = NULL,
    filenameSuffix = NULL,
    chunkSize = 1000,
    verbose = TRUE
) {
    if (isH5Liger(object))
        subsetH5LigerDataset(object, featureIdx = featureIdx, cellIdx = cellIdx,
                             useSlot = useSlot, newH5 = newH5,
                             filename = filename,
                             filenameSuffix = filenameSuffix,
                             chunkSize = chunkSize, verbose = verbose)
    else subsetMemLigerDataset(object, featureIdx = featureIdx,
                               cellIdx = cellIdx, useSlot = useSlot)
}

#' @export
#' @rdname subsetLigerDataset
subsetH5LigerDataset <- function(
    object,
    featureIdx = NULL,
    cellIdx = NULL,
    useSlot = NULL,
    newH5 = "auto",
    filename = NULL,
    filenameSuffix = NULL,
    chunkSize = 1000,
    verbose = TRUE
) {
    if (newH5 == "auto") {
        cellIdx <- .idxCheck(object, cellIdx, "cell")
        if (length(cellIdx) > 8000) newH5 <- TRUE
        else newH5 <- FALSE
    }
    if (isTRUE(newH5)) {
        newObj <- subsetH5LigerDatasetToH5(
            object, filename = filename, cellIdx = cellIdx,
            featureIdx = featureIdx, filenameSuffix = filenameSuffix,
            useSlot = useSlot, chunkSize = chunkSize, verbose = verbose)
    } else if (isFALSE(newH5)) {
        newObj <- subsetH5LigerDatasetToMem(
            object, cellIdx = cellIdx, featureIdx = featureIdx,
            useSlot = useSlot, chunkSize = chunkSize, verbose = verbose)
    }
    newObj
}

subsetH5LigerDatasetToMem <- function(
    object,
    featureIdx = NULL,
    cellIdx = NULL,
    useSlot = NULL,
    chunkSize = 1000,
    verbose = TRUE
) {
    if (!inherits(object, "ligerDataset")) {
        warning("`object` is not a ligerDataset obejct. Nothing to be done.")
        return(object)
    }
    if (!isH5Liger(object)) {
        warning("`object` is not HDF5 based. Nothing to be done.")
        return(object)
    }
    modal <- .classModalDict[[class(object)]]
    cellIdx <- .idxCheck(object, cellIdx, "cell")
    featureIdx <- .idxCheck(object, featureIdx, "feature")
    useSlot <- .checkLDSlot(object, useSlot)
    # Process raw data ####
    if ("raw.data" %in% useSlot & !is.null(raw.data(object))) {
        if (isTRUE(verbose)) .log("Subsetting `raw.data`", level = 2)
        rawData <- H5Apply(
            object, init = NULL, useData = "raw.data", chunkSize = chunkSize,
            verbose = verbose,
            FUN = function(chunk, sparseXIdx, chunkCellIdx, values) {
                subset <- chunk[featureIdx, chunkCellIdx %in% cellIdx]
                values <- cbind(values, subset)
            }
        )
        rownames(rawData) <- rownames(object)[featureIdx]
        colnames(rawData) <- colnames(object)[cellIdx]
    } else {
        rawData <- NULL
    }
    # Process norm.data ####
    if ("norm.data" %in% useSlot & !is.null(norm.data(object))) {
        if (isTRUE(verbose)) .log("Subsetting `norm.data`", level = 2)
        normData <- H5Apply(
            object, init = NULL, useData = "norm.data", chunkSize = chunkSize,
            verbose = verbose,
            FUN = function(chunk, sparseXIdx, chunkCellIdx, values) {
                subset <- chunk[featureIdx, chunkCellIdx %in% cellIdx]
                values <- cbind(values, subset)
            }
        )
        rownames(normData) <- rownames(object)[featureIdx]
        colnames(normData) <- colnames(object)[cellIdx]
    } else {
        normData <- NULL
    }
    # Process scaled data ####
    if ("scale.data" %in% useSlot & !is.null(scale.data(object))) {
        if (isTRUE(verbose)) .log("Subsetting `scale.data`", level = 2)
        scaledFeatureIdx <- NULL
        if (getH5File(object)$exists("scale.data.featureIdx")) {
            scaledFeatureIdx <- getH5File(object)[["scale.data.featureIdx"]][]
        } else if ("selected" %in% names(feature.meta(object))) {
            scaledFeatureIdx <- which(feature.meta(object)$selected)
        } else {
            warning("Unable to know what features are included scaled data. ",
                    "Skipped.")
        }
        if (!is.null(scaledFeatureIdx) && length(scaledFeatureIdx) > 0) {
            if (length(scaledFeatureIdx) == scale.data(object)$dims[1]) {
                # Before: scaledFeatureIdx is based on all features, selecting
                # what features are in the scale.data
                # After: based on variable features, selecting what features
                # from variable features are being involved in featureIdx.
                scaledFeatureIdx2 <- which(scaledFeatureIdx %in% featureIdx)
                if (isTRUE(verbose))
                    .log(length(scaledFeatureIdx2),
                         " features used in scale.data were selected. ",
                         level = 3)
                scaleData <- scale.data(object)[scaledFeatureIdx2, cellIdx]
                rownames(scaleData) <- rownames(object)[scaledFeatureIdx2]
                colnames(scaleData) <- colnames(object)[cellIdx]
            } else {
                scaleData <- NULL
                warning("Row dimension of scale.data does not match with ",
                        "feature selection. Unable to subset from H5.")
            }
        }
    } else {
        scaleData <- NULL
    }
    # `NULL[idx1, idx2]` returns `NULL`
    # V: k x genes
    H <- object@H[cellIdx, , drop = FALSE]
    V <- object@V[, featureIdx, drop = FALSE]
    A <- object@A
    B <- object@B[featureIdx, , drop = FALSE]
    U <- object@U#[cellIdx, , drop = FALSE]
    # New object construction ####
    createLigerDataset(raw.data = rawData, modal = modal, norm.data = normData,
                       scale.data = scaleData, H = H, V = V, A = A, B = B,
                       U = U,
                       feature.meta = feature.meta(object)[featureIdx, ,
                                                           drop = FALSE])
}

subsetH5LigerDatasetToH5 <- function(
        object,
        featureIdx = NULL,
        cellIdx = NULL,
        useSlot = NULL,
        filename = NULL,
        filenameSuffix = NULL,
        chunkSize = 1000,
        verbose = TRUE
) {
    # Input checks ####
    if (!inherits(object, "ligerDataset")) {
        warning("`object` is not a ligerDataset obejct. Nothing to be done.")
        return(object)
    }
    if (!isH5Liger(object)) {
        warning("`object` is not HDF5 based. Nothing to be done.")
        return(object)
    }
    modal <- .classModalDict[[class(object)]]
    cellIdx <- .idxCheck(object, cellIdx, "cell")
    featureIdx <- .idxCheck(object, featureIdx, "feature")
    useSlot <- .checkLDSlot(object, useSlot)
    if (is.null(filename) & is.null(filenameSuffix)) {
        oldFN <- h5file.info(object, "filename")
        filename <- paste0(oldFN, ".subset_",
                           format(Sys.time(), "%y%m%d_%H%M%S"),
                           ".h5")
    } else if (is.null(filename & !is.null(filenameSuffix))) {
        oldFN <- h5file.info(object, "filename")
        filename <- paste0(oldFN, ".", filenameSuffix, ".h5")
    }
    # Create new H5 file ####
    if (file.exists(filename)) {
        newH5File <- hdf5r::H5File$new(filename, mode = "r+")
    } else {
        newH5File <- hdf5r::H5File$new(filename, mode = "w")
    }
    if (isTRUE(verbose)) .log("New H5 file at: ", filename)
    newH5Meta <- h5file.info(object)
    newH5Meta$H5File <- newH5File
    newH5Meta$filename <- filename
    if (newH5Meta$format.type != "AnnData") {
        safeH5Create(newH5File, newH5Meta$barcodes.name,
                     dims = length(cellIdx), dtype = "char")
        newH5File[[newH5Meta$barcodes.name]][1:length(cellIdx)] <-
            colnames(object)[cellIdx]
    } else {
        # TODO: AnnData style barcodes storage.
    }

    safeH5Create(newH5File, newH5Meta$genes.name,
                 dims = length(featureIdx), dtype = "char")
    newH5File[[newH5Meta$genes.name]][1:length(featureIdx)] <-
        rownames(object)[featureIdx]
    # Process Raw Data ####
    if ("raw.data" %in% useSlot & !is.null(raw.data(object))) {
        # 1. Create paths to store i, p, x of sparse matrix
        if (isTRUE(verbose)) .log("Subsetting `raw.data`", level = 2)
        safeH5Create(newH5File, newH5Meta$indices.name, dims = 1, dtype = "int")
        i.h5d <- newH5File[[newH5Meta$indices.name]]
        safeH5Create(newH5File, newH5Meta$indptr.name, dims = 1, dtype = "int")
        p.h5d <- newH5File[[newH5Meta$indptr.name]]
        safeH5Create(newH5File, newH5Meta$raw.data, dims = 1, dtype = "double")
        x.h5d <- newH5File[[newH5Meta$raw.data]]
        # 2. Go into chunks
        subsetSizes <- list(ix = 0, p = 1)
        H5Apply(
            object, init = subsetSizes, useData = "raw.data",
            chunkSize = chunkSize, verbose = verbose,
            FUN = function(chunk, sparseXIdx, chunkCellIdx, values) {
                if (sum(chunkCellIdx %in% cellIdx) > 0) {
                    subset <- chunk[featureIdx, chunkCellIdx %in% cellIdx]
                    # Length of `subset@i` and `subset@x` should always be the
                    # same.
                    #
                    # Here, `values$ix` & `values$p` store the sizes of
                    # currently existing value in H5 file. `ix.new` and `p.new`
                    # is the new size after a round of concatenation. These are
                    # all single numbers.
                    #
                    # `i`, `p` and `x` are the constructor vectors of the subset
                    # sparse matrix.
                    i <- subset@i
                    ix.new <- values$ix + length(i)
                    hdf5r::extendDataSet(i.h5d, ix.new)
                    i.h5d[(values$ix + 1):ix.new] <- i
                    # Exactly same way as processing `i` above
                    x <- subset@x
                    hdf5r::extendDataSet(x.h5d, ix.new)
                    x.h5d[(values$ix + 1):ix.new] <- x
                    # Different way for processing `p`
                    p <- subset@p
                    p.new <- values$p + length(p) - 1
                    hdf5r::extendDataSet(p.h5d, p.new)
                    p.h5d[(values$p + 1):p.new] <-
                        p[2:length(p)] + p.h5d[values$p]
                    values$ix <- ix.new
                    values$p <- p.new
                }
                return(values)
            }
        )
    }
    # Process Normalized Data ####
    if ("norm.data" %in% useSlot & !is.null(norm.data(object))) {
        if (isTRUE(verbose)) .log("Subsetting `norm.data`", level = 2)
        safeH5Create(newH5File, newH5Meta$norm.data, dims = 1, dtype = "double")
        x.h5d <- newH5File[[newH5Meta$norm.data]]
        ipProcessedBefore <- exists("i.h5d") & exists("p.h5d")
        if (!ipProcessedBefore) {
            # H5D objects are created if have worked on raw.data. Otherwise,
            # have to recreate them for norm.data. (norm.data and raw.data share
            # the same `i` and `p` vectors as in the form of sparse matrices.)
            safeH5Create(newH5File, newH5Meta$indices.name,
                         dims = 1, dtype = "int")
            i.h5d <- newH5File[[newH5Meta$indices.name]]
            safeH5Create(newH5File, newH5Meta$indptr.name,
                         dims = 1, dtype = "int")
            p.h5d <- newH5File[[newH5Meta$indptr.name]]
        }
        subsetSizes <- list(ix = 0, p = 1)
        H5Apply(
            object, init = subsetSizes, useData = "norm.data",
            chunkSize = chunkSize, verbose = verbose,
            FUN = function(chunk, sparseXIdx, chunkCellIdx, values) {
                if (sum(chunkCellIdx %in% cellIdx) > 0) {
                    subset <- chunk[featureIdx, chunkCellIdx %in% cellIdx]
                    x <- subset@x
                    ix.new <- values$ix + length(x)
                    hdf5r::extendDataSet(x.h5d, ix.new)
                    x.h5d[(values$ix + 1):ix.new] <- x
                    if (!ipProcessedBefore) {
                        i <- subset@i
                        hdf5r::extendDataSet(i.h5d, ix.new)
                        i.h5d[(values$ix + 1):ix.new] <- i
                        p <- subset@p
                        p.new <- values$p + length(p) - 1
                        hdf5r::extendDataSet(p.h5d, p.new)
                        p.h5d[(values$p + 1):p.new] <-
                            p[2:length(p)] + p.h5d[values$p]
                        values$p <- p.new
                    }
                    values$ix <- ix.new
                }
                return(values)
            }
        )
    }
    # Process Scaled Data ####
    if ("scale.data" %in% useSlot & !is.null(scale.data(object))) {
        if (isTRUE(verbose)) .log("Subsetting `scale.data`", level = 2)
        scaledFeatureIdx <- NULL
        if (getH5File(object)$exists("scale.data.featureIdx")) {
            scaledFeatureIdx <- getH5File(object)[["scale.data.featureIdx"]][]
        } else if ("selected" %in% names(feature.meta(object))) {
            scaledFeatureIdx <- which(feature.meta(object)$selected)
        } else {
            warning("Unable to know what features are included scaled data. ",
                    "Skipped.")
        }
        if (!is.null(scaledFeatureIdx) && length(scaledFeatureIdx) > 0) {
            if (length(scaledFeatureIdx) == scale.data(object)$dims[1]) {
                scaledFeatureIdx2 <- which(scaledFeatureIdx %in% featureIdx)
                ng <- length(scaledFeatureIdx2)
                nc <- length(cellIdx)
                if (isTRUE(verbose))
                    .log(length(scaledFeatureIdx2),
                         " features used in scale.data were selected. ",
                         level = 3)
                safeH5Create(newH5File, dataPath = newH5Meta$scale.data,
                             dims = c(ng, nc), dtype = "double",
                             chunkSize = c(ng, nc))
                newH5File[[newH5Meta$scale.data]][1:ng,1:nc] <-
                    scale.data(object)[scaledFeatureIdx2, cellIdx]
                safeH5Create(newH5File, dataPath = "scale.data.featureIdx",
                             dims = ng, dtype = "int", chunkSize = ng)
                newH5File[["scale.data.featureIdx"]][1:ng] <-
                    which(featureIdx %in% scaledFeatureIdx)
            } else {
                warning("Row dimension of scale.data does not match with ",
                        "feature selection. Unable to subset from H5.")
            }
        }
    }
    newH5File$close()
    if (!"raw.data" %in% useSlot) newH5Meta$raw.data <- NULL
    if (!"norm.data" %in% useSlot) newH5Meta$norm.data <- NULL
    if (!"raw.data" %in% useSlot & !"norm.data" %in% slot) {
        newH5Meta$indices.name <- NULL
        newH5Meta$indptr.name <- NULL
    }
    if (!"scale.data" %in% useSlot) newH5Meta$scale.data <- NULL
    # New object construction ####
    newObj <- createLigerDataset.h5(
        h5file = filename,
        raw.data = newH5Meta$raw.data, norm.data = newH5Meta$norm.data,
        scale.data = newH5Meta$scale.data,
        barcodes.name = newH5Meta$barcodes.name,
        genes.name = newH5Meta$genes.name,
        indices.name = newH5Meta$indices.name,
        indptr.name = newH5Meta$indptr.name,
        modal = modal,
        feature.meta = feature.meta(object)[featureIdx, , drop = FALSE]
    )
    newObj@H <- object@H[cellIdx, , drop = FALSE]
    newObj@V <- object@V[, featureIdx, drop = FALSE]
    newObj@A <- object@A
    newObj@B <- object@B[featureIdx, , drop = FALSE]
    newObj@U <- object@U

    newObj
}

#' @export
#' @rdname subsetLigerDataset
subsetMemLigerDataset <- function(object, featureIdx = NULL, cellIdx = NULL,
                                  useSlot = NULL) {
    if (!inherits(object, "ligerDataset")) {
        warning("`object` is not a ligerDataset obejct. Nothing to be done.")
        return(object)
    }
    if (isH5Liger(object)) {
        stop("`object` is HDF5 based. Use `subsetH5LigerDataset()` instead.")
    }
    if (is.null(cellIdx) && is.null(featureIdx)) return(object)
    modal <- .classModalDict[[class(object)]]
    featureIdx <- .idxCheck(object, featureIdx, "feature")
    cellIdx <- .idxCheck(object, cellIdx, "cell")
    slot.use <- .checkLDSlot(object, useSlot)
    subsetData <- list(modal = modal)
    if ("raw.data" %in% slot.use) {
        subsetData$raw.data <- raw.data(object)[featureIdx, cellIdx,
                                                drop = FALSE]
    }
    if ("norm.data" %in% slot.use) {
        subsetData$norm.data <- norm.data(object)[featureIdx, cellIdx,
                                                  drop = FALSE]
    }
    if ("scale.data" %in% slot.use) {
        if (!is.null(scale.data(object))) {
            scaleFeatureIdx <-
                which(rownames(object) %in% rownames(scale.data(object)))
            subsetData$scale.data <-
                scale.data(object)[scaleFeatureIdx %in% featureIdx, cellIdx,
                                   drop = FALSE]
        }
        if (!is.null(object@scale.unshared.data)) {
            scaleUnsFeatureIdx <-
                which(rownames(object) %in%
                          rownames(object@scale.unshared.data))
            subsetData$scale.unshared.data <-
                object@scale.unshared.data[scaleUnsFeatureIdx %in% featureIdx,
                                           cellIdx,
                                           drop = FALSE]
        }
    }
    if (is.null(useSlot)) {
        # Users do not specify, subset the whole object
        subsetData <- c(subsetData,
                        list(H = object@H[cellIdx, , drop = FALSE],
                             V = object@V[, featureIdx, drop = FALSE],
                             A = object@A,
                             B = object@B[featureIdx, , drop = FALSE],
                             U = object@U,
                             feature.meta = object@feature.meta[featureIdx, ,
                                                                drop = FALSE]
                        ))
        # Additional subsetting for sub-classes, if applicable
        if (modal == "atac") {
            subsetData$peak <- peak(object)[, cellIdx, drop = FALSE]
        }
    }

    do.call("createLigerDataset", subsetData)
}
