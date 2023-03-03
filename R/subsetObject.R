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
#' \code{"rawData"}, \code{"normData"} and \code{"scaleData"}. Default
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
        returnObject = TRUE,
        ...
) {
    if (inherits(object, "ligerDataset")) {
        object <- subsetLigerDataset(
            object = object, featureIdx = featureIdx, cellIdx = cellIdx,
            useSlot = useSlot, chunkSize = chunkSize, verbose = verbose,
            newH5 = newH5, returnObject = returnObject, ...)
        return(object)
    }
    if (!inherits(object, "liger")) {
        warning("`object` is not a liger obejct. Nothing to be done.")
        return(object)
    }
    # Check subscription parameters ####
    cellIdx <- .idxCheck(object, cellIdx, "cell")
    orderedCellIdx <- sort(cellIdx)
    barcodes <- colnames(object)[orderedCellIdx]
    datasetVar <- object$dataset[orderedCellIdx]
    useDatasets <- as.vector(unique(datasetVar))
    # feature idx need different check from ligerDataset's .idxCheck
    if (!is.null(featureIdx)) {
        if (!is.character(featureIdx)) {
            stop("Feature subscription from liger object can only take ",
                 "character vector.")
        }
        genesList <- lapply(datasets(object)[useDatasets], rownames)
        allGenes <- unique(unlist(genesList, use.names = FALSE))
        if (!all(featureIdx %in% allGenes)) {
            notFound <- featureIdx[!featureIdx %in% allGenes]
            warning(length(notFound), " out of ", length(featureIdx),
                    " given features were not found in the union of all ",
                    "features of used datasets")
        }
        featureIdx <- featureIdx[featureIdx %in% allGenes]
        if (length(featureIdx) == 0)
            stop("No feature can be retrieved")
    }
    # Subset each involved dataset and create new liger object

    datasets.new <- list()
    for (d in useDatasets) {
        if (isTRUE(verbose)) .log("Subsetting dataset: ",  d)
        ld <- dataset(object, d)
        featureIdxDataset <- featureIdx
        if (isFALSE(returnObject))
            featureIdxDataset <- featureIdx[featureIdx %in% rownames(ld)]
        ld <- subsetLigerDataset(
            object = ld, featureIdx = featureIdxDataset,
            cellIdx = barcodes[datasetVar == d], useSlot = useSlot,
            chunkSize = chunkSize, verbose = verbose, newH5 = newH5,
            returnObject = returnObject)
        datasets.new[[d]] <- ld
    }
    if (isTRUE(returnObject)) {
        if (!is.null(featureIdx)) {
            W <- object@W[featureIdx, , drop = FALSE]
            varFeature <- varFeatures(object)[varFeatures(object) %in%
                                                   featureIdx]
        } else {
            W <- object@W
            varFeature <- varFeatures(object)
        }
        return(methods::new(
            "liger",
            datasets = datasets.new,
            cellMeta = cellMeta(object, cellIdx = orderedCellIdx,
                                  drop = FALSE),
            varFeatures = varFeature,
            W = W,
            H.norm = object@H.norm[orderedCellIdx, , drop = FALSE],
            uns = object@uns,
            commands = object@commands
        ))
    } else {
        return(datasets.new)
    }
}

#' Retrieve a single matrix of cells from a slot
#' @description Only retrieve data from specific slot to reduce memory used by
#' a whole \linkS4class{liger} object of the subset. Useful for plotting.
#' Internally used by \code{\link{plotCellScatter}} and
#' \code{\link{plotCellViolin}}.
#' @param object \linkS4class{liger} object
#' @param feature Gene names, factor index or cell metadata variable names.
#' Should be available in specified \code{slot}.
#' @param slot Exactly choose from \code{"rawData"}, \code{"normData"},
#' \code{"scaleData"}, \code{"H"}, \code{"H.norm"} or \code{"cellMeta"}.
#' @param cellIdx Any valid type of index that subset from all cells. Default
#' \code{NULL} uses all cells.
#' @param ... Additional arguments passed to \code{\link{subsetLiger}} when
#' \code{slot} is one of \code{"rawData"}, \code{"normData"} or
#' \code{"scaleData"}.
#' @return A matrix object where rows are cells and columns are specified
#' features.
#' @export
retrieveCellFeature <- function(
        object,
        feature,
        slot = c("rawData", "normData", "scaleData",
                 "H", "H.norm", "cellMeta"),
        cellIdx = NULL,
        ...
) {
    slot <- match.arg(slot)
    cellIdx <- .idxCheck(object, idx = cellIdx, orient = "cell")
    if (is.null(cellIdx)) cellIdx <- seq(ncol(object))
    if (slot %in% c("rawData", "normData", "scaleData")) {
        subsetData <- subsetLiger(object, featureIdx = feature,
                                  cellIdx = cellIdx, useSlot = slot,
                                  returnObject = FALSE, ...)
        # value is expected to be a list structured like this
        # value
        # |--dataset1
        #    |--------slot1 subset matrix
        #    |--------slot2 subset matrix (Though multi-slot not supported here)
        # |--dataset2
        # ......
        value <- list()
        for (d in names(object)) {
            value[[d]] <- subsetData[[d]][[slot]]
        }
        # Condition for scaleData
        if (all(sapply(value, function(x) dim(x)[1] == 0)))
            stop("No feature could be retrieved. ",
                 "Please check feature names or slot")
        value <- lapply(value, as.matrix)
        value <- as.data.frame(t(mergeDenseAll(value)))
        orderedCellIdx <- sort(cellIdx)
        revIdx <- sapply(cellIdx, function(x) which(orderedCellIdx == x))
        value <- value[revIdx, , drop = FALSE]
    } else if (slot == "H") {
        value <- Reduce(cbind, getMatrix(object, "H"))
        value <- as.data.frame(t(value[feature, cellIdx, drop = FALSE]))
    } else if (slot == "H.norm") {
        value <- getMatrix(object, "H.norm")[cellIdx, feature, drop = FALSE]
    } else {
        value <- cellMeta(object, feature, cellIdx = cellIdx,
                           as.data.frame = TRUE, drop = FALSE)
    }
    return(value)
}

#' Subset ligerDataset object
#' @description This function subsets a \linkS4class{ligerDataset} object with
#' valid feature and cell indices. For HDF5 based object, options are available
#' for subsetting data into memory or a new on-disk H5 file. Feature and cell
#' subscription is always based on the size of rawData. Therefore, the feature
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
#' \code{"rawData"}, \code{"normData"} and \code{"scaleData"}. Default
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
    verbose = TRUE,
    returnObject = TRUE,
    ...
) {
    if (isH5Liger(object))
        subsetH5LigerDataset(object, featureIdx = featureIdx, cellIdx = cellIdx,
                             useSlot = useSlot, newH5 = newH5,
                             filename = filename,
                             filenameSuffix = filenameSuffix,
                             chunkSize = chunkSize, verbose = verbose,
                             returnObject = returnObject, ...)
    else subsetMemLigerDataset(object, featureIdx = featureIdx,
                               cellIdx = cellIdx, useSlot = useSlot,
                               returnObject = returnObject, ...)
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
    verbose = TRUE,
    returnObject = TRUE
) {
    if (newH5 == "auto") {
        cellIdx <- .idxCheck(object, cellIdx, "cell")
        if (length(cellIdx) > 8000) newH5 <- TRUE
        else newH5 <- FALSE
    }
    if (isTRUE(newH5)) {
        if (isTRUE(returnObject))
            warning("Cannot set `returnObject = FALSE` when subsetting",
                    "H5 based ligerDataset to new H5 file.")
        newObj <- subsetH5LigerDatasetToH5(
            object, filename = filename, cellIdx = cellIdx,
            featureIdx = featureIdx, filenameSuffix = filenameSuffix,
            useSlot = useSlot, chunkSize = chunkSize, verbose = verbose)
    } else if (isFALSE(newH5)) {
        newObj <- subsetH5LigerDatasetToMem(
            object, cellIdx = cellIdx, featureIdx = featureIdx,
            useSlot = useSlot, chunkSize = chunkSize, verbose = verbose,
            returnObject = returnObject)
    }
    newObj
}

subsetH5LigerDatasetToMem <- function(
    object,
    featureIdx = NULL,
    cellIdx = NULL,
    useSlot = NULL,
    returnObject = TRUE,
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
    modal <- modalOf(object)
    cellIdx <- .idxCheck(object, cellIdx, "cell")
    featureIdx <- .idxCheck(object, featureIdx, "feature")
    # Having `useSlot` as a record of user specification, to know whether it's
    # a `NULL` but involve everything, or just a few slots.
    slotInvolved <- .checkLDSlot(object, useSlot)
    value <- list()
    # Process rawData ####
    if ("rawData" %in% slotInvolved & !is.null(rawData(object))) {
        if (isTRUE(verbose)) .log("Subsetting `rawData`", level = 2)
        rawData <- H5Apply(
            object, init = NULL, useData = "rawData", chunkSize = chunkSize,
            verbose = verbose,
            FUN = function(chunk, sparseXIdx, chunkCellIdx, values) {
                subset <- chunk[featureIdx, chunkCellIdx %in% cellIdx, drop = FALSE]
                values <- cbind(values, subset)
            }
        )
        rownames(rawData) <- rownames(object)[featureIdx]
        colnames(rawData) <- colnames(object)[cellIdx]
        value$rawData <- rawData
    }

    # Process normData ####
    if ("normData" %in% slotInvolved & !is.null(normData(object))) {
        if (isTRUE(verbose)) .log("Subsetting `normData`", level = 2)
        normData <- H5Apply(
            object, init = NULL, useData = "normData", chunkSize = chunkSize,
            verbose = verbose,
            FUN = function(chunk, sparseXIdx, chunkCellIdx, values) {
                subset <- chunk[featureIdx, chunkCellIdx %in% cellIdx, drop = FALSE]
                values <- cbind(values, subset)
            }
        )
        rownames(normData) <- rownames(object)[featureIdx]
        colnames(normData) <- colnames(object)[cellIdx]
        value$normData <- normData
    }

    # Process scaled data ####
    if ("scaleData" %in% slotInvolved & !is.null(scaleData(object))) {
        if (isTRUE(verbose)) .log("Subsetting `scaleData`", level = 2)
        scaledFeatureIdx <- NULL
        if (getH5File(object)$exists("scaleData.featureIdx")) {
            scaledFeatureIdx <- getH5File(object)[["scaleData.featureIdx"]][]
        } else if ("selected" %in% names(featureMeta(object))) {
            scaledFeatureIdx <- which(featureMeta(object)$selected)
        } else {
            warning("Unable to know what features are included scaled data. ",
                    "Skipped.")
        }
        if (!is.null(scaledFeatureIdx) && length(scaledFeatureIdx) > 0) {
            if (length(scaledFeatureIdx) == scaleData(object)$dims[1]) {
                # Before: scaledFeatureIdx is based on all features, selecting
                # what features are in the scaleData
                # After: based on variable features, selecting what features
                # from variable features are being involved in featureIdx.
                scaledFeatureIdx2 <- unlist(lapply(featureIdx, function(i) {
                    which(scaledFeatureIdx == i)}), use.names = FALSE)
                if (isTRUE(verbose))
                    .log(length(scaledFeatureIdx2),
                         " features used in scaleData were selected. ",
                         level = 3)
                scaleData <- scaleData(object)[scaledFeatureIdx2, cellIdx,
                                                drop = FALSE]
                rownames(scaleData) <- rownames(object)[scaledFeatureIdx][scaledFeatureIdx2]
                colnames(scaleData) <- colnames(object)[cellIdx]
            } else {
                scaleData <- NULL
                warning("Row dimension of scaleData does not match with ",
                        "feature selection. Unable to subset from H5.")
            }
        }
        value$scaleData <- scaleData
    }
    # `NULL[idx1, idx2]` returns `NULL`
    # V: k x genes
    if (is.null(useSlot)) {
        if (exists("scaledFeatureIdx2")) sfi <- scaledFeatureIdx2
        else sfi <- NULL
        value$H <- object@H[, cellIdx, drop = FALSE]
        value$V <- object@V[sfi, , drop = FALSE]
        value$A <- object@A
        value$B <- object@B[sfi, , drop = FALSE]
        value$U <- object@U
        value$featureMeta <- featureMeta(object)[featureIdx, , drop = FALSE]
        # Additional subsetting for sub-classes, if applicable
        if (modal == "atac") {
            value$rawPeak <- rawPeak(object)[, cellIdx, drop = FALSE]
            value$normPeak <- normPeak(object)[, cellIdx, drop = FALSE]
        }
    }
    if (isTRUE(returnObject)) {
        value$modal <- modal
        do.call(createLigerDataset, value)
    } else {
        value
    }
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
    modal <- modalOf(object)
    cellIdx <- .idxCheck(object, cellIdx, "cell")
    featureIdx <- .idxCheck(object, featureIdx, "feature")
    useSlot <- .checkLDSlot(object, useSlot)
    if (is.null(filename) & is.null(filenameSuffix)) {
        oldFN <- h5fileInfo(object, "filename")
        filename <- paste0(oldFN, ".subset_",
                           format(Sys.time(), "%y%m%d_%H%M%S"),
                           ".h5")
    } else if (is.null(filename & !is.null(filenameSuffix))) {
        oldFN <- h5fileInfo(object, "filename")
        filename <- paste0(oldFN, ".", filenameSuffix, ".h5")
    }
    # Create new H5 file ####
    if (file.exists(filename)) {
        newH5File <- hdf5r::H5File$new(filename, mode = "r+")
    } else {
        newH5File <- hdf5r::H5File$new(filename, mode = "w")
    }
    if (isTRUE(verbose)) .log("New H5 file at: ", filename)
    newH5Meta <- h5fileInfo(object)
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
    if ("rawData" %in% useSlot & !is.null(rawData(object))) {
        # 1. Create paths to store i, p, x of sparse matrix
        if (isTRUE(verbose)) .log("Subsetting `rawData`", level = 2)
        safeH5Create(newH5File, newH5Meta$indices.name, dims = 1, dtype = "int")
        i.h5d <- newH5File[[newH5Meta$indices.name]]
        safeH5Create(newH5File, newH5Meta$indptr.name, dims = 1, dtype = "int")
        p.h5d <- newH5File[[newH5Meta$indptr.name]]
        safeH5Create(newH5File, newH5Meta$rawData, dims = 1, dtype = "double")
        x.h5d <- newH5File[[newH5Meta$rawData]]
        # 2. Go into chunks
        subsetSizes <- list(ix = 0, p = 1)
        H5Apply(
            object, init = subsetSizes, useData = "rawData",
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
    if ("normData" %in% useSlot & !is.null(normData(object))) {
        if (isTRUE(verbose)) .log("Subsetting `normData`", level = 2)
        safeH5Create(newH5File, newH5Meta$normData, dims = 1, dtype = "double")
        x.h5d <- newH5File[[newH5Meta$normData]]
        ipProcessedBefore <- exists("i.h5d") & exists("p.h5d")
        if (!ipProcessedBefore) {
            # H5D objects are created if have worked on rawData. Otherwise,
            # have to recreate them for normData. (normData and rawData share
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
            object, init = subsetSizes, useData = "normData",
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
    if ("scaleData" %in% useSlot & !is.null(scaleData(object))) {
        if (isTRUE(verbose)) .log("Subsetting `scaleData`", level = 2)
        scaledFeatureIdx <- NULL
        if (getH5File(object)$exists("scaleData.featureIdx")) {
            scaledFeatureIdx <- getH5File(object)[["scaleData.featureIdx"]][]
        } else if ("selected" %in% names(featureMeta(object))) {
            scaledFeatureIdx <- which(featureMeta(object)$selected)
        } else {
            warning("Unable to know what features are included scaled data. ",
                    "Skipped.")
        }
        if (!is.null(scaledFeatureIdx) && length(scaledFeatureIdx) > 0) {
            if (length(scaledFeatureIdx) == scaleData(object)$dims[1]) {
                scaledFeatureIdx2 <- unlist(lapply(featureIdx, function(i)
                    which(scaledFeatureIdx == i)), use.names = FALSE)
                ng <- length(scaledFeatureIdx2)
                nc <- length(cellIdx)
                if (isTRUE(verbose))
                    .log(length(scaledFeatureIdx2),
                         " features used in scaleData were selected. ",
                         level = 3)
                safeH5Create(newH5File, dataPath = newH5Meta$scaleData,
                             dims = c(ng, nc), dtype = "double",
                             chunkSize = c(ng, 1000))
                newH5File[[newH5Meta$scaleData]][1:ng,1:nc] <-
                    scaleData(object)[scaledFeatureIdx2, cellIdx, drop = FALSE]
                safeH5Create(newH5File, dataPath = "scaleData.featureIdx",
                             dims = ng, dtype = "int", chunkSize = ng)
                newH5File[["scaleData.featureIdx"]][1:ng] <-
                    .getOrderedSubsetIdx(featureIdx, scaledFeatureIdx)
            } else {
                warning("Row dimension of scaleData does not match with ",
                        "feature selection. Unable to subset from H5.")
            }
        }
    }
    newH5File$close()
    if (!"rawData" %in% useSlot) newH5Meta$rawData <- NULL
    if (!"normData" %in% useSlot) newH5Meta$normData <- NULL
    if (!"rawData" %in% useSlot & !"normData" %in% useSlot) {
        newH5Meta$indices.name <- NULL
        newH5Meta$indptr.name <- NULL
    }
    if (!"scaleData" %in% useSlot) newH5Meta$scaleData <- NULL
    # New object construction ####
    newObj <- createLigerDataset.h5(
        h5file = filename,
        rawData = newH5Meta$rawData, normData = newH5Meta$normData,
        scaleData = newH5Meta$scaleData,
        barcodes.name = newH5Meta$barcodes.name,
        genes.name = newH5Meta$genes.name,
        indices.name = newH5Meta$indices.name,
        indptr.name = newH5Meta$indptr.name,
        modal = modal,
        featureMeta = featureMeta(object)[featureIdx, , drop = FALSE]
    )
    newObj@H <- object@H[, cellIdx, drop = FALSE]
    newObj@V <- object@V[featureIdx, , drop = FALSE]
    newObj@A <- object@A
    newObj@B <- object@B[featureIdx, , drop = FALSE]
    newObj@U <- object@U

    newObj
}

#' @export
#' @rdname subsetLigerDataset
subsetMemLigerDataset <- function(object, featureIdx = NULL, cellIdx = NULL,
                                  useSlot = NULL, returnObject = TRUE) {
    if (!inherits(object, "ligerDataset")) {
        warning("`object` is not a ligerDataset obejct. Nothing to be done.")
        return(object)
    }
    if (isH5Liger(object)) {
        stop("`object` is HDF5 based. Use `subsetH5LigerDataset()` instead.")
    }
    if (is.null(cellIdx) && is.null(featureIdx)) return(object)
    modal <- modalOf(object)
    featureIdx <- .idxCheck(object, featureIdx, "feature")
    cellIdx <- .idxCheck(object, cellIdx, "cell")
    slotInvolved <- .checkLDSlot(object, useSlot)
    subsetData <- list()
    if ("rawData" %in% slotInvolved) {
        subsetData$rawData <- rawData(object)[featureIdx, cellIdx,
                                                drop = FALSE]
    }
    if ("normData" %in% slotInvolved) {
        subsetData$normData <- normData(object)[featureIdx, cellIdx,
                                                  drop = FALSE]
    }
    if ("scaleData" %in% slotInvolved) {
        if (!is.null(scaleData(object))) {
            scaleFeatureIdx <- rownames(scaleData(object)) %in%
                rownames(object)[featureIdx]
            subsetData$scaleData <-
                scaleData(object)[scaleFeatureIdx, cellIdx, drop = FALSE]
        }
        if (!is.null(object@scaleUnsharedData)) {
            scaleUnsFeatureIdx <- rownames(object@scaleUnsharedData) %in%
                rownames(object)[featureIdx]
            subsetData$scaleUnsharedData <-
                object@scaleUnsharedData[scaleUnsFeatureIdx, cellIdx,
                                           drop = FALSE]
        }
    }
    if (is.null(useSlot)) {
        # Users do not specify, subset the whole object
        if (exists("scaleFeatureIdx")) sfi <- scaleFeatureIdx
        else sfi <- seq(scaleData(object))
        subsetData <- c(subsetData,
                        list(H = object@H[, cellIdx, drop = FALSE],
                             V = object@V[sfi, , drop = FALSE],
                             A = object@A,
                             B = object@B[sfi, , drop = FALSE],
                             U = object@U,
                             featureMeta = object@featureMeta[featureIdx, ,
                                                                drop = FALSE]
                        ))
        # Additional subsetting for sub-classes, if applicable
        if (modal == "atac") {
            subsetData$rawPeak <- rawPeak(object)[, cellIdx, drop = FALSE]
            subsetData$normPeak <- normPeak(object)[, cellIdx, drop = FALSE]
        }
    }

    if (isTRUE(returnObject)) {
        subsetData$modal <- modal
        return(do.call("createLigerDataset", subsetData))
    }
    else return(subsetData)
}

.getOrderedSubsetIdx <- function(allNames, subsetNames) {
    # subsetNames must be real subset, but can be in a different order from
    # original allNames

    # Label the order of original allNames
    idx <- seq_along(allNames)
    names(idx) <- allNames
    # Subscribe with named vector, so the value (label for original order) get
    # ordered by subscription
    subsetIdx <- idx[subsetNames]
    subsetIdx <- subsetIdx[!is.na(subsetIdx)]
    names(subsetIdx) <- NULL
    subsetIdx
}

