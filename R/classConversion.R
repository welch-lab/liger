setClass("ligerDataset")
setClassUnion("matrixLike", c("matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# From other things to liger class ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Converting other classes of data to a liger object
#' @export
#' @param object Object.
#' @param datasetVar If one or more datasets exist in the given object, specify
#' a name for the single dataset; select a variable from existing metadata (e.g.
#' colData column); specify a vector/factor that assign the dataset belonging.
#' Default \code{NULL} gathers things into one dataset and names it "sample".
#' @param modal Modality setting for each dataset.
#' @param ... Additional arguments passed to \code{\link{createLiger}}
#' @return a \linkS4class{liger} object.
#' @rdname as.liger
as.liger <- function(object, ...) UseMethod("as.liger", object)

#' @rdname as.liger
#' @export
#' @method as.liger dgCMatrix
as.liger.dgCMatrix <- function(
        object,
        datasetVar = NULL,
        modal = NULL,
        ...
) {
    rawDataList <- list(sample = object)
    if (!is.null(datasetVar)) {
        if (length(datasetVar) == 1) names(rawDataList) <- datasetVar
        else if (length(datasetVar) == ncol(rawDataList[[1]])) {
            if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
            datasetVar <- droplevels(datasetVar)
            rawDataList <- lapply(levels(datasetVar), function(var) {
                rawDataList[[1]][, datasetVar == var, drop = FALSE]
            })
            names(rawDataList) <- levels(datasetVar)
        } else {
            stop("Unable to understand `datasetVar`")
        }
    } else {
        datasetVar <- "sample"
    }
    modal <- .checkArgLen(modal, length(datasetVar))
    createLiger(rawData = rawDataList, modal = modal, ...)
}

#' @rdname as.liger
#' @export
#' @method as.liger SingleCellExperiment
as.liger.SingleCellExperiment <- function(
        object,
        datasetVar = NULL,
        modal = NULL,
        ...
) {
    if (!requireNamespace("SingleCellExperiment", quietly = "TRUE"))
        stop("Package \"SingleCellExperiment\" needed for this function ",
             "to work. Please install it by command:\n",
             "BiocManager::install('SingleCellExperiment')",
             call. = FALSE)
    rawDataList <- list(SingleCellExperiment::counts(object))
    sampleCol <- NULL
    # setNames <- function(x, n) {
    #     if (!is.null(x) && length(x) > 0) names(x) <- n
    #     return(x)
    # }
    if (is.null(datasetVar)) {
        # One dataset, no name so by default sce
        names(rawDataList) <- "sce"
        datasetVar <- factor("sce")
        # rawDataList <- setNames(rawDataList, "sce")
    } else {
        if (length(datasetVar) == 1) {
            if (datasetVar %in%
                colnames(SummarizedExperiment::colData(object))) {
                sampleCol <- datasetVar
                # Split by variable in colData
                datasetVar <- SummarizedExperiment::colData(object)[[datasetVar]]
                if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
                datasetVar <- droplevels(datasetVar)
                rawDataList <- lapply(levels(datasetVar), function(var) {
                    rawDataList[[1]][, datasetVar == var, drop = FALSE]
                })
                names(rawDataList) <- levels(datasetVar)
            } else {
                # One dataset, use given name
                names(rawDataList) <- datasetVar
                datasetVar <- factor(datasetVar)
            }
        } else if (length(datasetVar) == ncol(rawDataList[[1]])) {
            # Split by user given variable
            if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
            datasetVar <- droplevels(datasetVar)
            rawDataList <- lapply(levels(datasetVar), function(var) {
                rawDataList[[1]][, datasetVar == var, drop = FALSE]
            })
            names(rawDataList) <- levels(datasetVar)
        } else {
            stop("Invalid `datasetVar` specification.")
        }

    }
    modal <- .checkArgLen(modal, nlevels(datasetVar))
    lig <- createLiger(rawData = rawDataList, modal = modal, ...)

    cellMetadata <- SummarizedExperiment::colData(object)
    cellMeta(lig, colnames(cellMetadata)) <- cellMetadata
    for (rd in SingleCellExperiment::reducedDimNames(object)) {
        cellMeta(lig, rd) <- SingleCellExperiment::reducedDim(object, rd)
    }

    return(lig)
}

#' @param assay Name of assay to use. Default \code{NULL} uses current active
#' assay.
#' @rdname as.liger
#' @method as.liger Seurat
#' @export
as.liger.Seurat <- function(
        object,
        datasetVar = "orig.ident",
        modal = NULL,
        assay = NULL,
        ...
) {
    mat <- .getSeuratData(object, layer = "counts", slot = "counts",
                          assay = assay)
    datasetVar <- object[[datasetVar]][,1]
    if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
    datasetVar <- droplevels(datasetVar)
    rawDataList <- lapply(levels(datasetVar), function(d) {
        mat[, datasetVar == d, drop = FALSE]
    })
    names(rawDataList) <- levels(datasetVar)
    modal <- .checkArgLen(modal, nlevels(datasetVar))
    createLiger(rawDataList, modal, object[[]], ...)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# From other things to ligerDataset class ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Converting other classes of data to a as.ligerDataset object
#' @export
#' @param object Object.
#' @param modal Modality setting for each dataset.
#' @param ... Additional arguments passed to \code{\link{createLigerDataset}}
#' @return a \linkS4class{liger} object.
#' @rdname as.ligerDataset
#' @examples
#' ctrl <- dataset(pbmc, "ctrl")
#' ctrl
#' as.ligerDataset(ctrl, modal = "atac")
#' rawCounts <- rawData(ctrl)
#' class(rawCounts)
#' as.ligerDataset(rawCounts)
as.ligerDataset <- function(object, ...) UseMethod("as.ligerDataset", object)

#' @rdname as.ligerDataset
#' @export
#' @method as.ligerDataset ligerDataset
as.ligerDataset.ligerDataset <- function(
        object,
        modal = c("default", "rna", "atac", "spatial"),
        ...
) {
    modal <- match.arg(modal)
    newClass <- .modalClassDict[[modal]]
    if (inherits(object, newClass)) return(object)
    slotFromClass <- methods::slotNames(class(object))
    slotToClass <- methods::slotNames(newClass)
    if (any(!slotFromClass %in% slotToClass))
        warning("Will remove information in the following slots when ",
                "converting class from `", class(object), "` to `", newClass,
                "`: ", paste(slotFromClass[!slotFromClass %in% slotToClass],
                             collapse = ", "))
    newCallArgs <- list(Class = newClass)
    for (s in slotFromClass) {
        if (s %in% slotToClass)
            newCallArgs[[s]] <- methods::slot(object, s)
    }
    do.call("new", newCallArgs)
}

#' @rdname as.ligerDataset
#' @export
#' @method as.ligerDataset default
as.ligerDataset.default <- function(
        object,
        modal = c("default", "rna", "atac", "spatial"),
        ...
) {
    modal <- match.arg(modal)
    createLigerDataset(object, modal, ...)
}

#' @rdname as.ligerDataset
#' @export
#' @method as.ligerDataset matrix
as.ligerDataset.matrix <- function(
        object,
        modal = c("default", "rna", "atac", "spatial"),
        ...
) {
    modal <- match.arg(modal)
    object <- methods::as(object, "CsparseMatrix")
    createLigerDataset(object, modal, ...)
}

#' @rdname as.ligerDataset
#' @export
#' @method as.ligerDataset Seurat
#' @param assay Name of assay to use. Default \code{NULL} uses current active
#' assay.
as.ligerDataset.Seurat <- function(
        object,
        modal = c("default", "rna", "atac", "spatial"),
        assay = NULL,
        ...
) {
    modal <- match.arg(modal)
    mat <- .getSeuratData(object, "counts", "counts", assay = assay)
    createLigerDataset(rawData = mat, modal = modal, ...)
}


#' @rdname as.ligerDataset
#' @export
#' @method as.ligerDataset SingleCellExperiment
as.ligerDataset.SingleCellExperiment <- function(
        object,
        modal = c("default", "rna", "atac", "spatial"),
        ...
) {
    if (!requireNamespace("SingleCellExperiment", quietly = "TRUE"))
        stop("Package \"SingleCellExperiment\" needed for this function ",
             "to work. Please install it by command:\n",
             "BiocManager::install('SingleCellExperiment')",
             call. = FALSE)
    modal <- match.arg(modal)
    mat <- SingleCellExperiment::counts(object)
    createLigerDataset(rawData = mat, modal = modal, ...)
}

# #' @rdname as.ligerDataset
# #' @export
# #' @method as.ligerDataset anndata._core.anndata.AnnData
# as.ligerDataset.anndata._core.anndata.AnnData <- function(
#         object,
#         modal = c("default", "rna", "atac", "spatial"),
#         ...
# ) {
#     modal <- match.arg(modal)
#     message("Python object AnnData input. ")
# }

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# From ligerDataset class to other things ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#setAs("ligerDataset", "SingleCellExperiment", function(from) {
#    requireNamespace("SingleCellExperiment")
#    assays <- list()
#    if (!is.null(raw.data(from)))
#        assays <- c(assays, list(counts = raw.data(from)))
#    if (!is.null(normData(from)))
#        assays <- c(assays, list(normcounts = normData(from)))
#    SingleCellExperiment::SingleCellExperiment(
#        assays = assays
#    )
#})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# From liger class to other things ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Convert between liger and Seurat object
#' @description
#' For converting a \linkS4class{liger} object to a Seurat object, the
#' \code{rawData}, \code{normData}, and \code{scaleData} from each dataset,
#' the \code{cellMeta}, \code{H.norm} and \code{varFeatures} slot will be
#' included. Compatible with V4 and V5. It is not recommended to use this
#' conversion if your \linkS4class{liger} object contains datasets from
#' various modalities.
#' @param object A \linkS4class{liger} object to be converted
#' @param assay Name of assay to store the data. Default \code{NULL} detects by
#' dataset modality. If the object contains various modality, default to
#' \code{"liger"}.
#' @param identByDataset Logical, whether to combine dataset variable and default
#' cluster labeling to set the Idents. Default \code{FALSE}.
#' @param by.dataset [Deprecated]. Use \code{identByDataset} instead.
#' @param nms [Defunct] Will be ignored because new object structure does not
#' have related problem.
#' @param renormalize [Defunct] Will be ignored because since Seurat V5, layers
#' of data can exist at the same time and it is better to left it for users to
#' do it by theirselves.
#' @param use.liger.genes [Defunct] Will be ignored and will always set LIGER
#' variable features to the place.
#' @export
ligerToSeurat <- function(
        object,
        assay = NULL,
        identByDataset = FALSE,
        # Rename or defunct
        nms = NULL,
        renormalize = NULL,
        use.liger.genes = NULL,
        by.dataset = identByDataset
) {
    .checkObjVersion(object)
    .deprecateArgs(
        list(by.dataset = "identByDataset"),
        defunct = c("nms", "renormalize", "use.liger.genes"))
    if (is.null(assay)) {
        allModal <- toupper(modalOf(object))
        if (all(allModal == allModal[1])) {
            assay <- allModal[1]
        } else {
            assay <- "liger"
        }
        if (assay == "default") assay <- "RNA"
    }

    rawDataList <- getMatrix(object, "rawData", returnList = TRUE)
    if (any(sapply(rawDataList, is.null))) counts <- NULL
    else counts <- mergeSparseAll(rawDataList)
    normDataList <- getMatrix(object, "normData", returnList = TRUE)
    if (any(sapply(normDataList, is.null))) data <- NULL
    else data <- mergeSparseAll(normDataList)
    scaleDataList <- getMatrix(object, "scaleData", returnList = TRUE)
    if (any(sapply(scaleDataList, is.null))) scale.data <- NULL
    else scale.data <- mergeSparseAll(scaleDataList)

    orig.ident <- object$dataset
    idents <- defaultCluster(object)
    if (is.null(idents)) {
        idents <- orig.ident
    } else {
        if (isTRUE(identByDataset)) {
            idents <- factor(paste0(as.character(orig.ident), "_",
                                    as.character(idents)))
        }
    }

    # Split normal data.frame compatible info and dimReds
    metadata <- data.frame(row.names = colnames(object))
    dimReds <- list()
    for (i in seq_along(cellMeta(object))) {
        varname <- names(cellMeta(object))[i]
        var <- cellMeta(object)[[i]]
        if (is.null(dim(var))) metadata[[varname]] <- var
        else dimReds[[varname]] <- var
    }
    srt <- Seurat::CreateSeuratObject(counts = counts, assay = assay,
                                      meta.data = metadata)

    srt$orig.ident <- orig.ident
    Seurat::Idents(srt) <- idents

    if (!is.null(data)) {
        srt <- .setSeuratData(srt, layer = "ligerNormData", slot = "data",
                              value = data, assay = assay, denseIfNeeded = FALSE)
    }
    if (!is.null(scale.data)) {
        srt <- .setSeuratData(srt, layer = "ligerScaleData", slot = "scale.data",
                              value = scale.data, assay = assay,
                              denseIfNeeded = TRUE)
    }
    # Attempt to get H.norm primarily. If it is NULL, then turn to H
    h <- getMatrix(object, "H.norm") %||%
        getMatrix(pbmc, "H", returnList = TRUE)
    if (is.list(h)) {
        # Not from H.norm but H list. Only get merged when all datasets have H.
        if (any(sapply(h, is.null))) h <- NULL
        else h <- t(Reduce(cbind, h))
    }
    if (!is.null(h)) {
        hDR <- SeuratObject::CreateDimReducObject(
            embeddings = h,
            loadings = getMatrix(object, "W", returnList = FALSE),
            assay = assay,
            misc = list(
                H = getMatrix(object, "H", returnList = TRUE),
                V = getMatrix(object, "V", returnList = TRUE)
            ),
            key = "iNMF_"
        )
        srt[["inmf"]] <- hDR
    }

    for (var in names(dimReds)) {
        dimred <- SeuratObject::CreateDimReducObject(
            embeddings = dimReds[[var]],
            assay = assay,
            key = paste0(var, "_")
        )
        srt[[var]] <- dimred
    }
    Seurat::VariableFeatures(srt) <- varFeatures(object)
    return(srt)
}


seuratToLiger <- as.liger.Seurat
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# From old version to new version ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert old liger object to latest version
#' @param object \code{liger} object from rliger version <1.99.0
#' @param dimredName The name of variable in \code{cellMeta} slot to store the
#' dimensionality reduction matrix, which originally located in
#' \code{tsne.coords} slot. Default \code{"tsne.coords"}.
#' @param clusterName The name of variable in \code{cellMeta} slot to store the
#' clustering assignment, which originally located in \code{clusters} slot.
#' Default \code{"clusters"}.
#' @param h5FilePath Named list, to specify the path to the H5 file of each
#' dataset if location has been changed. Default \code{NULL} looks at the file
#' paths stored in object.
#' @export
#' @examples
#' \dontrun{
#' # Suppose you have a liger object of old version (<1.99.0)
#' newLig <- convertOldLiger(oldLig)
#' }
convertOldLiger <- function(
        object,
        dimredName = "tsne.coords",
        clusterName = "clusters",
        h5FilePath = NULL
) {
    ver120 <- package_version("1.99.0")
    if (object@version >= ver120) return(object)
    if (inherits(object@raw.data[[1]], "H5File")) {
        ldList <- convertOldLiger.H5(object, h5FilePath = h5FilePath)
    } else {
        ldList <- convertOldLiger.mem(object)
    }
    cellMeta <- object@cell.data
    varFeatures <- object@var.genes
    cellID <- unlist(lapply(ldList, colnames), use.names = FALSE)
    # 4. Wrap up liger object
    cellMeta <- S4Vectors::DataFrame(cellMeta)
    # TODO: check default prototype of tsne.coords and clusters.
    dimred <- object@tsne.coords[rownames(cellMeta), , drop = FALSE]
    colnames(dimred) <- seq_len(ncol(dimred))
    cellMeta[[dimredName]] <- dimred
    cellMeta[[clusterName]] <- object@clusters[rownames(cellMeta)]
    rownames(cellMeta) <- cellID
    hnorm <- object@H.norm
    rownames(hnorm) <- cellID
    newObj <- createLiger(ldList, W = t(object@W), H.norm = hnorm,
                          varFeatures = varFeatures, cellMeta = cellMeta,
                          addPrefix = FALSE, removeMissing = FALSE)
    return(newObj)
}

convertOldLiger.mem <- function(object) {
    dataLists <- list()
    if (.hasSlot(object, "raw.data")) dataLists$rawData <- object@raw.data
    if (.hasSlot(object, "norm.data")) dataLists$normData <- object@norm.data
    if (.hasSlot(object, "scale.data")) dataLists$scaleData <- object@scale.data
    if (.hasSlot(object, "H")) dataLists$H <- object@H
    if (.hasSlot(object, "V")) dataLists$V <- object@V
    if (.hasSlot(object, "U")) dataLists$U <- object@U
    # 1. Deal with cell metadata which establish a correct mapping of cell
    # barcodes and datasets belonging
    allDatasets <- Reduce(union, lapply(dataLists, names))
    cellMeta <- object@cell.data
    cellMetaDatasets <- unique(as.vector(cellMeta$dataset))
    if (!identical(sort(allDatasets), sort(cellMetaDatasets))) {
        # Datasets entry for matrices don't match with cell metadata
        # Only take the intersection
        allDatasets <- intersect(allDatasets, cellMetaDatasets)
        cellMeta <- cellMeta[cellMeta[["dataset"]] %in% allDatasets, ]
    }

    # Split `dataLists` by dataset
    datasetLists <- list()
    for (d in allDatasets) {
        for (slot in names(dataLists)) {
            datasetLists[[d]][[slot]] <- dataLists[[slot]][[d]]
        }
    }

    # For each existing dataset
    ldList <- list()
    for (d in allDatasets) {
        # "BC" for barcodes
        # 2. Check and clean up cell barcodes and feature idx issue
        cellMetaBC <- rownames(cellMeta)[cellMeta$dataset == d]
        features <- NULL
        varFeatures <- object@var.genes
        dataList <- datasetLists[[d]]

        # Check cell barcodes
        bcPassing <- .checkIDIdentical(
            ref = cellMetaBC,
            onCol = dataList[c("rawData", "normData")],
            onRow = dataList[c("scaleData", "H")]
        )

        # Check raw, norm data features
        if (!is.null(dataList$rawData)) features <- rownames(dataList$rawData)
        else features <- rownames(dataList$normData)
        if (is.null(features)) {
            warning("Cannot detect feature names for dataset \"", d, "\". ",
                    "Skipped.")
            next
        }
        ftPassing <- .checkIDIdentical(
            ref = features,
            onRow = dataList[c("rawData", "normData")]
        )

        # Check var features
        if (!is.null(dataList$V) &&
            is.null(colnames(dataList$V)) &&
            !is.null(varFeatures)) {
                ## This should not happen but unfortunately, old `V`s might not
                ## have var features as their colnames
                colnames(dataList$V) <- varFeatures
        }
        hvgPassing <- .checkIDIdentical(
            ref = varFeatures,
            onCol = dataList[c("scaleData", "V", "U")]
        )

        # Remove data that has inconsistent information
        passing <- .combinePassingSignal(names(dataList),
                                         bcPassing, ftPassing, hvgPassing)
        dataList <- dataList[passing]
        for (s in c("scaleData", "H", "V", "U")) {
            if (!is.null(dataList[[s]])) {
                dataList[[s]] <- t(dataList[[s]])
            }
        }
        # 3. Construct ligerDataset objects for each dataset
        ldList[[d]] <- do.call(createLigerDataset, dataList)
        colnames(ldList[[d]]) <- paste0(d, "_", colnames(ldList[[d]]))
    }
    return(ldList)
}

convertOldLiger.H5 <- function(object, h5FilePath = NULL) {
    .log("Please use caution when restoring an H5 based liger object, because ",
         "old version does not have solid restriction on cell/feature ",
         "identifier matching. New rliger assumes all data was produced ",
         "with standard old rliger workflow.")
    dataLists <- list()
    if (.hasSlot(object, "H")) dataLists$H <- object@H
    if (.hasSlot(object, "V")) dataLists$V <- object@V
    if (.hasSlot(object, "U")) dataLists$U <- object@U

    # 1. Deal with cell metadata which establish a correct mapping of cell
    # barcodes and datasets belonging
    allDatasets <- Reduce(union, lapply(dataLists, names))
    cellMeta <- object@cell.data
    cellMetaDatasets <- unique(as.vector(cellMeta$dataset))
    if (!identical(sort(allDatasets), sort(cellMetaDatasets))) {
        # Datasets entry for matrices don't match with cell metadata
        # Only take the intersection
        allDatasets <- intersect(allDatasets, cellMetaDatasets)
        cellMeta <- cellMeta[cellMeta[["dataset"]] %in% allDatasets, ]
    }

    # Split `dataLists` by dataset
    datasetLists <- list()
    for (d in allDatasets) {
        for (slot in names(dataLists)) {
            datasetLists[[d]][[slot]] <- dataLists[[slot]][[d]]
        }
    }

    # For each existing dataset
    ldList <- list()
    for (d in allDatasets) {
        # "BC" for barcodes
        # 2. Check and clean up cell barcodes and feature idx issue
        cellMetaBC <- rownames(cellMeta)[cellMeta$dataset == d]
        #features <- NULL
        varFeatures <- object@var.genes
        dataList <- datasetLists[[d]]

        # Check cell barcodes
        bcPassing <- .checkIDIdentical(
            ref = cellMetaBC,
            #onCol = dataList[c("rawData", "normData")],
            onRow = dataList[c("H")]
        )

        # Check raw, norm data features
        # if (!is.null(dataList$rawData)) features <- rownames(dataList$rawData)
        # else features <- rownames(dataList$normData)
        # if (is.null(features)) {
        #     warning("Cannot detect feature names for dataset \"", d, "\". ",
        #             "Skipped.")
        #     next
        # }
        # ftPassing <- .checkIDIdentical(
        #     ref = features,
        #     onRow = dataList[c("rawData", "normData")]
        # )

        # Check var features
        if (!is.null(dataList$V) &&
            is.null(colnames(dataList$V)) &&
            !is.null(varFeatures)) {
            ## This should not happen but unfortunately, old `V`s might not
            ## have var features as their colnames
            colnames(dataList$V) <- varFeatures
        }
        hvgPassing <- .checkIDIdentical(
            ref = varFeatures,
            onCol = dataList[c("V", "U")]
        )

        # Remove data that has inconsistent information
        passing <- .combinePassingSignal(names(dataList),
                                         bcPassing, hvgPassing)
        dataList <- dataList[passing]
        for (s in c("H", "V", "U")) {
            if (!is.null(dataList[[s]])) {
                dataList[[s]] <- t(dataList[[s]])
            }
        }
        # 3. Construct H5 ligerDataset objects for each dataset
        if (!is.null(h5FilePath[[d]])) h5Path <- h5FilePath[[d]]
        else h5Path <- object@h5file.info[[d]]$file.path
        if (!hdf5r::is_hdf5(name = h5Path)) {
            stop("File path for dataset \"", d, "\" not found or is not an H5 ",
                 "file: ", h5Path)
        }
        h5Format <- object@h5file.info[[d]]$format.type
        ldList[[d]] <- do.call(createH5LigerDataset, c(
            list(h5file = h5Path, formatType = h5Format),
            dataList
        ))
        colnames(ldList[[d]]) <- paste0(d, "_", colnames(ldList[[d]]))

        # 4. Check for potential existing processed result
        newSlotNameMap <- list(norm.data = "normData",
                               "scale.data" = "scaleData",
                               "scale.unshared.data" = "scaleUnsharedData")
        for (s in c("norm.data", "scale.data", "scale.unshared.data")) {
            h5file <- getH5File(ldList[[d]])
            if (h5file$link_exists(s)) {
                h5fileInfo(ldList[[d]], newSlotNameMap[[s]], check = FALSE) <- s
            }
        }
    }

    return(ldList)
}

.checkIDIdentical <- function(
    ref,
    onCol = NULL,
    onRow = NULL
) {
    # ref - a character vector as a reference
    # onCol - a list of matrix, where the colnames should match to reference
    # onRow - a list of matrix, where the rownames should match to reference
    colPassing <- rep(TRUE, length(onCol))
    names(colPassing) <- names(onCol)
    for (slot in names(onCol)) {
        if (is.na(slot)) next
        if (!identical(colnames(onCol[[slot]]), ref)) {
            warning("Inconsistent column ID in slot `", slot, "`.")
            colPassing[slot] <- FALSE
        }
    }
    rowPassing <- rep(TRUE, length(onRow))
    names(rowPassing) <- names(onRow)
    for (slot in names(onRow)) {
        if (is.na(slot)) next
        if (!identical(rownames(onRow[[slot]]), ref)) {
            warning("Inconsistent row ID in slot `", slot, "`.")
            rowPassing[slot] <- FALSE
        }
    }
    return(unlist(list(colPassing, rowPassing)))
}

.combinePassingSignal <- function(slotNames, ...) {
    passings <- list(...)
    passings <- lapply(passings, function(x) {
        x <- x[slotNames]
        names(x) <- slotNames
        x[is.na(x)] <- TRUE
        x
    })
    Reduce("&", passings)
}

.hasSlot <- function(object, name) {
    tryCatch(
        expr = {
            methods::slot(object, name)
            return(TRUE)
        },
        error = function(e) {
            .log("Skipped slot `", name, "` which is not available.")
            return(FALSE)
        },
        warining = function(w) {
            .log(w)
            return(FALSE)
        }
    )
}
