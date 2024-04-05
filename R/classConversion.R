#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# From other things to liger class ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname as.liger
#' @export
#' @method as.liger dgCMatrix
as.liger.dgCMatrix <- function(
        object,
        datasetVar = NULL,
        modal = NULL,
        ...
) {
    datasetVar <- datasetVar %||% "sample"
    datasetVar <- .checkArgLen(datasetVar, ncol(object), repN = TRUE,
                               class = c("factor", "character"))
    if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
    datasetVar <- droplevels(datasetVar)

    rawDataList <- splitRmMiss(object, datasetVar)
    modal <- .checkArgLen(modal, length(rawDataList), class = "character")
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
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) # nocov start
        cli::cli_abort(
            "Package {.pkg SingleCellExperiment} is needed for this function to work.
            Please install it by command:
            {.code BiocManager::install('SingleCellExperiment')}"
        )
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
        cli::cli_abort(
            "Package {.pkg SummarizedExperiment} is needed for this function to work.
            Please install it by command:
            {.code BiocManager::install('SummarizedExperiment')}"
        ) # nocov end
    raw <- SingleCellExperiment::counts(object)

    if (is.null(datasetVar)) {
        if ("sample" %in% colnames(SummarizedExperiment::colData(object))) {
            datasetVar <- SummarizedExperiment::colData(object)[["sample"]]
        } else {
            datasetVar <- "SCE"
        }
    } else if (length(datasetVar) == 1) {
        if (datasetVar %in% colnames(SummarizedExperiment::colData(object))) {
            datasetVar <- SummarizedExperiment::colData(object)[[datasetVar]]
        }
    }
    datasetVar <- .checkArgLen(datasetVar, ncol(object), repN = TRUE,
                               class = c("factor", "character"))
    if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
    datasetVar <- droplevels(datasetVar)
    raw <- splitRmMiss(raw, datasetVar)
    modal <- .checkArgLen(modal, length(raw), class = "character")
    lig <- createLiger(raw, modal = modal, ...)
    colDataCopy <- SummarizedExperiment::colData(object)
    for (cdn in colnames(colDataCopy)) {
        if (cdn %in% names(cellMeta(lig))) {
            same <- identical(colDataCopy[[cdn]], cellMeta(lig, cdn))
            if (same) next
            cdnNew <- paste0("SCE_", cdn)
            cli::cli_alert_warning(
                "Variable name {.val {cdn}} in colData of SingleCellExperiment conflicts with liger default variables. Modified to {.val {cdnNew}}.")
        } else {
            cdnNew <- cdn
        }
        cellMeta(lig, cdnNew) <- colDataCopy[[cdn]]
    }
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
        datasetVar = NULL,
        modal = NULL,
        assay = NULL,
        ...
) {
    raw <- .getSeuratData(object, layer = "counts", slot = "counts",
                          assay = assay)
    if (!is.list(raw)) {
        if (is.null(datasetVar)) {
            if ("orig.ident" %in% colnames(object[[]])) {
                datasetVar <- object[["orig.ident", drop = TRUE]]
            } else {
                datasetVar <- "Seurat"
            }
        } else if (length(datasetVar) == 1) {
            if (datasetVar %in% colnames(object[[]])) {
                datasetVar <- object[[datasetVar, drop = TRUE]]
            }
        }
        datasetVar <- .checkArgLen(datasetVar, ncol(object), repN = TRUE,
                                   class = c("factor", "character"))
        if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
        datasetVar <- droplevels(datasetVar)
        raw <- splitRmMiss(raw, datasetVar, rmMiss = FALSE)
    } else {
        names(raw) <- gsub("counts.", "", names(raw))
    }

    datasetVar <- datasetVar %||% "Seurat"
    modal <- .checkArgLen(modal, length(raw), class = "character")
    lig <- createLiger(raw, modal = modal, ...)
    colnames(object) <- colnames(lig)
    for (cdn in colnames(object[[]])) {
        if (cdn %in% names(cellMeta(lig))) {
            same <- identical(object[[cdn, drop = TRUE]], cellMeta(lig, cdn))
            if (same) next
            cdnNew <- paste0("Seurat_", cdn)
            cli::cli_alert_warning(
                "Variable name {.val {cdn}} in meta.data of Seurat conflicts with liger default variables. Modified to {.val {cdnNew}}.")
        } else {
            cdnNew <- cdn
        }
        cellMeta(lig, cdnNew) <- object[[cdn, drop = TRUE]]
    }
    for (rd in SeuratObject::Reductions(object)) {
        mat <- object[[rd]][[]]
        colnames(mat) <- seq_len(ncol(mat))
        dimRed(lig, rd) <- mat
    }
    return(lig)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# From other things to ligerDataset class ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname as.ligerDataset
#' @export
#' @method as.ligerDataset ligerDataset
as.ligerDataset.ligerDataset <- function(
        object,
        modal = c("default", "rna", "atac", "spatial", "meth"),
        ...
) {
    modal <- match.arg(modal)
    newClass <- .modalClassDict[[modal]]
    if (inherits(object, newClass)) return(object)
    slotFromClass <- methods::slotNames(class(object))
    slotToClass <- methods::slotNames(newClass)
    if (any(!slotFromClass %in% slotToClass))
        cli::cli_alert_warning(
            "Will remove information in the following slots when converting class
            from {.cls {class(object)}} to {.cls {newClass}}: {.val {slotFromClass[!slotFromClass %in% slotToClass]}}")
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
        modal = c("default", "rna", "atac", "spatial", "meth"),
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
        modal = c("default", "rna", "atac", "spatial", "meth"),
        ...
) {
    modal <- match.arg(modal)
    createLigerDataset(object, modal, ...)
}

#' @rdname as.ligerDataset
#' @export
#' @method as.ligerDataset Seurat
#' @param assay Name of assay to use. Default \code{NULL} uses current active
#' assay.
as.ligerDataset.Seurat <- function(
        object,
        modal = c("default", "rna", "atac", "spatial", "meth"),
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
        modal = c("default", "rna", "atac", "spatial", "meth"),
        ...
) {
    if (!requireNamespace("SingleCellExperiment", quietly = "TRUE")) # nocov start
        cli::cli_abort(
            "Package {.pkg SingleCellExperiment} is needed for this function to work.
            Please install it by command:
            {.code BiocManager::install('SingleCellExperiment')}"
        ) # nocov end
    modal <- match.arg(modal)
    mat <- SingleCellExperiment::counts(object)
    createLigerDataset(rawData = mat, modal = modal, ...)
}

###### AnnData object presented as H5AD file on disk is already supported with
###### H5-based ligerDataset object without having to create a AnnData object
###### in a running Python session first.

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
#' \code{"LIGER"}. Default dataset modality setting is understood as
#' \code{"RNA"}.
#' @param identByDataset Logical, whether to combine dataset variable and
#' default cluster labeling to set the Idents. Default \code{FALSE}.
#' @param merge Logical, whether to merge layers of different datasets into one.
#' Not recommended. Default \code{FALSE}.
#' @param by.dataset [Deprecated]. Use \code{identByDataset} instead.
#' @param nms [Defunct] Will be ignored because new object structure does not
#' have related problem.
#' @param renormalize [Defunct] Will be ignored because since Seurat V5, layers
#' of data can exist at the same time and it is better to left it for users to
#' do it by themselves.
#' @param use.liger.genes [Defunct] Will be ignored and will always set LIGER
#' variable features to the place.
#' @export
#' @rdname ligerToSeurat
#' @return Always returns Seurat object(s) of the latest version. By default a
#' Seurat object with split layers, e.g. with layers like "counts.ctrl" and
#' "counts.stim". If \code{merge = TRUE}, return a single Seurat object with
#' layers for all datasets merged.
#' @examples
#' seu <- ligerToSeurat(pbmc)
ligerToSeurat <- function(
        object,
        assay = NULL,
        identByDataset = FALSE,
        merge = FALSE,
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
            assay <- "LIGER"
        }
        if (assay == "DEFAULT") assay <- "RNA"
    }

    rawDataList <- getMatrix(object, "rawData", returnList = TRUE)
    rawDataList <- rawDataList[!sapply(rawDataList, is.null)]
    if (isTRUE(merge)) rawDataList <- mergeSparseAll(rawDataList)
    if (!length(rawDataList)) {
        cli::cli_abort("rawData not found.")
    }
    Assay <- SeuratObject::CreateAssay5Object(rawDataList)

    normDataList <- getMatrix(object, "normData", returnList = TRUE)
    normDataList <- normDataList[!sapply(normDataList, is.null)]
    if (isTRUE(merge)) {
        normed <- mergeSparseAll(normDataList)
        SeuratObject::LayerData(Assay, layer = "ligerNormData") <- normed
    } else {
        for (i in seq_along(normDataList)) {
            layerName <- paste0("ligerNormData.", names(normDataList)[i])
            SeuratObject::LayerData(Assay, layer = layerName) <- normDataList[[i]]
        }
    }

    scaleDataList <- getMatrix(object, "scaleData", returnList = TRUE)
    scaleDataList <- scaleDataList[!sapply(scaleDataList, is.null)]
    if (isTRUE(merge)) {
        scaled <- mergeSparseAll(scaleDataList)
        SeuratObject::LayerData(Assay, layer = "ligerScaleData") <- scaled
    } else {
        for (i in seq_along(scaleDataList)) {
            layerName <- paste0("ligerScaleData.", names(scaleDataList)[i])
            SeuratObject::LayerData(Assay, layer = layerName) <- scaleDataList[[i]]
        }
    }

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
    metadata <- .DataFrame.as.data.frame(cellMeta(object))
    dimReds <- dimReds(object)
    srt <- Seurat::CreateSeuratObject(counts = Assay, assay = assay,
                                      meta.data = metadata)

    srt$orig.ident <- orig.ident
    Seurat::Idents(srt) <- idents

    # Attempt to get H.norm primarily. If it is NULL, then turn to H
    h <- getMatrix(object, "H.norm") %||%
        getMatrix(object, "H", returnList = TRUE)
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

#' @rdname as.liger
#' @export
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
convertOldLiger <- function( # nocov start
        object,
        dimredName,
        clusterName = "clusters",
        h5FilePath = NULL
) {
    ver1990 <- package_version("1.99.0")
    if (object@version == ver1990) {
        return(rliger2_to_rliger_namespace(object, dimredName = dimredName))
    }
    if (object@version > ver1990) return(object)
    tryCatch(
        {
            if (inherits(object@raw.data[[1]], "H5File")) {
                ldList <- convertOldLiger.H5(object, h5FilePath = h5FilePath)
            } else {
                ldList <- convertOldLiger.mem(object)
            }
        },
        error = function(e) {
            print(e)
            cli::cli_alert_danger(
                "Conversion failed. Please check the error message above."
            )
            cli::cli_alert_info(
                "For 'inconsistent ID' error, please use an old version of {.pkg rliger} and manually fix the rownames/colnames matching."
            )
            cli::cli_alert("{.code dimnames()} of raw.data and norm.data must be identical for each dataset.")
            cli::cli_alert("{.code rownames()} of scale.data and H must be identical to the colnames of raw.data, for each dataset.")
            cli::cli_alert("{.code colnames()} of scale.data, V and U (if available) must be identical to the var.genes.")
        }
    )

    cellMeta <- object@cell.data
    varFeatures <- object@var.genes
    cellID <- unlist(lapply(ldList, colnames), use.names = FALSE)
    # 4. Wrap up liger object
    cellMeta <- S4Vectors::DataFrame(cellMeta)
    oldID <- rownames(cellMeta)
    # TODO: check default prototype of tsne.coords and clusters.
    dimred <- object@tsne.coords[oldID, , drop = FALSE]
    colnames(dimred) <- paste0(dimredName, "_", seq_len(ncol(dimred)))
    cellMeta$barcode <- oldID
    cellMeta[[clusterName]] <- object@clusters[rownames(cellMeta)]
    rownames(cellMeta) <- cellID
    hnorm <- object@H.norm
    rownames(hnorm) <- cellID
    newObj <- createLiger(ldList, W = t(object@W), H.norm = hnorm,
                          varFeatures = varFeatures, cellMeta = cellMeta,
                          addPrefix = FALSE, removeMissing = FALSE)
    dimRed(newObj, dimredName) <- dimred
    defaultCluster(newObj) <- clusterName
    defaultDimRed(newObj) <- dimredName
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
            cli::cli_abort(
                "Cannot detect feature names for dataset {.val {d}}."
            )
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
        if (!all(startsWith(colnames(ldList[[d]]), d)))
            colnames(ldList[[d]]) <- paste0(d, "_", colnames(ldList[[d]]))
    }
    return(ldList)
}

convertOldLiger.H5 <- function(object, h5FilePath = NULL) {
    cli::cli_alert_warning(
        "Please use caution when restoring an H5 based liger object, because
        the old version does not have any solid restriction on cell/feature
        identifier matching. New rliger > 1.99 assumes all data was produced
        with standard old rliger workflow.", wrap = TRUE
    )
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
            cli::cli_abort("File path for dataset {.val {d}} not found or is not an H5 file: {.file {h5Path}}")
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
            cli::cli_alert_danger("Inconsistent column ID in slot {slot}.")
            colPassing[slot] <- FALSE
        }
    }
    rowPassing <- rep(TRUE, length(onRow))
    names(rowPassing) <- names(onRow)
    for (slot in names(onRow)) {
        if (is.na(slot)) next
        if (!identical(rownames(onRow[[slot]]), ref)) {
            cli::cli_alert_danger("Inconsistent row ID in slot {slot}.")
            rowPassing[slot] <- FALSE
        }
    }
    return(c(colPassing, rowPassing))
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
            cli::cli_alert_info("Skipped slot {name} which is not available.")
            return(FALSE)
        },
        warining = function(w) {
            cli::cli_alert_warning(w)
            return(FALSE)
        }
    )
} # nocov end
