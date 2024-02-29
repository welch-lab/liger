setClass("ligerDataset")
setClassUnion("matrixLike", c("matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# From other things to liger class ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Converting other classes of data to a liger object
#' @description
#' This function converts data stored in SingleCellExperiment (SCE), Seurat
#' object or a merged sparse matrix (dgCMatrix) into a liger object. This is
#' designed for a container object or matrix that already contains multiple
#' datasets to be integerated with LIGER. For individual datasets, please use
#' \code{\link{createLiger}} instead.
#' @export
#' @param object Object.
#' @param datasetVar Specify the dataset belonging by: 1. Select a variable from
#' existing metadata in the object (e.g. colData column); 2. Specify a
#' vector/factor that assign the dataset belonging. 3. Give a single character
#' string which means that all data is from one dataset (must not be a metadata
#' variable, otherwise it is understood as 1.). Default \code{NULL} gathers
#' things into one dataset and names it "sample" for dgCMatrix, attempts
#' to find variable "sample" from SCE or "orig.ident" from Seurat.
#' @param modal Modality setting for each dataset. See
#' \code{\link{createLiger}}.
#' @param ... Additional arguments passed to \code{\link{createLiger}}
#' @details
#' For Seurat V5 structure, it is highly recommended that users make use of its
#' split layer feature, where things like "counts", "data", and "scale.data"
#' can be held for each dataset in the same Seurat object, e.g. with
#' "count.ctrl", "count.stim", not merged. If a Seurat object with split layers
#' is given, \code{datasetVar} will be ignored and the layers will be directly
#' used.
#' @return a \linkS4class{liger} object.
#' @rdname as.liger
#' @examples
#' # dgCMatrix (common sparse matrix class), usually obtained from other
#' # container object, and contains multiple samples merged in one.
#' matList <- rawData(pbmc)
#' multiSampleMatrix <- mergeSparseAll(matList)
#' # The `datasetVar` argument expects the variable assigning the sample source
#' pbmc2 <- as.liger(multiSampleMatrix, datasetVar = pbmc$dataset)
#' pbmc2
#'
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = multiSampleMatrix)
#' )
#' sce$sample <- pbmc$dataset
#' pbmc3 <- as.liger(sce, datasetVar = "sample")
#' pbmc3
#'
#' seu <- SeuratObject::CreateSeuratObject(multiSampleMatrix)
#' # Seurat creates variable "orig.ident" by identifying the cell barcode
#' # prefixes, which is indeed what we need in this case. Users might need
#' # to be careful and have it confirmed first.
#' pbmc4 <- as.liger(seu, datasetVar = "orig.ident")
#' pbmc4
#'
#' # As per Seurat V5 updates with layered data, specifically helpful udner the
#' # scenario of dataset integration. "counts" and etc for each datasets can be
#' # split into layers.
#' seu5 <- seu
#' seu5[["RNA"]] <- split(seu5[["RNA"]], pbmc$dataset)
#' print(SeuratObject::Layers(seu5))
#' pbmc5 <- as.liger(seu5)
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
    datasetVar <- datasetVar %||% "sample"
    datasetVar <- .checkArgLen(datasetVar, ncol(object), repN = TRUE,
                               class = c("factor", "character"))
    if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
    datasetVar <- droplevels(datasetVar)

    rawDataList <- splitRmMiss(object, datasetVar)
    modal <- .checkArgLen(modal, length(rawDataList))
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
        stop("Package \"SingleCellExperiment\" needed for this function ",
             "to work. Please install it by command:\n",
             "BiocManager::install('SingleCellExperiment')",
             call. = FALSE)
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
        stop("Package \"SummarizedExperiment\" needed for this function ",
             "to work. Please install it by command:\n",
             "BiocManager::install('SummarizedExperiment')",
             call. = FALSE) # nocov end
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
    modal <- .checkArgLen(modal, length(raw))
    lig <- createLiger(raw, modal = modal, ...)
    colDataCopy <- SummarizedExperiment::colData(object)
    for (cdn in colnames(colDataCopy)) {
        if (cdn %in% names(cellMeta(lig))) {
            same <- identical(colDataCopy[[cdn]], cellMeta(lig, cdn))
            if (same) next
            cdnNew <- paste0("SCE_", cdn)
            warning("Variable name \"", cdn, "\" in colData of SingleCellExperiment ",
                    "conflicts with liger default variables. Modified to ", cdnNew, ".")
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
        raw <- splitRmMiss(raw, datasetVar)
    } else {
        names(raw) <- gsub("counts.", "", names(raw))
    }

    datasetVar <- datasetVar %||% "Seurat"
    modal <- .checkArgLen(modal, length(raw))
    lig <- createLiger(raw, modal = modal, ...)
    colnames(object) <- colnames(lig)
    for (cdn in colnames(object[[]])) {
        if (cdn %in% names(cellMeta(lig))) {
            same <- identical(object[[cdn, drop = TRUE]], cellMeta(lig, cdn))
            if (same) next
            cdnNew <- paste0("Seurat_", cdn)
            warning("Variable name \"", cdn, "\" in meta.data of Seurat ",
                    "conflicts with liger default variables. Modified to ", cdnNew, ".")
        } else {
            cdnNew <- cdn
        }
        cellMeta(lig, cdnNew) <- object[[cdn, drop = TRUE]]
    }
    for (rd in SeuratObject::Reductions(object)) {
        mat <- object[[rd]][[]]
        colnames(mat) <- seq_len(ncol(mat))
        cellMeta(lig, rd) <- mat
    }
    return(lig)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# From other things to ligerDataset class ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Converting other classes of data to a as.ligerDataset object
#' @description
#' Works for converting a matrix or container object to a single ligerDataset,
#' and can also convert the modality preset of a ligerDataset. When used with
#' a dense matrix object, it automatically converts the matrix to sparse form
#' (\code{\link[Matrix]{dgCMatrix-class}}). When used with container objects
#' such as Seurat or SingleCellExperiment, it is highly recommended that the
#' object contains only one dataset/sample which is going to be integrated with
#' LIGER. For multi-sample objects, please use \code{\link{as.liger}} with
#' dataset source variable specified.
#' @export
#' @param object Object.
#' @param modal Modality setting for each dataset. Choose from \code{"default"},
#' \code{"rna"}, \code{"atac"}, \code{"spatial"}, \code{"meth"}.
#' @param ... Additional arguments passed to \code{\link{createLigerDataset}}
#' @return a \linkS4class{liger} object.
#' @rdname as.ligerDataset
#' @examples
#' ctrl <- dataset(pbmc, "ctrl")
#' ctrl
#' # Convert the modality preset
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
        modal = c("default", "rna", "atac", "spatial", "meth"),
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
        stop("Package \"SingleCellExperiment\" needed for this function ",
             "to work. Please install it by command:\n",
             "BiocManager::install('SingleCellExperiment')",
             call. = FALSE) # nocov end
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
        stop("rawData not found.")
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
    metadata <- data.frame(row.names = colnames(object))
    dimReds <- list()
    for (i in seq_along(cellMeta(object))) {
        varname <- names(cellMeta(object))[i]
        var <- cellMeta(object)[[i]]
        if (is.null(dim(var))) metadata[[varname]] <- var
        else dimReds[[varname]] <- var
    }
    srt <- Seurat::CreateSeuratObject(counts = Assay, assay = assay,
                                      meta.data = metadata)

    srt$orig.ident <- orig.ident
    Seurat::Idents(srt) <- idents

    # if (!is.null(data)) {
    #     srt <- .setSeuratData(srt, layer = "ligerNormData", slot = "data",
    #                           value = data, assay = assay, denseIfNeeded = FALSE)
    # }
    # if (!is.null(scale.data)) {
    #     srt <- .setSeuratData(srt, layer = "ligerScaleData", slot = "scale.data",
    #                           value = scale.data, assay = assay,
    #                           denseIfNeeded = TRUE)
    # }
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
} # nocov end
