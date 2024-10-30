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
        raw <- splitRmMiss(raw, datasetVar)
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
#' @param by.dataset `r lifecycle::badge("superseded")`. Use
#' \code{identByDataset} instead.
#' @param nms `r lifecycle::badge("defunct")` Will be ignored because new object
#' structure does not have related problem.
#' @param renormalize `r lifecycle::badge("defunct")` Will be ignored because
#' since Seurat V5, layers of data can exist at the same time and it is better
#' to left it for users to do it by themselves.
#' @param use.liger.genes `r lifecycle::badge("defunct")` Will be ignored and
#' will always set LIGER variable features to the place.
#' @export
#' @rdname ligerToSeurat
#' @return Always returns Seurat object(s) of the latest version. By default a
#' Seurat object with split layers, e.g. with layers like "counts.ctrl" and
#' "counts.stim". If \code{merge = TRUE}, return a single Seurat object with
#' layers for all datasets merged.
#' @examples
#' if (requireNamespace("SeuratObject", quietly = TRUE) &&
#'     requireNamespace("Seurat", quietly = TRUE)) {
#'     seu <- ligerToSeurat(pbmc)
#' }
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

#' Update old liger object to up-to-date structure
#' @description
#' Due to massive updates since rliger 2.0, old liger object structures are no
#' longer compatible with the current package. This function will update the
#' object to the latest structure.
#' @param object An object of any version of rliger
#' @param dimredName Name of the dimension reduction embedding to be stored.
#' Please see Details section.
#' @param clusterName Name of the clustering assignment variable to be stored.
#' Please see Details section.
#' @param h5FilePath Named character vector for all H5 file paths. Not required
#' for object run with in-memory analysis. For object containing H5-based
#' analysis (e.g. online iNMF), this must be supplied if the H5 file location is
#' different from that at creation time.
#' @details
#' Old liger object (<1.99.0) stores only one embedding at slot
#' \code{tsne.coords}. \code{dimredName} must be specified as a single
#' character. Pre-release version (1.99.0) stores multiple embeddings in
#' \code{cellMeta}. \code{dimredName} must be exact existing variable names in
#' \code{cellMeta} slot.
#'
#' Old liger object stores clustering assignment in slot \code{clusters}.
#' \code{clusterName} must be specified as a single character. Pre-release
#' version does not require this.
#' @return Updated liger object.
#' @export
#' @examples
#' \dontrun{
#' # Suppose you have a liger object of old version (<1.99.0)
#' newLig <- updateLigerObject(oldLig,
#'                             dimredName = "UMAP",
#'                             clusterName = "louvain")
#' }
updateLigerObject <- function(
        object,
        dimredName,
        clusterName = "clusters",
        h5FilePath = NULL
) {
    objClass <- class(object)
    if (objClass[1] != "liger") {
        # DUE TO THE FACT THAT RLIGER PACKAGE HAS BEEN PREVIOUSLY NAMED BY
        # 'liger' AND 'rliger2', `inherits()` WON'T WORK FOR OBJECTS CREATED
        # WITH OLD NAMESPACES WHICH IS NOT INSTALLED.
        cli::cli_abort("object is not of class {.cls liger}.")
    }
    attr(class(object), "package") <- "rliger"
    # Instead of looking at the "version" number encoded in the object,
    # More safe to check for slots
    if (methods::.hasSlot(object, "raw.data") &&
        methods::.hasSlot(object, "norm.data") &&
        methods::.hasSlot(object, "scale.data")) {
        # Old structure
        cli::cli_alert_info("Detected {.cls liger} object with old-style structure.")
        object <- updateOldLiger(object, dimredName = dimredName,
                                 clusterName = clusterName,
                                 h5FilePath = h5FilePath)
    } else if (methods::.hasSlot(object, "datasets") &&
               methods::.hasSlot(object, "cellMeta")) {
        if (methods::.hasSlot(object, "dimReds")) {
            # Current structure
            cli::cli_alert_info("Detected {.cls liger} object with up-to-date structure.")
            if (isH5Liger(object)) {
                object <- restoreH5Liger(object, filePath = h5FilePath)
            }
            methods::slot(object, "version") <- package_version(version)
        } else {
            # 1.99.0 dev version structure
            cli::cli_alert_info("Detected {.cls liger} object with pre-release structure.")
            object <- updateRliger2NS(object, dimredName = dimredName)
        }
    }

    return(object)
}

updateOldLiger <- function(
        object,
        dimredName,
        clusterName = "clusters",
        h5FilePath = NULL
) {
    # Instead of having strict check points everywhere, this updating function
    # focuses on constructing new structure with ligerDataset class, and then
    # fill in as much information as possible. Old rliger package didn't have
    # many checks so the given object might be somehow broken. We don't want to
    # really refuse to update when users might really need it.

    # First detect the set of datasets from all expression data slots and
    # metadata slot, so that we don't miss any dataset partially removed for
    # some reason.
    if (missing(dimredName)) {
        cli::cli_alert_warning("{.field dimredName} not specified. Setting default {.val embedding}.")
        cli::cli_alert_info("It can be renamed with command: {.code names(dimReds(obj)) <- 'newName'}")
        dimredName <- "embedding"
    }
    datasetNames <- Reduce(union, list(
        names(methods::slot(object, "raw.data")),
        names(methods::slot(object, "norm.data"))
    ))

    if (!is.null(h5FilePath)) {
        if (is.null(names(h5FilePath)) ||
            !is.character(h5FilePath)) {
            cli::cli_abort("{.field h5FilePath} must be a named character vector.")
        }
    }
    ldList <- list()
    ldmiscList <- list()
    for (i in seq_along(datasetNames)) {
        dn <- datasetNames[i]
        cli::cli_alert_info("Constructing dataset {.val {dn}}")
        components <- list(Class = "ligerDataset")
        # Old liger structure does not have modality presets, so
        # Class = "ligerDataset"
        h5meta <- NULL
        if (methods::.hasSlot(object, "h5file.info")) {
            h5meta <- methods::slot(object, "h5file.info")[[dn]]
        }

        if (!is.null(h5meta)) {
            # H5 mode
            dsh5path <- NULL
            if (!dn %in% names(h5FilePath)) {
                dsh5path <- h5meta$file.path
            } else {
                dsh5path <- h5FilePath[[dn]]
            }
            if (!file.exists(dsh5path)) {
                cli::cli_alert_danger("H5 file not found for {.val {dn}}: {.file {dsh5path}}. Skipped.")
                next
            }
            components <- tryCatch({
                h5file <- hdf5r::H5File$new(dsh5path, mode = "r+")
                formatType <- h5meta$format.type
                if (formatType == "10X") {
                    barcodesPath <- "matrix/barcodes"
                    genePath <- "matrix/features/name"
                    xPath <- "matrix/data"
                    iPath <- "matrix/indices"
                    pPath <- "matrix/indptr"
                } else if (formatType == "AnnData") {
                    barcodesPath <- sprintf("obs/%s", hdf5r::h5attr(h5file[['obs']], "_index"))
                    genePath <- sprintf("raw/var/%s", hdf5r::h5attr(h5file[['raw/var']], "_index"))
                    xPath <- "raw/X/data"
                    iPath <- "raw/X/indices"
                    pPath <- "raw/X/indptr"
                } else {
                    cli::cli_abort("Unsupported H5 format type {.val {formatType}}.")
                }
                rn <- h5file[[genePath]][]
                cn <- h5file[[barcodesPath]][]
                normPath <- if ("norm.data" %in% names(h5file)) "norm.data" else NULL
                if ("scale.data" %in% names(h5file)) {
                    cli::cli_alert_warning("Old H5 based scaled data format is no longer supported.")
                    cli::cli_alert_info("Please rerun {.code scaleNotCenter()} after loading.")
                }
                h5metaNew <- list(
                    H5File = h5file,
                    filename = h5file$filename,
                    formatType = formatType,
                    indicesName = iPath,
                    indptrName = pPath,
                    barcodesName = barcodesPath,
                    genesName = genePath,
                    rawData = xPath,
                    normData = normPath,
                    scaleData = NULL
                )
                components$rawData <- h5file[[xPath]]
                if (!is.null(normPath)) components$normData <- h5file[[normPath]]
                components$h5fileInfo <- h5metaNew
                components$colnames <- cn
                components$rownames <- rn
                components
            }, error = function(e) {
                cli::cli_alert_danger("Failed to extract necessary H5 information for {.val {dn}}. Skipped.")
                cli::cli_alert_danger(e$message)
                return(list())
            })
            if (length(components) == 0)
                next
        } else {
            # in-memory mode
            ldmiscList[[dn]] <- list()
            components$rawData <- methods::slot(object, "raw.data")[[dn]]
            components$normData <- methods::slot(object, "norm.data")[[dn]]
            if (!is.null(components$rawData)) {
                cn <- colnames(components$rawData)
                rn <- rownames(components$rawData)
            } else if (!is.null(components$normData)) {
                cn <- colnames(components$normData)
                rn <- rownames(components$normData)
            } else {
                cli::cli_alert_danger("Failed because raw and normalized matrices are both missing for cell and feature name detection. Skipped.")
                next
            }
            components$colnames <- cn
            components$rownames <- rn

            # Transpose and sparsify scale data
            ldScaleData <- methods::slot(object, "scale.data")[[dn]]
            if (!.skip(ldScaleData, "scaled data")) {
                ldScaleData <- t(methods::as(ldScaleData, "CsparseMatrix"))
                if (!all(rownames(ldScaleData) %in% rn)) {
                    cli::cli_alert_warning("Scaled data contains features that do not exist in raw data. Including into unstructured misc data.")
                    ldmiscList[[dn]]$scaleData <- ldScaleData
                }
            }

            if (methods::.hasSlot(object, "scale.unshared.data")) {
                ldScaleUnsharedData <- methods::slot(object, "scale.unshared.data")[[dn]]
                if (!.skip(ldScaleUnsharedData, "scaled data of unshared features")) {
                    ldScaleUnsharedData <- t(methods::as(ldScaleUnsharedData, "CsparseMatrix"))
                    components$scaleUnsharedData <- ldScaleUnsharedData
                }
            }
        }

        if (methods::.hasSlot(object, "var.unshared.features")) {
            vuf <- methods::slot(object, "var.unshared.features")[[dn]]
            if (!.skip(vuf, "unshared features")) {
                components$varUnsharedFeatures <- vuf
            }
        }

        ldH <- methods::slot(object, "H")[[dn]]
        if (!.skip(ldH, "factorized H matrix")) {
            components$H <- t(ldH)
        }

        ldV <- methods::slot(object, "V")[[dn]]
        if (!.skip(ldV, "factorized V matrix")) {
            ldV <- t(ldV)
            if (is.null(rownames(ldV)) ||
                !all(rownames(ldV) %in% rn)) {
                cli::cli_alert_warning("Invalid feature names in factorized V matrix. Including into unstructured misc data.")
                ldmiscList[[dn]]$V <- ldV
            } else {
                components$V <- ldV
            }
        }

        if (methods::.hasSlot(object, "A")) {
            ldA <- methods::slot(object, "A")[[dn]]
            if (!.skip(ldA, "onlineINMF intermediate A matrix")) {
                components$A <- t(ldA)
            }
        }
        if (methods::.hasSlot(object, "B")) {
            ldB <- methods::slot(object, "B")[[dn]]
            if (!.skip(ldB, "onlineINMF intermediate B matrix")) {
                components$B <- t(ldB)
            }
        }
        if (methods::.hasSlot(object, "U")) {
            ldU <- methods::slot(object, "U")[[dn]]
            if (!.skip(ldU, "UINMF factorized U matrix")) {
                components$U <- t(ldU)
            }
        }
        components$featureMeta <- S4Vectors::DataFrame(row.names = rn)

        ldList[[dn]] <- do.call("new", components)
        cli::cli_alert_success("Constructed dataset {.val {dn}}")
    }
    ligerComp <- list(Class = "liger",
                      datasets = ldList,
                      uns = list(defaultCluster = clusterName,
                                 defaultDimRed = dimredName,
                                 datasetMist = ldmiscList),
                      version = utils::packageVersion("rliger"))
    # All other all-cell individual elements
    vg <- methods::slot(object, "var.genes")
    if (!.skip(vg, "variable genes")) {
        allisec <- Reduce(intersect, lapply(ldList, rownames))
        if (any(!vg %in% allisec)) {
            diff <- setdiff(vg, allisec)
            cli::cli_alert_warning("{.val {length(diff)}} variable gene{?s} not found in all datasets ({.val {diff}}). Including all into unstructured misc data.")
            ligerComp$uns$varFeatures <- diff
        } else {
            ligerComp$varFeatures <- vg
        }
    }

    hnorm <- methods::slot(object, "H.norm")
    if (!.skip(hnorm, "aligned H matrix")) {
        ligerComp$H.norm <- hnorm
    }

    w <- methods::slot(object, "W")
    if (!.skip(w, "factorized W matrix")) {
        ligerComp$W <- t(w)
    }

    aggData <- methods::slot(object, "agg.data")
    if (!.skip(aggData, "data aggregated within clusters")) {
        cli::cli_alert_warning("Data aggregated within clusters is included into unstructured misc data.")
        ligerComp$uns$agg.data <- aggData
    }

    param <- methods::slot(object, "parameters")
    if (!.skip(param, "parameters")) {
        cli::cli_alert_warning("Parameters are included into unstructured misc data.")
        ligerComp$uns$parameters <- param
    }

    snf <- methods::slot(object, "snf")
    if (!.skip(snf, "SNF intermediate information")) {
        cli::cli_alert_warning("SNF info is included into unstructured misc data.")
        ligerComp$uns$snf <- snf
    }

    dimredList <- list()
    tsnecoords <- methods::slot(object, "tsne.coords")
    if (!.skip(tsnecoords, "dimension reduction")) {
        dimredList[[dimredName]] <- tsnecoords
        ligerComp$dimReds <- dimredList
    }

    # Finally for metadata
    metadata <- methods::slot(object, "cell.data")
    metadata$dataset <- NULL

    alignCluster <- methods::slot(object, "alignment.clusters")
    if (!.skip(alignCluster, "alignment clusters")) {
        metadata <- tryCatch({
            metadata[names(alignCluster), "alignCluster"] <- alignCluster
            metadata
        }, error = function(e) {
            cli::cli_alert_danger("Failed to add alignment clusters to cellMeta. Included into unstructured misc data.")
            cli::cli_alert_danger(e$message)
            ligerComp$uns$alignment_clusters <- alignCluster
            return(metadata)
        })
    }

    clusters <- methods::slot(object, "clusters")
    if (!.skip(clusters, "cluster assignment")) {
        metadata <- tryCatch({
            metadata[names(clusters), clusterName] <- clusters
            metadata
        }, error = function(e) {
            cli::cli_alert_danger("Failed to add cluster assignment to cellMeta. Included into unstructured misc data.")
            cli::cli_alert_danger(e$message)
            ligerComp$uns[[clusterName]] <- clusters
            return(metadata)
        })
    }

    # Make some effort to fix this part
    bcOrignal <- rownames(metadata)
    bcNew <- unlist(lapply(ldList, colnames), use.names = FALSE)
    if (any(!bcNew %in% bcOrignal)) {
        diff <- setdiff(bcNew, bcOrignal)
        cli::cli_alert_warning("{.val {length(diff)}} cell{?s} from expression matrices not found in metadata. Filling in NA{?s}.")
        fill <- data.frame(row.names = diff)
        for (i in seq(ncol(metadata))) {
            fill[[colnames(metadata)[i]]] <- NA
        }
        metadata <- rbind(metadata, fill)
    }
    if (any(!bcOrignal %in% bcNew)) {
        diff <- setdiff(bcOrignal, bcNew)
        cli::cli_alert_warning("{.val {length(diff)}} cell{?s} from metadata not found in expression matrices. Cutting out into unstructured misc data.")
        ligerComp$uns$cell.data.missing <- metadata[diff, ]
    }
    metadata <- metadata[bcNew, ]
    metadata <- cbind(
        dataset = factor(
            rep(names(ldList), sapply(ldList, ncol)),
            levels = names(ldList)
        ),
        metadata
    )
    ligerComp$cellMeta <- S4Vectors::DataFrame(metadata)

    new <- tryCatch({
        do.call("new", ligerComp)
    },
    error = function(e) {
        cli::cli_alert_danger("Failed to create new liger object basing on existing information")
        cli::cli_alert_danger(e$message)
        cli::cli_alert_info("Returning list structure of all extracted information.")
        return(ligerComp)
    })

    return(new)
}

.skip <- function(component, name) {
    if (length(component) == 0) {
        cli::cli_alert_warning("No {name} found. Skipped.")
        return(TRUE)
    } else {
        return(FALSE)
    }
}

updateRliger2NS <- function(
        object,
        dimredName,
        clusterName = NULL,
        h5FilePath = NULL
) {
    # Slots that contain things of rliger2 namespace
    datasetList <- list()
    drList <- list()
    commandList <- list()
    cm <- methods::slot(object, "cellMeta")
    if (missing(dimredName)) {
        cli::cli_alert_warning("{.field dimredName} not specified. Keeping them all in {.field cellMeta} if they exist.")
        multiDimVarIdx <- sapply(cm, function(v) !is.null(dim(v)))
        guessArg <- paste0('"', names(cm)[multiDimVarIdx], '"', collapse = ", ")
        cli::cli_alert_info("Guess you need: {.code dimredName = c({guessArg})}")
        dimredName <- NULL
    }
    for (i in dimredName) {
        if (i %in% colnames(cm)) {
            drList[[i]] <- cm[[i]]
            cm[[i]] <- NULL
        } else {
            cli::cli_alert_danger("Requested dimension reduction {.val {i}} does not exist. Skipped")
        }
    }

    for (i in names(methods::slot(object, "datasets"))) {
        ld <- methods::slot(object, "datasets")[[i]]
        basics <- list(
            rawData = methods::slot(ld, "rawData"),
            normData = methods::slot(ld, "normData"),
            scaleData = methods::slot(ld, "scaleData"),
            H = methods::slot(ld, "H"),
            V = methods::slot(ld, "V"),
            A = methods::slot(ld, "A"),
            B = methods::slot(ld, "B"),
            varUnsharedFeatures = methods::slot(ld, "varUnsharedFeatures"),
            scaleUnsharedData = methods::slot(ld, "scaleUnsharedData"),
            U = methods::slot(ld, "U"),
            h5fileInfo = methods::slot(ld, "h5fileInfo"),
            featureMeta = methods::slot(ld, "featureMeta"),
            colnames = methods::slot(ld, "colnames"),
            rownames = methods::slot(ld, "rownames")
        )
        if (!is.null(basics$scaleData)) {
            basics$scaleData <- methods::as(basics$scaleData, "CsparseMatrix")
        }
        if (!is.null(basics$scaleUnsharedData)) {
            basics$scaleUnsharedData <- methods::as(basics$scaleUnsharedData,
                                                    "CsparseMatrix")
        }
        if (class(ld)[1] == "ligerATACDataset") {
            basics <- c(basics, list(
                Class = "ligerATACDataset",
                rawPeak = methods::slot(ld, "rawPeak"),
                normPeak = methods::slot(ld, "normPeak")
            ))
        } else if (class(ld)[1] == "ligerRNADataset") {
            basics <- c(basics, list(
                Class = "ligerRNADataset"
            ))
        } else if (class(ld)[1] == "ligerMethDataset") {
            basics <- c(basics, list(
                Class = "ligerMethDataset"
            ))
        } else if (class(ld)[1] == "ligerSpatialDataset") {
            basics <- c(basics, list(
                Class = "ligerSpatialDataset",
                coordinate = methods::slot(ld, "coordinate")
            ))
        } else {
            basics <- c(basics, list(
                Class = "ligerDataset"
            ))
        }
        datasetList[[i]] <- do.call("new", basics)
    }

    for (i in names(methods::slot(object, "commands"))) {
        cmd <- methods::slot(object, "commands")[[i]]
        commandList[[i]] <- methods::new(
            "ligerCommand",
            funcName = methods::slot(cmd, "funcName"),
            time = methods::slot(cmd, "time"),
            call = methods::slot(cmd, "call"),
            parameters = methods::slot(cmd, "parameters"),
            objSummary = methods::slot(cmd, "objSummary"),
            ligerVersion = methods::slot(cmd, "ligerVersion"),
            dependencyVersion = methods::slot(cmd, "dependencyVersion")
        )
    }

    new <- methods::new(
        "liger",
        datasets = datasetList,
        cellMeta = cm,
        varFeatures = methods::slot(object, "varFeatures"),
        W = methods::slot(object, "W"),
        H.norm = methods::slot(object, "H.norm"),
        uns = methods::slot(object, "uns"),
        commands = commandList,
        version = utils::packageVersion("rliger")
    )

    for (i in names(drList)) {
        dimRed(new, i) <- drList[[i]]
    }

    new
}


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
convertOldLiger <- updateLigerObject

# convertOldLiger <- function( # nocov start
#         object,
#         dimredName,
#         clusterName = "clusters",
#         h5FilePath = NULL
# ) {
#     ver1990 <- package_version("1.99.0")
#     if (object@version == ver1990) {
#         return(rliger2_to_rliger_namespace(object, dimredName = dimredName))
#     }
#     if (object@version > ver1990) return(object)
#     tryCatch(
#         {
#             if (inherits(object@raw.data[[1]], "H5File")) {
#                 ldList <- convertOldLiger.H5(object, h5FilePath = h5FilePath)
#             } else {
#                 ldList <- convertOldLiger.mem(object)
#             }
#         },
#         error = function(e) {
#             print(e)
#             cli::cli_alert_danger(
#                 "Conversion failed. Please check the error message above."
#             )
#             cli::cli_alert_info(
#                 "For 'inconsistent ID' error, please use an old version of {.pkg rliger} and manually fix the rownames/colnames matching."
#             )
#             cli::cli_alert("{.code dimnames()} of raw.data and norm.data must be identical for each dataset.")
#             cli::cli_alert("{.code rownames()} of scale.data and H must be identical to the colnames of raw.data, for each dataset.")
#             cli::cli_alert("{.code colnames()} of scale.data, V and U (if available) must be identical to the var.genes.")
#         }
#     )
#
#     cellMeta <- object@cell.data
#     varFeatures <- object@var.genes
#     cellID <- unlist(lapply(ldList, colnames), use.names = FALSE)
#     # 4. Wrap up liger object
#     cellMeta <- S4Vectors::DataFrame(cellMeta)
#     oldID <- rownames(cellMeta)
#     # TODO: check default prototype of tsne.coords and clusters.
#     dimred <- object@tsne.coords[oldID, , drop = FALSE]
#     colnames(dimred) <- paste0(dimredName, "_", seq_len(ncol(dimred)))
#     cellMeta$barcode <- oldID
#     cellMeta[[clusterName]] <- object@clusters[rownames(cellMeta)]
#     rownames(cellMeta) <- cellID
#     hnorm <- object@H.norm
#     rownames(hnorm) <- cellID
#     newObj <- createLiger(ldList, W = t(object@W), H.norm = hnorm,
#                           varFeatures = varFeatures, cellMeta = cellMeta,
#                           addPrefix = FALSE, removeMissing = FALSE)
#     dimRed(newObj, dimredName) <- dimred
#     defaultCluster(newObj) <- clusterName
#     defaultDimRed(newObj) <- dimredName
#     return(newObj)
# }
#
# convertOldLiger.mem <- function(object) {
#     dataLists <- list()
#     if (.hasSlot(object, "raw.data")) dataLists$rawData <- object@raw.data
#     if (.hasSlot(object, "norm.data")) dataLists$normData <- object@norm.data
#     if (.hasSlot(object, "scale.data")) dataLists$scaleData <- object@scale.data
#     if (.hasSlot(object, "H")) dataLists$H <- object@H
#     if (.hasSlot(object, "V")) dataLists$V <- object@V
#     if (.hasSlot(object, "U")) dataLists$U <- object@U
#     # 1. Deal with cell metadata which establish a correct mapping of cell
#     # barcodes and datasets belonging
#     allDatasets <- Reduce(union, lapply(dataLists, names))
#     cellMeta <- object@cell.data
#     cellMetaDatasets <- unique(as.vector(cellMeta$dataset))
#     if (!identical(sort(allDatasets), sort(cellMetaDatasets))) {
#         # Datasets entry for matrices don't match with cell metadata
#         # Only take the intersection
#         allDatasets <- intersect(allDatasets, cellMetaDatasets)
#         cellMeta <- cellMeta[cellMeta[["dataset"]] %in% allDatasets, ]
#     }
#
#     # Split `dataLists` by dataset
#     datasetLists <- list()
#     for (d in allDatasets) {
#         for (slot in names(dataLists)) {
#             datasetLists[[d]][[slot]] <- dataLists[[slot]][[d]]
#         }
#     }
#
#     # For each existing dataset
#     ldList <- list()
#     for (d in allDatasets) {
#         # "BC" for barcodes
#         # 2. Check and clean up cell barcodes and feature idx issue
#         cellMetaBC <- rownames(cellMeta)[cellMeta$dataset == d]
#         features <- NULL
#         varFeatures <- object@var.genes
#         dataList <- datasetLists[[d]]
#
#         # Check cell barcodes
#         bcPassing <- .checkIDIdentical(
#             ref = cellMetaBC,
#             onCol = dataList[c("rawData", "normData")],
#             onRow = dataList[c("scaleData", "H")]
#         )
#
#         # Check raw, norm data features
#         if (!is.null(dataList$rawData)) features <- rownames(dataList$rawData)
#         else features <- rownames(dataList$normData)
#         if (is.null(features)) {
#             cli::cli_abort(
#                 "Cannot detect feature names for dataset {.val {d}}."
#             )
#         }
#         ftPassing <- .checkIDIdentical(
#             ref = features,
#             onRow = dataList[c("rawData", "normData")]
#         )
#
#         # Check var features
#         if (!is.null(dataList$V) &&
#             is.null(colnames(dataList$V)) &&
#             !is.null(varFeatures)) {
#                 ## This should not happen but unfortunately, old `V`s might not
#                 ## have var features as their colnames
#                 colnames(dataList$V) <- varFeatures
#         }
#         hvgPassing <- .checkIDIdentical(
#             ref = varFeatures,
#             onCol = dataList[c("scaleData", "V", "U")]
#         )
#
#         # Remove data that has inconsistent information
#         passing <- .combinePassingSignal(names(dataList),
#                                          bcPassing, ftPassing, hvgPassing)
#         dataList <- dataList[passing]
#         for (s in c("scaleData", "H", "V", "U")) {
#             if (!is.null(dataList[[s]])) {
#                 dataList[[s]] <- t(dataList[[s]])
#             }
#         }
#         # 3. Construct ligerDataset objects for each dataset
#         ldList[[d]] <- do.call(createLigerDataset, dataList)
#         if (!all(startsWith(colnames(ldList[[d]]), d)))
#             colnames(ldList[[d]]) <- paste0(d, "_", colnames(ldList[[d]]))
#     }
#     return(ldList)
# }
#
# convertOldLiger.H5 <- function(object, h5FilePath = NULL) {
#     cli::cli_alert_warning(
#         "Please use caution when restoring an H5 based liger object, because
#         the old version does not have any solid restriction on cell/feature
#         identifier matching. New rliger > 1.99 assumes all data was produced
#         with standard old rliger workflow.", wrap = TRUE
#     )
#     dataLists <- list()
#     if (.hasSlot(object, "H")) dataLists$H <- object@H
#     if (.hasSlot(object, "V")) dataLists$V <- object@V
#     if (.hasSlot(object, "U")) dataLists$U <- object@U
#
#     # 1. Deal with cell metadata which establish a correct mapping of cell
#     # barcodes and datasets belonging
#     allDatasets <- Reduce(union, lapply(dataLists, names))
#     cellMeta <- object@cell.data
#     cellMetaDatasets <- unique(as.vector(cellMeta$dataset))
#     if (!identical(sort(allDatasets), sort(cellMetaDatasets))) {
#         # Datasets entry for matrices don't match with cell metadata
#         # Only take the intersection
#         allDatasets <- intersect(allDatasets, cellMetaDatasets)
#         cellMeta <- cellMeta[cellMeta[["dataset"]] %in% allDatasets, ]
#     }
#
#     # Split `dataLists` by dataset
#     datasetLists <- list()
#     for (d in allDatasets) {
#         for (slot in names(dataLists)) {
#             datasetLists[[d]][[slot]] <- dataLists[[slot]][[d]]
#         }
#     }
#
#     # For each existing dataset
#     ldList <- list()
#     for (d in allDatasets) {
#         # "BC" for barcodes
#         # 2. Check and clean up cell barcodes and feature idx issue
#         cellMetaBC <- rownames(cellMeta)[cellMeta$dataset == d]
#         #features <- NULL
#         varFeatures <- object@var.genes
#         dataList <- datasetLists[[d]]
#
#         # Check cell barcodes
#         bcPassing <- .checkIDIdentical(
#             ref = cellMetaBC,
#             #onCol = dataList[c("rawData", "normData")],
#             onRow = dataList[c("H")]
#         )
#
#         # Check raw, norm data features
#         # if (!is.null(dataList$rawData)) features <- rownames(dataList$rawData)
#         # else features <- rownames(dataList$normData)
#         # if (is.null(features)) {
#         #     warning("Cannot detect feature names for dataset \"", d, "\". ",
#         #             "Skipped.")
#         #     next
#         # }
#         # ftPassing <- .checkIDIdentical(
#         #     ref = features,
#         #     onRow = dataList[c("rawData", "normData")]
#         # )
#
#         # Check var features
#         if (!is.null(dataList$V) &&
#             is.null(colnames(dataList$V)) &&
#             !is.null(varFeatures)) {
#             ## This should not happen but unfortunately, old `V`s might not
#             ## have var features as their colnames
#             colnames(dataList$V) <- varFeatures
#         }
#         hvgPassing <- .checkIDIdentical(
#             ref = varFeatures,
#             onCol = dataList[c("V", "U")]
#         )
#
#         # Remove data that has inconsistent information
#         passing <- .combinePassingSignal(names(dataList),
#                                          bcPassing, hvgPassing)
#         dataList <- dataList[passing]
#         for (s in c("H", "V", "U")) {
#             if (!is.null(dataList[[s]])) {
#                 dataList[[s]] <- t(dataList[[s]])
#             }
#         }
#         # 3. Construct H5 ligerDataset objects for each dataset
#         if (!is.null(h5FilePath[[d]])) h5Path <- h5FilePath[[d]]
#         else h5Path <- object@h5file.info[[d]]$file.path
#         if (!hdf5r::is_hdf5(name = h5Path)) {
#             cli::cli_abort("File path for dataset {.val {d}} not found or is not an H5 file: {.file {h5Path}}")
#         }
#         h5Format <- object@h5file.info[[d]]$format.type
#         ldList[[d]] <- do.call(createH5LigerDataset, c(
#             list(h5file = h5Path, formatType = h5Format),
#             dataList
#         ))
#         colnames(ldList[[d]]) <- paste0(d, "_", colnames(ldList[[d]]))
#
#         # 4. Check for potential existing processed result
#         newSlotNameMap <- list(norm.data = "normData",
#                                "scale.data" = "scaleData",
#                                "scale.unshared.data" = "scaleUnsharedData")
#         for (s in c("norm.data", "scale.data", "scale.unshared.data")) {
#             h5file <- getH5File(ldList[[d]])
#             if (h5file$link_exists(s)) {
#                 h5fileInfo(ldList[[d]], newSlotNameMap[[s]], check = FALSE) <- s
#             }
#         }
#     }
#
#     return(ldList)
# }
#
# .checkIDIdentical <- function(
#     ref,
#     onCol = NULL,
#     onRow = NULL
# ) {
#     # ref - a character vector as a reference
#     # onCol - a list of matrix, where the colnames should match to reference
#     # onRow - a list of matrix, where the rownames should match to reference
#     colPassing <- rep(TRUE, length(onCol))
#     names(colPassing) <- names(onCol)
#     for (slot in names(onCol)) {
#         if (is.na(slot)) next
#         if (!identical(colnames(onCol[[slot]]), ref)) {
#             cli::cli_alert_danger("Inconsistent column ID in slot {slot}.")
#             colPassing[slot] <- FALSE
#         }
#     }
#     rowPassing <- rep(TRUE, length(onRow))
#     names(rowPassing) <- names(onRow)
#     for (slot in names(onRow)) {
#         if (is.na(slot)) next
#         if (!identical(rownames(onRow[[slot]]), ref)) {
#             cli::cli_alert_danger("Inconsistent row ID in slot {slot}.")
#             rowPassing[slot] <- FALSE
#         }
#     }
#     return(c(colPassing, rowPassing))
# }
#
# .combinePassingSignal <- function(slotNames, ...) {
#     passings <- list(...)
#     passings <- lapply(passings, function(x) {
#         x <- x[slotNames]
#         names(x) <- slotNames
#         x[is.na(x)] <- TRUE
#         x
#     })
#     Reduce("&", passings)
# }

# .hasSlot <- function(object, name) {
#     tryCatch(
#         expr = {
#             methods::slot(object, name)
#             return(TRUE)
#         },
#         error = function(e) {
#             cli::cli_alert_info("Skipped slot {name} which is not available.")
#             return(FALSE)
#         },
#         warining = function(w) {
#             cli::cli_alert_warning(w)
#             return(FALSE)
#         }
#     )
# } # nocov end
