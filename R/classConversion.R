setClass("ligerDataset")
setClassUnion("matrixLike", c("matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# From other things to ligerDataset class ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert an object of various classes to a ligerDataset object
#' @rdname as.ligerDataset
#' @description S4 method to convert objects of various types to a
#' \linkS4class{ligerDataset} object or a modality specific sub-class of
#' \linkS4class{ligerDataset} object. Supported classes include a matrix like
#' object, a \code{SingleCellExperiment} object, a \code{Seurat} object, and
#' \code{AnnData} object. This method also supports modality setting to a
#' \linkS4class{ligerDataset} inherited object.
#' @param x An object to be converted
#' @param modal The modality of this dataset. Default \code{"default"} for RNA.
#' Can choose from \code{"rna"}, \code{"atac"}.
#' @return A ligerDataset object by default, or a modality specific sub-class of
#' ligerDataset object according to \code{modal}.
#' @export
#' @examples
#' data("pbmc", package = "rliger")
#' ctrl <- dataset(pbmc, "ctrl")
#' ctrl
#' as.ligerDataset(ctrl, modal = "atac")
#' rawCounts <- rawData(ctrl)
#' class(rawCounts)
#' as.ligerDataset(rawCounts)
setGeneric("as.ligerDataset",
           function(x, modal = c("default", "rna", "atac")) {
               standardGeneric("as.ligerDataset")
           })

#' @rdname as.ligerDataset
#' @export
setMethod(
    "as.ligerDataset",
    signature(x = "SingleCellExperiment",
              modal = "ANY"),
    function(x, modal = c("default", "rna", "atac")) {
        if (!requireNamespace("SingleCellExperiment", quietly = "TRUE"))
            stop("Package \"SingleCellExperiment\" needed for this function ",
                 "to work. Please install it by command:\n",
                 "BiocManager::install('SingleCellExperiment')",
                 call. = FALSE)
        if ("counts" %in% SummarizedExperiment::assayNames(x))
            raw.data <- SingleCellExperiment::counts(x)
        else raw.data <- NULL
        if ("logcounts" %in% SummarizedExperiment::assayNames(x))
            normData <- SingleCellExperiment::logcounts(x)
        else normData <- NULL
        if ("counts" %in% SummarizedExperiment::assayNames(x))
            raw.data <- SingleCellExperiment::counts(x)
        createLigerDataset(raw.data = raw.data, normData = normData,
                           modal = modal)
    }
)

#' @rdname as.ligerDataset
#' @export
setMethod(
    "as.ligerDataset",
    signature(x = "ligerDataset",
              modal = "ANY"),
    function(x, modal = c("default", "rna", "atac")) {
        modal <- match.arg(modal)
        newClass <- .modalClassDict[[modal]]
        if (class(x) == newClass) return(x)
        slotFromClass <- methods::slotNames(class(x))
        slotToClass <- methods::slotNames(newClass)
        if (any(!slotFromClass %in% slotToClass))
            warning("Will remove information in the following slots when ",
                    "converting class from `", class(x), "` to `", newClass,
                    "`: ", paste(slotFromClass[!slotFromClass %in% slotToClass],
                                 collapse = ", "))
        newCallArgs <- list(Class = newClass)
        for (s in slotFromClass) {
            if (s %in% slotToClass)
                newCallArgs[[s]] <- methods::slot(x, s)
        }
        do.call("new", newCallArgs)
    }
)

#' @rdname as.ligerDataset
#' @export
setMethod(
    "as.ligerDataset",
    signature(x = "matrixLike",
              modal = "ANY"),
    function(x, modal = c("default", "rna", "atac")) {
        createLigerDataset(x, modal)
    }
)

#' @rdname as.ligerDataset
#' @export
setMethod(
    "as.ligerDataset",
    signature(x = "Seurat",
              modal = "ANY"),
    function(x, modal = c("default", "rna", "atac")) {
        if (!requireNamespace("Seurat", quietly = "TRUE"))
            stop("Package \"Seurat\" needed for this function to work. ",
                 "Please install it by command:\n",
                 "BiocManager::install('Seurat')",
                 call. = FALSE)
        counts <- Seurat::GetAssayData(x, "counts")
        normData <- Seurat::GetAssayData(x, "data")
        if (identical(counts, normData)) normData <- NULL
        scale.data <- Seurat::GetAssayData(x, "scale.data")
        if (sum(dim(scale.data)) == 0) scale.data <- NULL
        createLigerDataset(raw.data = counts, normData = normData,
                           scale.data = scale.data, modal = modal)
    }
)

setClass("anndata._core.anndata.AnnData")
#' @rdname as.ligerDataset
#' @export
setMethod(
    "as.ligerDataset",
    signature(x = "anndata._core.anndata.AnnData",
              modal = "ANY"),
    function(x, modal = c("default", "rna", "atac")) {
        modal <- match.arg(modal)
        message("Python object AnnData input. ")
    }
)

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
#' @export
#' @examples
#' \dontrun{
#' # Suppose you have a liger object of old version (<1.99.0)
#' newLig <- convertOldLiger(oldLig)
#' }
convertOldLiger <- function(
        object,
        dimredName = "tsne.coords",
        clusterName = "clusters"
) {
    ver120 <- package_version("1.99.0")
    if (object@version >= ver120) return(object)
    if (inherits(object@raw.data[[1]], "H5File")) {
        convertOldLiger.H5(object, dimredName = dimredName,
                           clusterName = clusterName)
    } else {
        convertOldLiger.mem(object, dimredName = dimredName,
                            clusterName = clusterName)
    }
}

convertOldLiger.mem <- function(object, dimredName = "tsne.coords",
                                clusterName = "clusters") {
    dataLists <- list(
        rawData = object@raw.data,
        normData = object@norm.data,
        scaleData = object@scale.data,
        H = object@H,
        V = object@V,
        U = object@U
    )
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
                          varFeatures = varFeatures, cellMeta = cellMeta)
    return(newObj)
}

convertOldLiger.H5 <- function(object, dimredName = "tsne_coords",
                               clusterName = "clusters") {
    stop("Not implemented yet")
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
