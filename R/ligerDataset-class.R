setClassUnion("dgCMatrix_OR_NULL", c("dgCMatrix", "NULL"))
setClassUnion("matrix_OR_NULL", c("matrix", "NULL"))
setClassUnion("matrixLike", c("matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"))
setClassUnion("matrixLike_OR_NULL", c("matrixLike", "NULL"))
# It is quite hard to handle "H5D here, which is indeed defined as an R6 class.
# I'm not sure if this is a proper solution
setOldClass("H5D")
suppressWarnings(setClassUnion("dgCMatrix_OR_H5D_OR_NULL", c("dgCMatrix", "H5D", "NULL")))
setClassUnion("matrix_OR_H5D_OR_NULL", c("matrix", "H5D", "NULL"))

setClassUnion("index",
              members = c("logical", "numeric", "character"))

#' ligerDataset class
#'
#' Object for storing dastaset specific information. Will be embedded within a
#' higher level \linkS4class{liger} object
#' @slot rawData Raw data.
#' @slot normData Normalized data
#' @slot scaleData Scaled data, usually with subset variable features
#' @slot scaleUnsharedData Scaled data of features not shared with other
#' datasets
#' @slot varUnsharedFeatures Variable features not shared with other datasets
#' @slot V matrix
#' @slot A matrix
#' @slot B matrix
#' @slot H matrix
#' @slot U matrix
#' @slot h5fileInfo list
#' @slot featureMeta Feature metadata, DataFrame
#' @slot colnames character
#' @slot rownames character
#' @importClassesFrom S4Vectors DataFrame
#' @exportClass ligerDataset
ligerDataset <- setClass(
    "ligerDataset",
    representation(
        rawData = "dgCMatrix_OR_H5D_OR_NULL",
        normData = "dgCMatrix_OR_H5D_OR_NULL",
        scaleData = "matrix_OR_H5D_OR_NULL",
        scaleUnsharedData = "ANY",
        varUnsharedFeatures = "character",
        H = "matrix_OR_NULL",
        V = "matrix_OR_NULL",
        A = "matrix_OR_NULL",
        B = "matrix_OR_NULL",
        U = "matrix_OR_NULL",
        h5fileInfo = "list",
        featureMeta = "DataFrame",
        colnames = "character",
        rownames = "character"
    )
)

# ------------------------------------------------------------------------------
# Dataset creatinfg function ####
# ------------------------------------------------------------------------------

#' Create in-memory ligerDataset object
#' @param rawData,normData,scaleData A \code{\link[Matrix]{dgCMatrix-class}}
#' object for the raw or normalized expression count or a dense matrix of scaled
#' variable gene expression, respectively. Default \code{NULL} for all three but
#' at lease one has to be specified.
#' @param modal Name of modality for this dataset. Currently options of
#' \code{"default"}, \code{"rna"}, and \code{"atac"} are supported. Default
#' \code{"default"}.
#' @param featureMeta Data frame of feature metadata. Default \code{NULL}.
#' @param ... Additional slot data. See \linkS4class{ligerDataset} for detail.
#' Given values will be directly placed at corresponding slots.
#' @seealso \linkS4class{ligerDataset}, \linkS4class{ligerATACDataset},
#' \code{\link{createH5LigerDataset}}
#' @export
createLigerDataset <- function(
        rawData = NULL,
        modal = c("default", "rna", "atac"),
        normData = NULL,
        scaleData = NULL,
        featureMeta = NULL,
        ...
) {
    modal <- match.arg(modal)
    args <- as.list(environment())
    additional <- list(...)
    # TODO h5 file support
    # Necessary initialization of slots
    if (is.null(rawData) & is.null(normData) & is.null(scaleData)) {
        stop("At least one type of expression data (rawData, normData or ",
             "scaleData) has to be provided")
    }
    # Look for proper colnames and rownames
    cn <- NULL
    rn <- NULL
    for (i in c("rawData", "normData", "scaleData")) {
        cn <- colnames(args[[i]])
        if (!is.null(cn)) break
    }
    if (!is.null(rawData)) {
        rn <- rownames(rawData)
        if (!inherits(rawData, "dgCMatrix"))
            rawData <- methods::as(rawData, "CsparseMatrix")
    }
    if (!is.null(normData)) {
        if (is.null(rn)) rn <- rownames(normData)
        if (!inherits(normData, "dgCMatrix"))
            rawData <- methods::as(normData, "CsparseMatrix")
    }
    if (!is.null(scaleData)) {
        if (is.null(rn)) rn <- rownames(scaleData)
    }
    if (is.null(h5fileInfo)) h5fileInfo <- list()
    if (is.null(featureMeta))
        featureMeta <- S4Vectors::DataFrame(row.names = rn)
    else if (!inherits(featureMeta, "DFrame"))
        featureMeta <- S4Vectors::DataFrame(featureMeta)
    # Create ligerDataset
    allData <- list(.modalClassDict[[modal]],
                    rawData = rawData, normData = normData,
                    scaleData = scaleData, featureMeta = featureMeta,
                    colnames = cn, rownames = rn)
    allData <- c(allData, additional)
    x <- do.call("new", allData)
    return(x)
}

.modalClassDict <- list(
    default = "ligerDataset",
    rna = "ligerRNADataset",
    atac = "ligerATACDataset"
)

.classModalDict <- list(
    ligerDataset = "default",
    ligerRNADataset = "rna",
    ligerATACDataset = "atac"
)

#' Create on-disk ligerDataset Object
#' @param h5file Filename of a H5 file
#' @param formatType Select preset of H5 file structure. Current available
#' options are \code{"10X"} and \code{"AnnData"}.
#' @param rawData,indicesName,indptrName The path in a H5 file for the raw
#' sparse matrix data. These three types of data stands for the \code{x},
#' \code{i}, and \code{p} slots of a \code{\link[Matrix]{dgCMatrix-class}}
#' object. Default \code{NULL} uses \code{formatType} preset.
#' @param normData The path in a H5 file for the "x" vector of the normalized
#' sparse matrix. Default \code{NULL}.
#' @param scaleData The path in a H5 file for the dense 2D scaled matrix.
#' Default \code{NULL}.
#' @param genesName,barcodesName The path in a H5 file for the gene names and
#' cell barcodes. Default \code{NULL} uses \code{formatType} preset.
#' @param modal Name of modality for this dataset. Currently options of
#' \code{"default"}, \code{"rna"}, and \code{"atac"} are supported. Default
#' \code{"default"}.
#' @param featureMeta Data frame for feature metadata. Default \code{NULL}.
#' @param ... Additional slot data. See \linkS4class{ligerDataset} for detail.
#' Given values will be directly placed at corresponding slots.
createH5LigerDataset <- function(
        h5file,
        formatType = NULL,
        rawData = NULL,
        normData = NULL,
        scaleData = NULL,
        barcodesName = NULL,
        genesName = NULL,
        indicesName = NULL,
        indptrName = NULL,
        modal = c("default", "rna", "atac"),
        featureMeta = NULL,
        ...
) {
    if (!hdf5r::is_hdf5(h5file)) {
        stop("Please specify an HDF5 filename to argument `h5file`.")
    }
    modal <- match.arg(modal)
    additional <- list(...)
    h5file <- hdf5r::H5File$new(h5file, mode = "r+")
    if (!is.null(formatType) &&
        formatType %in% c("10X", "AnnData")) {
        if (formatType == "10X") {
            barcodesName <- "matrix/barcodes"
            barcodes <- h5file[[barcodesName]][]
            rawData <- "matrix/data"
            indicesName <- "matrix/indices"
            indptrName <- "matrix/indptr"
            genesName <- "matrix/features/name"
            genes <- h5file[[genesName]][]
        } else if (formatType == "AnnData") {
            barcodesName <- "obs"
            barcodes <- h5file[[barcodesName]][]$cell
            rawData <- "raw.X/data"
            indicesName <- "raw.X/indices"
            indptrName <- "raw.X/indptr"
            genesName <- "raw.var"
            genes <- h5file[[genesName]][]
        } else {
            stop("Specified `formatType` '", formatType,
                 "' is not supported for now.")
        }
    } else {
        barcodes <- h5file[[barcodesName]][]
        genes <- h5file[[genesName]][]
    }
    # The order of list elements matters. Put "paths" together so easier for
    # checking link existence.
    h5.meta <- list(
        H5File = h5file,
        filename = h5file$filename,
        formatType = formatType,
        indicesName = indicesName,
        indptrName = indptrName,
        barcodesName = barcodesName,
        genesName = genesName,
        rawData = rawData,
        normData = normData,
        scaleData = scaleData
    )
    if (!is.null(rawData)) rawData <- h5file[[rawData]]
    if (!is.null(normData)) normData <- h5file[[normData]]
    if (!is.null(scaleData)) scaleData <- h5file[[scaleData]]
    if (is.null(featureMeta))
        featureMeta <- S4Vectors::DataFrame(row.names = genes)
    else if (!inherits(featureMeta, "DFrame"))
        featureMeta <- S4Vectors::DataFrame(featureMeta)
    allData <- list(.modalClassDict[[modal]],
                    rawData = rawData,
                    normData = normData,
                    scaleData = scaleData,
                    #H = H, V = V, A = A, B = B, U = U,
                    h5fileInfo = h5.meta, featureMeta = featureMeta,
                    colnames = barcodes, rownames = genes)
    allData <- c(allData, additional)
    x <- do.call("new", allData)
    x
}

#' Check if a liger or ligerDataset object is made of HDF5 file
#' @param object A liger or ligerDataset object.
#' @param dataset If \code{object} is of liger class, check a specific dataset.
#' If \code{NULL}, Check if all datasets are made of HDF5 file. Default
#' \code{NULL}.
#' @return \code{TRUE} or \code{FALSE} for the specified check.
#' @export
isH5Liger <- function(object, dataset = NULL) {
    if (inherits(object, "ligerDataset")) {
        if (length(object@h5fileInfo) == 0) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    } else if (inherits(object, "liger")) {
        dataset <- .checkUseDatasets(object, dataset)
        allCheck <- unlist(lapply(datasets(object)[dataset], isH5Liger))
        return(all(allCheck))
    } else {
        warning("Given object is not of liger or ligerDataset class.")
        return(FALSE)
    }
}

#' Return preset modality of a ligerDataset object or that of all datasets in a
#' liger object
#' @param object a \linkS4class{ligerDataset} object or a \linkS4class{liger}
#' object
#' @return A single character of modality setting value for
#' \linkS4class{ligerDataset} \code{object}, or a named vector for
#' \linkS4class{liger} object, where the names are dataset names.
#' @export
modalOf <- function(object) {
    if (inherits(object, "ligerDataset")) {
        if (class(object) %in% names(.classModalDict))
            return(.classModalDict[[class(object)]])
        else {
            warning("DEVELOPERS, please add this ligerDataset sub-class to ",
                    "`.classModalDict`", immediate. = TRUE)
            return("UNKNOWN")
        }
    } else if (inherits(object, "liger")) {
        return(sapply(datasets(object), modalOf))
    }
}

# ------------------------------------------------------------------------------
# Object validity ####
# ------------------------------------------------------------------------------

.checkLigerDatasetBarcodes <- function(x) {
    # cell barcodes all consistant
    if (is.null(colnames(x))) {
        return(paste0("No valid cell barcode detected for ligerDataset.\n",
                      "Please create object with matrices with colnames."))
    }
    for (slot in c("rawData", "normData", "scaleData", "scaleUnsharedData",
                   "H", "rawPeak", "normPeak")) {
        if (!slot %in% methods::slotNames(x)) next
        data <- methods::slot(x, slot)
        if (!is.null(data)) {
            barcodes.slot <- colnames(data)
            if (!identical(colnames(x), barcodes.slot)) {
                return(paste0("Inconsistant cell identifiers in ", slot, "."))
            }
        }
    }
    for (slot in c("scaleData", "V")) {
        featuresToCheck <- rownames(methods::slot(x, slot))
        check <- !featuresToCheck %in% rownames(x)
        if (any(check)) {
            msg <- paste0("Features in ", slot, " not found from dataset: ",
                          paste(featuresToCheck[check], collapse = ", "))
            return(msg)
        }
    }
    TRUE
}

.checkH5LigerDatasetLink <- function(x) {
    restoreGuide <- "Please try running `restoreH5Liger(object)`."
    if (!"H5File" %in% names(h5fileInfo(x))) {
        return(paste("`h5fileInfo` incomplete.", restoreGuide))
    }
    h5file <- getH5File(x)
    if (is.null(h5file)) {
        return(paste("`H5File` is NULL in `h5fileInfo` slot.", restoreGuide))
    }
    if (!h5file$is_valid) {
        return(paste("`H5File` is invalid in `h5fileInfo` slot.", restoreGuide))
    }
    if (!is.null(rawData(x))) {
        if (!rawData(x)$is_valid) {
            return(paste("`rawData` slot is invalid.", restoreGuide))
        }
    }
    if (!is.null(normData(x))) {
        if (!normData(x)$is_valid) {
            return(paste("`normData` slot is invalid.", restoreGuide))
        }
    }
    if (!is.null(scaleData(x))) {
        if (!scaleData(x)$is_valid) {
            return(paste("`scaleData` slot is invalid.", restoreGuide))
        }
    }
    TRUE
}

.valid.ligerDataset <- function(object) {
    if (isH5Liger(object)) {
        message("Checking h5 ligerDataset validity")
        .checkH5LigerDatasetLink(object)
    } else {
        message("Checking in memory ligerDataset validity")
        .checkLigerDatasetBarcodes(object)
    }
    # TODO more checks
    # TODO debating on whether to have check of the matching between scaleData
    # features and selected variable features.
}

setValidity("ligerDataset", .valid.ligerDataset)

# ------------------------------------------------------------------------------
# Generics ####
# ------------------------------------------------------------------------------
#' Access data in a liger or ligerDataset object
#' @param x \linkS4class{ligerDataset} object or \linkS4class{liger} object
#' @param name,dataset Name or numeric index of a dataset when \code{x} is a
#' \linkS4class{liger} object, in order to only get the data from this specified
#' dataset. Default \code{NULL}.
#' @param value matrix like or \code{NULL}
#' @param check Whether to perform object validity check on setting new value.
#' @param info Name of the entry in \code{h5fileInfo} slot.
#' @param slot The slot name when using \code{getMatrix}.
#' @param returnList When using \code{getMatrix} and data from one dataset is
#' being returned, whether it should be contained in a list or return the data
#' directly.
#' @export
#' @rdname ligerDataset-access
setGeneric("rawData", function(x) standardGeneric("rawData"))

#' @export
#' @rdname ligerDataset-access
setGeneric("rawData<-", function(x, check = TRUE, value) standardGeneric("rawData<-"))

#' @export
#' @rdname ligerDataset-access
setGeneric("normData", function(x) standardGeneric("normData"))

#' @export
#' @rdname ligerDataset-access
setGeneric("normData<-", function(x, check = TRUE, value) standardGeneric("normData<-"))

#' @export
#' @rdname ligerDataset-access
setGeneric("scaleData", function(x, name = NULL) standardGeneric("scaleData"))

#' @export
#' @rdname ligerDataset-access
setGeneric("scaleData<-", function(x, check = TRUE, value) standardGeneric("scaleData<-"))

#' @export
#' @rdname ligerDataset-access
setGeneric("scaleUnsharedData", function(x, name = NULL) standardGeneric("scaleUnsharedData"))

#' @export
#' @rdname ligerDataset-access
setGeneric("scaleUnsharedData<-", function(x, check = TRUE, value) standardGeneric("scaleUnsharedData<-"))

#' @export
#' @rdname ligerDataset-access
setGeneric("getH5File", function(x, dataset = NULL) standardGeneric("getH5File"))

#' @export
#' @rdname ligerDataset-access
setGeneric("h5fileInfo", function(x, info = NULL) standardGeneric("h5fileInfo"))

#' @export
#' @rdname ligerDataset-access
setGeneric("h5fileInfo<-", function(x, info = NULL, check = TRUE, value) standardGeneric("h5fileInfo<-"))

#' @export
#' @rdname ligerDataset-access
setGeneric("featureMeta", function(x, check = NULL) standardGeneric("featureMeta"))

#' @export
#' @rdname ligerDataset-access
setGeneric("featureMeta<-", function(x, check = TRUE, value) standardGeneric("featureMeta<-"))

#' @export
#' @rdname ligerDataset-access
setGeneric("getMatrix", function(x, slot = "rawData", dataset = NULL, returnList = FALSE) standardGeneric("getMatrix"))
# ------------------------------------------------------------------------------
# Methods ####
# ------------------------------------------------------------------------------

setMethod(
    f = "show",
    signature(object = "ligerDataset"),
    definition = function(object) {
        # Use class(object) so that the inheriting classes can be shown properly
        cat("An object of class", class(object), "with",
            ncol(object), "cells\n")
        if (isH5Liger(object) &
            !isTRUE(methods::validObject(object, test = TRUE))) {
            warning("Link to HDF5 file fails. Please try running ",
                    "`restoreH5Liger(object)`.")
            return()
        }
        for (slot in c("rawData", "normData", "scaleData",
                       "scaleUnsharedData")) {
            data <- methods::slot(object, slot)
            if (!is.null(data)) {
                if (inherits(data, c("matrix", "dgCMatrix",
                                     "dgTMatrix", "dgeMatrix"))) {
                    cat(paste0(slot, ":"), nrow(data), "features\n")
                }
                if (inherits(data, "H5D")) {
                    cat(paste0(slot, ":"), paste(data$dims, collapse = " x "),
                        "values in H5D object\n")
                }
            }
        }
        # Information for sub-classes added below, in condition statements
        if ("rawPeak" %in% methods::slotNames(object)) {
            if (!is.null(rawPeak(object)))
                cat("rawPeak:", nrow(rawPeak(object)), "regions\n")
            if (!is.null(normPeak(object)))
                cat("normPeak:", nrow(normPeak(object)), "regions\n")
        }

        invisible(x = NULL)
    }
)

#' @rdname dim-liger
#' @export
setMethod("dim", "ligerDataset", function(x) {
    nr <- length(x@rownames)
    nc <- length(x@colnames)
    c(nr, nc)
})

#' @rdname dim-liger
#' @export
setMethod("dimnames", "ligerDataset", function(x) {
    rn <- x@rownames
    cn <- x@colnames
    list(rn, cn)
})

#' @rdname dim-liger
#' @export
setReplaceMethod("dimnames", c("ligerDataset", "list"), function(x, value) {
    if (!isH5Liger(x)) {
        if (!is.null(rawData(x))) {
            rownames(x@rawData) <- value[[1L]]
            colnames(x@rawData) <- value[[2L]]
        }
        if (!is.null(normData(x))) {
            rownames(x@normData) <- value[[1L]]
            colnames(x@normData) <- value[[2L]]
        }
        if (!is.null(scaleData(x))) {
            colnames(x@scaleData) <- value[[2L]]
            rownames(x@scaleData) <-
                value[[1L]][x@rownames %in% rownames(x@scaleData)]
        }
        if (!is.null(x@scaleUnsharedData)) {
            colnames(x@scaleUnsharedData) <- value[[2L]]
            rownames(x@scaleUnsharedData) <-
                value[[1L]][x@rownames %in% rownames(x@scaleUnsharedData)]
        }
    }
    # else {
        #h5file <- getH5File(x)
        #h5meta <- h5fileInfo(x)
        #h5file[[h5meta$genesName]][1:length(value[[1L]])] <- value[[1L]]
        #h5file[[h5meta$barcodesName]][1:length(value[[2L]])] <- value[[2L]]
    #}
    if (!is.null(x@H))
        rownames(x@H) <- value[[2L]]
    if (!is.null(x@V))
        colnames(x@V) <- value[[1L]][x@rownames %in% rownames(x@V)]
    if (!is.null(x@B))
        rownames(x@B) <- value[[1L]][x@rownames %in% rownames(x@V)]
    if (!is.null(x@U)) colnames(x@U) <- value[[2L]]
    if ("rawPeak" %in% methods::slotNames(x)) {
        if (!is.null(rawPeak(x))) colnames(rawPeak(x)) <- value[[2L]]
    }
    if ("normPeak" %in% methods::slotNames(x)) {
        if (!is.null(normPeak(x))) colnames(normPeak(x)) <- value[[2L]]
    }
    x@rownames <- value[[1L]]
    x@colnames <- value[[2L]]
    return(x)
})

#' @export
#' @rdname subsetLigerDataset
setMethod(
    "[",
    signature(x = "ligerDataset", i = "index", j = "missing"),
    function(x, i, j, ...)
        subsetLigerDataset(x, featureIdx = i, cellIdx = NULL, ...)
)

#' @export
#' @rdname subsetLigerDataset
setMethod(
    "[",
    signature(x = "ligerDataset", i = "missing", j = "index"),
    function(x, i, j, ...)
        subsetLigerDataset(x, featureIdx = NULL, cellIdx = j, ...)
)

#' @export
#' @rdname subsetLigerDataset
setMethod(
    "[",
    signature(x = "ligerDataset", i = "index", j = "index"),
    function(x, i, j, ...)
        subsetLigerDataset(x, featureIdx = i, cellIdx = j, ...)
)

#' @export
#' @rdname ligerDataset-access
setMethod("rawData", "ligerDataset",
          function(x) x@rawData)

#' @export
#' @rdname ligerDataset-access
setReplaceMethod(
    "rawData",
    signature(x = "ligerDataset", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, check = TRUE, value) {
        if (isH5Liger(x))
            stop("Cannot replace slot with in-memory data for H5 based object.")
        x@rawData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-access
setReplaceMethod(
    "rawData",
    signature(x = "ligerDataset", check = "ANY", value = "H5D"),
    function(x, check = TRUE, value) {
        if (!isH5Liger(x))
            stop("Cannot replace slot with on-disk data for in-memory object.")
        x@rawData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-access
setMethod("normData", "ligerDataset",
          function(x) x@normData)

#' @export
#' @rdname ligerDataset-access
setReplaceMethod(
    "normData",
    signature(x = "ligerDataset", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, check = TRUE, value) {
        if (isH5Liger(x))
            stop("Cannot replace slot with in-memory data for H5 based object.")
        x@normData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-access
setReplaceMethod(
    "normData",
    signature(x = "ligerDataset", check = "ANY", value = "H5D"),
    function(x, check = TRUE, value) {
        if (!isH5Liger(x))
            stop("Cannot replace slot with on-disk data for in-memory object.")
        x@normData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-access
setMethod("scaleData", "ligerDataset",
          function(x, name = NULL) x@scaleData)

#' @export
#' @rdname ligerDataset-access
setMethod(
    "scaleData",
    signature(x = "liger", name = "character"),
    function(x, name) {
        scaleData(dataset(x, name))
    }
)

#' @export
#' @rdname ligerDataset-access
setMethod(
    "scaleData",
    signature(x = "liger", name = "numeric"),
    function(x, name) {
        scaleData(dataset(x, name))
    }
)

#' @export
#' @rdname ligerDataset-access
setReplaceMethod(
    "scaleData",
    signature(x = "ligerDataset", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, check = TRUE, value) {
        if (isH5Liger(x))
            stop("Cannot replace slot with in-memory data for H5 based object.")
        x@scaleData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-access
setReplaceMethod(
    "scaleData",
    signature(x = "ligerDataset", check = "ANY", value = "H5D"),
    function(x, check = TRUE, value) {
        if (!isH5Liger(x))
            stop("Cannot replace slot with on-disk data for in-memory object.")
        x@scaleData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)


#' @export
#' @rdname ligerDataset-access
setMethod("scaleUnsharedData", "ligerDataset",
          function(x, name = NULL) x@scaleUnsharedData)

#' @export
#' @rdname ligerDataset-access
setMethod(
    "scaleUnsharedData",
    signature(x = "liger", name = "character"),
    function(x, name) {
        scaleUnsharedData(dataset(x, name))
    }
)

#' @export
#' @rdname ligerDataset-access
setMethod(
    "scaleUnsharedData",
    signature(x = "liger", name = "numeric"),
    function(x, name) {
        scaleUnsharedData(dataset(x, name))
    }
)

#' @export
#' @rdname ligerDataset-access
setReplaceMethod(
    "scaleUnsharedData",
    signature(x = "ligerDataset", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, check = TRUE, value) {
        if (isH5Liger(x))
            stop("Cannot replace slot with in-memory data for H5 based object.")
        x@scaleUnsharedData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-access
setReplaceMethod(
    "scaleUnsharedData",
    signature(x = "ligerDataset", check = "ANY", value = "H5D"),
    function(x, check = TRUE, value) {
        if (!isH5Liger(x))
            stop("Cannot replace slot with on-disk data for in-memory object.")
        x@scaleUnsharedData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-access
setMethod("getH5File",
          signature = signature(x = "ligerDataset", dataset = "missing"),
          function(x, dataset = NULL) h5fileInfo(x, "H5File"))

#' @export
#' @rdname ligerDataset-access
setMethod("getH5File",
          signature = signature(x = "liger", dataset = "character"),
          function(x, dataset = NULL) {
              if (is.null(dataset)) dataset <- names(x)
              if (any(!dataset %in% names(x))) {
                  stop("Specified dataset name(s) not found: ",
                       paste(dataset[!dataset %in% names(x)], collapse = ", "))
              }
              results <- lapply(datasets(x)[dataset],
                                function(ld) h5fileInfo(ld, "H5File"))
              if (length(results) == 1) results <- results[[1]]
              results
          })

#' @export
#' @rdname ligerDataset-access
setMethod(
    "h5fileInfo",
    signature = signature(x = "ligerDataset", info = "ANY"),
    function(x, info = NULL) {
        if (is.null(info)) result <- x@h5fileInfo
        else {
            if (length(info) == 1) result <- x@h5fileInfo[[info]]
            else {
                if (any(!info %in% names(x@h5fileInfo))) {
                    stop("Specified h5file info not found: ",
                         paste(info[!info %in% names(x@h5fileInfo)],
                               collapse = ", "))
                }
                result <- x@h5fileInfo[info]
                names(result) <- info
            }
        }
        return(result)
    })

#' @export
#' @rdname ligerDataset-access
setReplaceMethod(
    "h5fileInfo",
    signature = signature(
        x = "ligerDataset",
        info = "ANY",
        check = "ANY",
        value = "ANY"
    ),
    function(x, info = NULL, check = TRUE, value) {
        if (is.null(info)) {
            x@h5fileInfo <- value
        } else {
            if (!is.character(info) | length(info) != 1)
                stop('`info` has to be a single character.')
            if (info %in% c('indicesName', 'indptrName', 'barcodesName',
                            'genesName', 'rawData', 'normData',
                            'scaleData')) {
                if (!getH5File(x)$exists(value)) {
                    stop("Specified info is invalid, '", info,
                         "' does not exists in the HDF5 file.")
                }
            }
            x@h5fileInfo[[info]] <- value
            if (info %in% c('rawData', 'normData', 'scaleData',
                            'scaleUnsharedData')) {
                x <- do.call(paste0(info, "<-"),
                             list(x = x,
                                  value = getH5File(x)[[h5fileInfo(x, info)]],
                                  check = check))
            }
        }
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-access
setMethod("featureMeta", signature(x = "ligerDataset", check = "ANY"),
          function(x, check = NULL) {
    x@featureMeta
})

#' @export
#' @rdname ligerDataset-access
setReplaceMethod(
    "featureMeta",
    signature(x = "ligerDataset", check = "ANY"),
    function(x, check = TRUE, value) {
        if (!inherits(value, "DFrame"))
            value <- S4Vectors::DataFrame(value)
        x@featureMeta <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-access
setMethod(
    "getMatrix", signature(x = "ligerDataset", dataset = "missing",
                           returnList = "missing"),
    function(x,
             slot = c("rawData", "normData", "scaleData",
                      "scaleUnsharedData", "H", "V", "U"),
             dataset = NULL) {
        # TODO: Currently directly find the data with slot, but need to
        # think about maintainability when we need to change slot name.
        slot <- match.arg(slot)
        methods::slot(x, slot)
    })

#' @export
#' @rdname ligerDataset-access
setMethod(
    "getMatrix", signature(x = "liger"),
    function(x,
             slot = c("rawData", "normData", "scaleData",
                      "scaleUnsharedData", "H", "V", "U", "W", "H.norm"),
             dataset = NULL,
             returnList = FALSE) {
        slot <- match.arg(slot)
        if (slot == "W") return(x@W)
        if (slot == "H.norm") return(x@H.norm)
        if (is.null(dataset)) {
            return(lapply(datasets(x), function(ld) getMatrix(ld, slot)))
        } else {
            if (length(dataset) == 1) {
                if (isTRUE(returnList)) {
                    result <- list(getMatrix(dataset(x, dataset), slot))
                    names(result) <- dataset
                    return(result)
                } else if (isFALSE(returnList))
                    return(getMatrix(dataset(x, dataset), slot))
            } else {
                lds <- datasets(x)[dataset]
                return(lapply(lds, function(ld) getMatrix(ld, slot)))
            }
        }
    })

#' cbind method for ligerDataset objects
#' @description The list of \code{datasets} slot, the rows of
#' \code{cellMeta} slot and the list of \code{commands} slot will be simply
#' concatenated. Variable features in \code{varFeatures} slot will be taken a
#' union. The \eqn{W} and \eqn{H.norm} matrices are not taken into account for
#' now.
#' @param x \linkS4class{ligerDataset} object
#' @param ... \linkS4class{ligerDataset} objects
#' @param deparse.level Not used here.
#' @return A new \linkS4class{liger} object containing all datasets in the order
#' of specification.
#' @name cbind.ligerDataset
NULL

#' @rdname cbind.ligerDataset
#' @export
#' @method cbind ligerDataset
cbind.ligerDataset <- function(x, ...,
                               deparse.level = 1) {
    args <- list(...)
    isLD <- sapply(args, function(x) inherits(x, "ligerDataset"))
    if (any(!isLD)) {
        warning("Discarding arguments that are not of ligerDataset class")
        args <- args[isLD]
    }
    if (!missing(x)) args <- c(list(x), args)
    isH5 <- sapply(args, isH5Liger)
    if (all(isH5)) .cbind.ligerDataset.h5(args)
    else if (!any(isH5)) .cbind.ligerDataset.mem(args)
    else
        stop("Cannot `cbind` a hybrid of H5 ligerDatasets and ",
             "in-memory ligerDatasets for now.")
}
