setClassUnion("dgCMatrix_OR_NULL", c("dgCMatrix", "NULL"))
setClassUnion("matrix_OR_NULL", c("matrix", "NULL"))
setClassUnion("matrixLike", c("matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"))
setClassUnion("matrixLike_OR_NULL", c("matrixLike", "NULL"))
# It is quite hard to handle "H5D here, which is indeed defined as an R6 class.
# I'm not sure if this is a proper solution
setOldClass("H5D")
setOldClass("H5Group")
suppressWarnings(setClassUnion("dgCMatrix_OR_H5D_OR_NULL", c("dgCMatrix", "H5D", "NULL")))
setClassUnion("matrix_OR_H5D_OR_NULL", c("matrix", "H5D", "NULL"))
setClassUnion("matrixLike_OR_H5D_OR_H5Group_OR_NULL", c("matrixLike", "H5D", "H5Group", "NULL"))
setClassUnion("index",
              members = c("logical", "numeric", "character"))

#' ligerDataset class
#'
#' Object for storing dastaset specific information. Will be embedded within a
#' higher level \linkS4class{liger} object
#' @docType class
#' @rdname ligerDataset-class
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
        scaleData = "matrixLike_OR_H5D_OR_H5Group_OR_NULL",
        scaleUnsharedData = "matrixLike_OR_H5D_OR_H5Group_OR_NULL",
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dataset creatinfg function ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check if a liger or ligerDataset object is made of HDF5 file
#' @param object A liger or ligerDataset object.
#' @param dataset If \code{object} is of liger class, check a specific dataset.
#' If \code{NULL}, Check if all datasets are made of HDF5 file. Default
#' \code{NULL}.
#' @return \code{TRUE} or \code{FALSE} for the specified check.
#' @export
#' @examples
#' isH5Liger(pbmc)
#' isH5Liger(pbmc, "ctrl")
#' ctrl <- dataset(pbmc, "ctrl")
#' isH5Liger(ctrl)
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
#' @examples
#' modalOf(pbmc)
#' ctrl <- dataset(pbmc, "ctrl")
#' modalOf(ctrl)
#' ctrl.atac <- as.ligerDataset(ctrl, modal = "atac")
#' modalOf(ctrl.atac)
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Object validity ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.checkLigerDatasetBarcodes <- function(x) {
    # TODO: Functions should have cyclomatic complexity of less than 15,
    # this has 17
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
        # message("Checking h5 ligerDataset validity")
        .checkH5LigerDatasetLink(object)
    } else {
        # message("Checking in memory ligerDataset validity")
        .checkLigerDatasetBarcodes(object)
    }
    # TODO more checks
    # TODO debating on whether to have check of the matching between scaleData
    # features and selected variable features.
}

setValidity("ligerDataset", .valid.ligerDataset)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S4 Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @param x,object A \code{ligerDataset} object.
#' @param dataset Not applicable for \code{ligerDataset} methods.
#' @param value See detail sections for requirements
#' @param check Whether to perform object validity check on setting new value.
#' @param info Name of the entry in \code{h5fileInfo} slot.
#' @param slot The slot name when using \code{getMatrix}.
#' @param returnList Not applicable for \code{ligerDataset} methods.
#' @param i,j Feature and cell index for \code{`[`} method. For \code{`[[`}
#' method, use a single variable name with \code{i} and \code{j} is not
#' applicable.
#' @param drop Not applicable.
#' @param ... See detailed sections for explanation.
#' @rdname ligerDataset-class
#' @export
#' @examples
#' ctrl <- dataset(pbmc, "ctrl")
#'
#' # Methods for base generics
#' ctrl
#' print(ctrl)
#' dim(ctrl)
#' ncol(ctrl)
#' nrow(ctrl)
#' colnames(ctrl)[1:5]
#' rownames(ctrl)[1:5]
#' ctrl[1:5, 1:5]
#'
#' # rliger generics
#' ## raw data
#' m <- rawData(ctrl)
#' class(m)
#' dim(m)
#' ## normalized data
#' pbmc <- normalize(pbmc)
#' ctrl <- dataset(pbmc, "ctrl")
#' m <- normData(ctrl)
#' class(m)
#' dim(m)
#' ## scaled data
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
#' ctrl <- dataset(pbmc, "ctrl")
#' m <- scaleData(ctrl)
#' class(m)
#' dim(m)
#' n <- scaleData(pbmc, "ctrl")
#' identical(m, n)
#' ## Any other matrices
#' pbmc <- online_iNMF(pbmc, k = 20, miniBatch_size = 100)
#' ctrl <- dataset(pbmc, "ctrl")
#' V <- getMatrix(ctrl, "V")
#' V[1:5, 1:5]
#' Vs <- getMatrix(pbmc, "V")
#' length(Vs)
#' names(Vs)
#' identical(Vs$ctrl, V)
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
                    if (length(data$dims) == 1) {
                        cat(paste0(slot, ":"), data$dims,
                            "non-zero values in H5D object\n")
                    } else {
                        cat(paste0(slot, ":"),
                            paste(data$dims, collapse = " x "),
                            "values in H5D object\n")
                    }
                }
                if (inherits(data, "H5Group")) {
                    cat(paste0(slot, ":"), data[["data"]]$dims,
                        "non-zero values in H5Group object\n")
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

#' @section Dimensionality:
#' For a \code{ligerDataset} object, the column orientation is assigned for
#' cells and rows are for features. Therefore, for \code{ligerDataset} objects,
#' \code{dim()} returns a numeric vector of two numbers which are number of
#' features and number of cells. \code{dimnames()} returns a list of two
#' character vectors, which are the feature names and the cell barcodes.
#'
#' For direct call of \code{dimnames<-} method, \code{value} should be a list
#' with a character vector of feature names as the first element and cell
#' identifiers as the second element. For \code{colnames<-} method, the
#' character vector of cell identifiers. For \code{rownames<-} method, the
#' character vector of feature names.
#' @rdname ligerDataset-class
#' @export
setMethod("dim", "ligerDataset", function(x) {
    nr <- length(x@rownames)
    nc <- length(x@colnames)
    c(nr, nc)
})

#' @rdname ligerDataset-class
#' @export
setMethod("dimnames", "ligerDataset", function(x) {
    rn <- x@rownames
    cn <- x@colnames
    list(rn, cn)
})

#' @rdname ligerDataset-class
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
                value[[1L]][.getVarFtIdx(x@rownames, rownames(x@scaleData))]
        }
        if (!is.null(x@scaleUnsharedData)) {
            colnames(x@scaleUnsharedData) <- value[[2L]]
            rownames(x@scaleUnsharedData) <-
                value[[1L]][.getVarFtIdx(x@rownames,
                                         rownames(x@scaleUnsharedData))]
        }
    }
    if (!is.null(x@H))
        colnames(x@H) <- value[[2L]]
    if (!is.null(x@V))
        rownames(x@V) <- value[[1L]][.getVarFtIdx(x@rownames, rownames(x@V))]
    if (!is.null(x@B))
        rownames(x@B) <- value[[1L]][.getVarFtIdx(x@rownames, rownames(x@B))]
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

.getVarFtIdx <- function(full, var) {
    # full - character vector of ligerDataset rownames
    # var - character vector of var features (in scaleData, V), might not follow
    # the original order of `full`
    # return numeric value that select ordered corresponding replacement from
    # `dimnames<-`'s value[[1]]
    fullNamedIdx <- seq_along(full)
    names(fullNamedIdx) <- full
    varNumIdx <- fullNamedIdx[var]
    names(varNumIdx) <- NULL
    return(varNumIdx)
}

#' @section Subsetting:
#' For more detail of subsetting a \code{liger} object or a
#' \linkS4class{ligerDataset} object, please check out \code{\link{subsetLiger}}
#' and \code{\link{subsetLigerDataset}}. Here, we set the S4 method
#' "single-bracket" \code{[} as a quick wrapper to subset a \code{ligerDataset}
#' object. \code{i} and \code{j} serves as feature and cell subscriptor,
#' respectively, which can be any valid index refering the available features
#' and cells in a dataset. \code{...} arugments are passed to
#' \code{subsetLigerDataset} so that advanced options are allowed.
#' @export
#' @rdname ligerDataset-class
setMethod(
    "[",
    signature(x = "ligerDataset", i = "index", j = "missing"),
    function(x, i, j, ...) {
        subsetLigerDataset(x, featureIdx = i, cellIdx = NULL, ...)
    }
)

#' @export
#' @rdname ligerDataset-class
setMethod(
    "[",
    signature(x = "ligerDataset", i = "missing", j = "index"),
    function(x, i, j, ...) {
        subsetLigerDataset(x, featureIdx = NULL, cellIdx = j, ...)
    }
)

#' @export
#' @rdname ligerDataset-class
setMethod(
    "[",
    signature(x = "ligerDataset", i = "index", j = "index"),
    function(x, i, j, ...) {
        subsetLigerDataset(x, featureIdx = i, cellIdx = j, ...)
    }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Raw, norm, scale data access ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @section Matrix access:
#' For \code{ligerDataset} object, \code{rawData()}, \code{normData},
#' \code{scaleData()} and \code{scaleUnsharedData()} methods are exported for
#' users to access the corresponding feature expression matrix. Replacement
#' methods are also available to modify the slots.
#'
#' For other matrices, such as the \eqn{H} and \eqn{V}, which are dataset
#' specific, please use \code{getMatrix()} method with specifying slot name.
#' Directly accessing slot with \code{@} is generally not recommended.
#' @export
#' @rdname ligerDataset-class
setGeneric("rawData", function(x) standardGeneric("rawData"))

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "rawData<-",
    function(x, check = TRUE, value) standardGeneric("rawData<-")
)

#' @export
#' @rdname ligerDataset-class
setMethod("rawData", "ligerDataset",
          function(x) x@rawData)

#' @export
#' @rdname ligerDataset-class
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
#' @rdname ligerDataset-class
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
#' @rdname ligerDataset-class
setGeneric("normData", function(x) standardGeneric("normData"))

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "normData<-",
    function(x, check = TRUE, value) standardGeneric("normData<-")
)

#' @export
#' @rdname ligerDataset-class
setMethod("normData", "ligerDataset",
          function(x) x@normData)

#' @export
#' @rdname ligerDataset-class
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
#' @rdname ligerDataset-class
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
#' @rdname ligerDataset-class
setGeneric(
    "scaleData",
    function(x, dataset = NULL) standardGeneric("scaleData")
)

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "scaleData<-",
    function(x, check = TRUE, value) standardGeneric("scaleData<-")
)

#' @export
#' @rdname ligerDataset-class
setMethod("scaleData", c("ligerDataset", "missing"),
          function(x, dataset = NULL) x@scaleData)

#' @export
#' @rdname liger-class
setMethod(
    "scaleData",
    signature(x = "liger", dataset = "character"),
    function(x, dataset) {
        scaleData(dataset(x, dataset))
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "scaleData",
    signature(x = "liger", dataset = "numeric"),
    function(x, dataset) {
        scaleData(dataset(x, dataset))
    }
)

#' @export
#' @rdname ligerDataset-class
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
#' @rdname ligerDataset-class
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
#' @rdname ligerDataset-class
setReplaceMethod(
    "scaleData",
    signature(x = "ligerDataset", check = "ANY", value = "H5Group"),
    function(x, check = TRUE, value) {
        if (!isH5Liger(x))
            stop("Cannot replace slot with on-disk data for in-memory object.")
        x@scaleData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "scaleUnsharedData",
    function(x, dataset = NULL) standardGeneric("scaleUnsharedData")
)

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "scaleUnsharedData<-",
    function(x, check = TRUE, value) standardGeneric("scaleUnsharedData<-")
)

#' @export
#' @rdname ligerDataset-class
setMethod("scaleUnsharedData", c("ligerDataset", "missing"),
          function(x, dataset = NULL) x@scaleUnsharedData)

#' @export
#' @rdname liger-class
setMethod(
    "scaleUnsharedData",
    signature(x = "liger", dataset = "character"),
    function(x, dataset) {
        scaleUnsharedData(dataset(x, dataset))
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "scaleUnsharedData",
    signature(x = "liger", dataset = "numeric"),
    function(x, dataset) {
        scaleUnsharedData(dataset(x, dataset))
    }
)

#' @export
#' @rdname ligerDataset-class
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
#' @rdname ligerDataset-class
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
#' @rdname ligerDataset-class
setReplaceMethod(
    "scaleUnsharedData",
    signature(x = "ligerDataset", check = "ANY", value = "H5Group"),
    function(x, check = TRUE, value) {
        if (!isH5Liger(x))
            stop("Cannot replace slot with on-disk data for in-memory object.")
        x@scaleUnsharedData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "getMatrix",
    function(x, slot = "rawData", dataset = NULL, returnList = FALSE) {
        standardGeneric("getMatrix")
    }
)

#' @export
#' @rdname ligerDataset-class
setMethod(
    "getMatrix", signature(x = "ligerDataset", dataset = "missing",
                           returnList = "missing"),
    function(x,
             slot = c("rawData", "normData", "scaleData",
                      "scaleUnsharedData", "H", "V", "U", "A", "B"),
             dataset = NULL) {
        # TODO: Currently directly find the data with slot, but need to
        # think about maintainability when we need to change slot name.
        slot <- match.arg(slot)
        methods::slot(x, slot)
    })

#' @export
#' @rdname liger-class
setMethod(
    "getMatrix", signature(x = "liger"),
    function(x,
             slot = c("rawData", "normData", "scaleData",
                      "scaleUnsharedData", "H", "V", "U", "A", "B",
                      "W", "H.norm"),
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# H5 related ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @section H5 file and information access:
#' A \code{ligerDataset} object has a slot called \code{h5fileInfo}, which is a
#' list object. The first element is called \code{$H5File}, which is an
#' \code{H5File} class object and is the connection to the input file. The
#' second element is \code{$filename} which stores the absolute path of the H5
#' file in the current machine. The third element \code{$formatType} stores the
#' name of preset being used, if applicable. The other following keys pair with
#' paths in the H5 file that point to specific data for constructing a feature
#' expression matrix.
#'
#' \code{h5fileInfo()} method access the list described above and simply
#' retrieves the corresponding value. When \code{info = NULL}, returns the whole
#' list. When \code{length(info) == 1}, returns the requested list value. When
#' more info requested, returns a subset list.
#'
#' The replacement method modifies the list elements and corresponding slot
#' value (if applicable) at the same time. For example, running
#' \code{h5fileInfo(obj, "rawData") <- newPath} not only updates the list, but
#' also updates the \code{rawData} slot with the \code{H5D} class data at
#' "newPath" in the \code{H5File} object.
#'
#' \code{getH5File()} is a wrapper and is equivalent to
#' \code{h5fileInfo(obj, "H5File")}.
#' @export
#' @rdname ligerDataset-class
setGeneric("h5fileInfo", function(x, info = NULL) standardGeneric("h5fileInfo"))

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "h5fileInfo<-",
    function(x, info = NULL, check = TRUE, value) {
        standardGeneric("h5fileInfo<-")
    }
)

#' @export
#' @rdname ligerDataset-class
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
#' @rdname ligerDataset-class
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
                stop("`info` has to be a single character.")
            if (info %in% c("indicesName", "indptrName", "barcodesName",
                            "genesName", "rawData", "normData",
                            "scaleData")) {
                if (!getH5File(x)$exists(value)) {
                    stop("Specified info is invalid, '", info,
                         "' does not exists in the HDF5 file.")
                }
            }
            x@h5fileInfo[[info]] <- value
            if (info %in% c("rawData", "normData", "scaleData",
                            "scaleUnsharedData")) {
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
#' @rdname ligerDataset-class
setGeneric("getH5File", function(x, dataset = NULL) standardGeneric("getH5File"))

#' @export
#' @rdname ligerDataset-class
setMethod("getH5File",
          signature = signature(x = "ligerDataset", dataset = "missing"),
          function(x, dataset = NULL) h5fileInfo(x, "H5File"))

#' @export
#' @rdname liger-class
setMethod("getH5File",
          signature = signature(x = "liger", dataset = "ANY"),
          function(x, dataset = NULL) {
              if (is.null(dataset)) dataset <- names(x)
              dataset <- .checkUseDatasets(x, dataset)
              results <- lapply(datasets(x)[dataset],
                                function(ld) h5fileInfo(ld, "H5File"))
              if (length(results) == 1) results <- results[[1]]
              results
          })

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Feature metadata ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @section Feature metadata access:
#' A slot \code{featureMeta} is included for each \code{ligerDataset} object.
#' This slot requires a \code{\link[S4Vectors]{DataFrame-class}} object, which
#' is the same as \code{cellMeta} slot of a \linkS4class{liger} object. However,
#' the associated S4 methods only include access to the whole table for now.
#' Internal information access follows the same way as data.frame operation.
#' For example, \code{featureMeta(ligerD)$nCell} or
#' \code{featureMeta(ligerD)[varFeatures(ligerObj), "gene_var"]}.
#' @export
#' @rdname ligerDataset-class
setGeneric("featureMeta", function(x, check = NULL) {
    standardGeneric("featureMeta")
})

#' @export
#' @rdname ligerDataset-class
setGeneric("featureMeta<-", function(x, check = TRUE, value) {
    standardGeneric("featureMeta<-")
})

#' @export
#' @rdname ligerDataset-class
setMethod("featureMeta", signature(x = "ligerDataset", check = "ANY"),
          function(x, check = NULL) {
    x@featureMeta
})

#' @export
#' @rdname ligerDataset-class
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 Method ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @section Concatenate ligerDataset:
#' \code{cbind()} method is implemented for concatenating \code{ligerDataset}
#' objects by cells. When applying, all feature expression matrix will be merged
#' with taking a union of all features for the rows.
#' @param deparse.level Not used here.
#' @rdname ligerDataset-class
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
    # See mergeObject.R
    if (all(isH5)) .cbind.ligerDataset.h5(args)
    else if (!any(isH5)) .cbind.ligerDataset.mem(args)
    else
        stop("Cannot `cbind` a hybrid of H5 ligerDatasets and ",
             "in-memory ligerDatasets for now.")
}
