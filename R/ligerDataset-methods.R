

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
        if (length(dataset) == 0) return(FALSE)
        allCheck <- unlist(lapply(datasets(object)[dataset], isH5Liger))
        return(all(allCheck))
    } else {
        cli::cli_alert_danger("Given object is not of {.cls liger} or {.cls ligerDataset} class.")
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
# S4 Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @param x,object A \code{ligerDataset} object.
#' @param dataset Not applicable for \code{ligerDataset} methods.
#' @param value See detail sections for requirements
#' @param check Whether to perform object validity check on setting new value.
#' @param info Name of the entry in \code{h5fileInfo} slot.
#' @param slot The slot name when using \code{getMatrix}.
#' @param returnList Not applicable for \code{ligerDataset} methods.
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
#' if (requireNamespace("RcppPlanc", quietly = TRUE)) {
#'     pbmc <- runOnlineINMF(pbmc, k = 20, minibatchSize = 100)
#'     ctrl <- dataset(pbmc, "ctrl")
#'     V <- getMatrix(ctrl, "V")
#'     V[1:5, 1:5]
#'     Vs <- getMatrix(pbmc, "V")
#'     length(Vs)
#'     names(Vs)
#'     identical(Vs$ctrl, V)
#' }
setMethod(
    f = "show",
    signature(object = "ligerDataset"),
    definition = function(object) {
        # Use class(object) so that the inheriting classes can be shown properly
        cat("An object of class", class(object), "with",
            ncol(object), "cells\n")
        if (isH5Liger(object) &
            !isTRUE(methods::validObject(object, test = TRUE))) {
            cli::cli_alert_danger("Link to HDF5 file fails. Please try running {.code restoreH5Liger(object)}.")
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
#' @section Subsetting:
#' For more detail of subsetting a \code{liger} object or a
#' \linkS4class{ligerDataset} object, please check out \code{\link{subsetLiger}}
#' and \code{\link{subsetLigerDataset}}. Here, we set the S3 method
#' "single-bracket" \code{[} as a quick wrapper to subset a \code{ligerDataset}
#' object. \code{i} and \code{j} serves as feature and cell subscriptor,
#' respectively, which can be any valid index refering the available features
#' and cells in a dataset. \code{...} arugments are passed to
#' \code{subsetLigerDataset} so that advanced options are allowed.
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
                value[[1L]][match(rownames(x@scaleData), x@rownames)]
        }
        if (!is.null(x@scaleUnsharedData)) {
            colnames(x@scaleUnsharedData) <- value[[2L]]
            rownames(x@scaleUnsharedData) <-
                value[[1L]][match(rownames(x@scaleUnsharedData), x@rownames)]
        }
    }
    if (!is.null(x@H))
        colnames(x@H) <- value[[2L]]
    if (!is.null(x@V))
        rownames(x@V) <- value[[1L]][match(rownames(x@V), x@rownames)]
    if (!is.null(x@B))
        rownames(x@B) <- value[[1L]][match(rownames(x@B), x@rownames)]
    if (!is.null(x@U)) colnames(x@U) <- value[[2L]]
    if ("rawPeak" %in% methods::slotNames(x)) {
        if (!is.null(rawPeak(x))) colnames(rawPeak(x)) <- value[[2L]]
    }
    if ("normPeak" %in% methods::slotNames(x)) {
        if (!is.null(normPeak(x))) colnames(normPeak(x)) <- value[[2L]]
    }
    if ("coordinate" %in% methods::slotNames(x)) {
        if (!is.null(coordinate(x))) rownames(coordinate(x)) <- value[[2L]]
    }
    x@rownames <- value[[1L]]
    x@colnames <- value[[2L]]
    return(x)
})


#' Subset ligerDataset object
#' @name sub-ligerDataset
#' @param x A \linkS4class{ligerDataset} object
#' @param i Numeric, logical index or character vector of feature names to
#' subscribe. Leave missing for all features.
#' @param j Numeric, logical index or character vector of cell IDs to subscribe.
#' Leave missing for all cells.
#' @param ... Additional arguments passed to \code{\link{subsetLigerDataset}}.
#' @export
#' @method [ ligerDataset
#' @return If \code{i} is given, the selected metadata will be returned; if it
#' is missing, the whole cell metadata table in
#' \code{S4Vectors::\link[S4Vectors]{DataFrame}} class will be returned.
#' @examples
#' ctrl <- dataset(pbmc, "ctrl")
#' ctrl[1:5, 1:5]
`[.ligerDataset` <- function(x, i, j, ...) {
    if (missing(i) && missing(j)) {
        return(x)
    }
    if (missing(i)) i <- NULL
    if (missing(j)) j <- NULL
    subsetLigerDataset(x, featureIdx = i, cellIdx = j, ...)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Raw, norm, scale data access ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' @export
#' @rdname ligerDataset-class
setMethod("rawData", "ligerDataset",
          function(x, dataset = NULL) x@rawData)

#' @export
#' @rdname ligerDataset-class
setReplaceMethod(
    "rawData",
    signature(x = "ligerDataset", dataset = "ANY", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, dataset = NULL, check = TRUE, value) {
        if (isH5Liger(x))
            cli::cli_abort("Cannot replace slot with in-memory data for H5 based object.")
        x@rawData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-class
setReplaceMethod(
    "rawData",
    signature(x = "ligerDataset", dataset = "ANY", check = "ANY", value = "H5D"),
    function(x, dataset = NULL, check = TRUE, value) {
        if (!isH5Liger(x))
            cli::cli_abort("Cannot replace slot with on-disk data for in-memory object.")
        x@rawData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)


#' @export
#' @rdname ligerDataset-class
setMethod("normData", "ligerDataset",
          function(x, dataset = NULL) x@normData)

#' @export
#' @rdname ligerDataset-class
setReplaceMethod(
    "normData",
    signature(x = "ligerDataset", dataset = "ANY", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, dataset = NULL, check = TRUE, value) {
        if (isH5Liger(x))
            cli::cli_abort("Cannot replace slot with in-memory data for H5 based object.")
        x@normData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-class
setReplaceMethod(
    "normData",
    signature(x = "ligerDataset", dataset = "ANY", check = "ANY", value = "H5D"),
    function(x, dataset = NULL, check = TRUE, value) {
        if (!isH5Liger(x))
            cli::cli_abort("Cannot replace slot with on-disk data for in-memory object.")
        x@normData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-class
setMethod("scaleData", c("ligerDataset", "missing"),
          function(x, dataset = NULL) x@scaleData)

#' @export
#' @rdname ligerDataset-class
setReplaceMethod(
    "scaleData",
    signature(x = "ligerDataset", dataset = "ANY", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, dataset = NULL, check = TRUE, value) {
        if (isH5Liger(x))
            cli::cli_abort("Cannot replace slot with in-memory data for H5 based object.")
        x@scaleData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-class
setReplaceMethod(
    "scaleData",
    signature(x = "ligerDataset", dataset = "ANY", check = "ANY", value = "H5D"),
    function(x, dataset = NULL, check = TRUE, value) {
        if (!isH5Liger(x))
            cli::cli_abort("Cannot replace slot with on-disk data for in-memory object.")
        x@scaleData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-class
setReplaceMethod(
    "scaleData",
    signature(x = "ligerDataset", dataset = "ANY", check = "ANY", value = "H5Group"),
    function(x, dataset = NULL, check = TRUE, value) {
        if (!isH5Liger(x))
            cli::cli_abort("Cannot replace slot with on-disk data for in-memory object.")
        x@scaleData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-class
setMethod("scaleUnsharedData", c("ligerDataset", "missing"),
          function(x, dataset = NULL) x@scaleUnsharedData)

#' @export
#' @rdname ligerDataset-class
setReplaceMethod(
    "scaleUnsharedData",
    signature(x = "ligerDataset", dataset = "missing", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, check = TRUE, value) {
        if (isH5Liger(x))
            cli::cli_abort("Cannot replace slot with in-memory data for H5 based object.")
        if (!is.null(value)) {
            bc <- x@colnames
            if (!all(endsWith(bc, colnames(value)))) {
                cli::cli_abort(
                    "Column names of {.var value} do not match to those of the object."
                )
            }
            colnames(value) <- bc
        }
        x@scaleUnsharedData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-class
setReplaceMethod(
    "scaleUnsharedData",
    signature(x = "ligerDataset", dataset = "missing", check = "ANY", value = "H5D"),
    function(x, check = TRUE, value) {
        if (!isH5Liger(x))
            cli::cli_abort("Cannot replace slot with on-disk data for in-memory object.")
        x@scaleUnsharedData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname ligerDataset-class
setReplaceMethod(
    "scaleUnsharedData",
    signature(x = "ligerDataset", dataset = "missing", check = "ANY", value = "H5Group"),
    function(x, check = TRUE, value) {
        if (!isH5Liger(x))
            cli::cli_abort("Cannot replace slot with on-disk data for in-memory object.")
        x@scaleUnsharedData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# H5 related ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                    cli::cli_abort(
                        "Specified {.code info} not found: {.val {info[!info %in% names(x@h5fileInfo)]}}"
                    )
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
                cli::cli_abort("{.var info} has to be a single character.")
            if (info %in% c("indicesName", "indptrName", "barcodesName",
                            "genesName", "rawData", "normData",
                            "scaleData")) {
                if (!getH5File(x)$exists(value)) {
                    cli::cli_abort("Specified {.var info} is invalid, {.field info} does not exist in the HDF5 file.")
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Feature metadata ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

#' @export
#' @rdname liger-class
setMethod("varUnsharedFeatures",
          signature(x = "ligerDataset", dataset = "missing"),
          function(x, dataset = NULL) x@varUnsharedFeatures)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "varUnsharedFeatures",
    signature(x = "ligerDataset", dataset = "missing", check = "ANY", value = "character"),
    function(x, dataset = NULL, check = TRUE, value) {
        x@varUnsharedFeatures <- value
        if (isTRUE(check)) {
            if (!all(value %in% rownames(x))) {
                cli::cli_alert_warning("Not all features passed are found.")
            }
        }
        return(x)
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
        cli::cli_alert_warning("Discarding arguments that are not of {.cls ligerDataset} class")
        args <- args[isLD]
    }
    if (!missing(x)) args <- c(list(x), args)
    isH5 <- sapply(args, isH5Liger)
    # See mergeObject.R
    if (all(isH5)) .cbind.ligerDataset.h5(args)
    else if (!any(isH5)) .cbind.ligerDataset.mem(args)
    else
        cli::cli_abort("Cannot {.fn cbind} a hybrid of H5 and in-memory {.cls ligerDataset}s for now.")
}
