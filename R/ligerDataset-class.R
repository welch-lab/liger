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
#' @slot raw.data Raw data.
#' @slot norm.data Normalized data
#' @slot scale.data Scaled data, usually with subset variable features
#' @slot scale.unshared.data Scaled data of features not shared with other
#' datasets
#' @slot var.unshared.features Variable features not shared with other datasets
#' @slot V matrix
#' @slot A matrix
#' @slot B matrix
#' @slot U matrix
#' @slot agg.data list
#' @slot h5file.info list
#' @slot feature.meta Feature metadata, DataFrame
#' @importClassesFrom S4Vectors DataFrame
#' @exportClass ligerDataset
ligerDataset <- setClass(
    "ligerDataset",
    representation(
        raw.data = "dgCMatrix_OR_H5D_OR_NULL",
        norm.data = "dgCMatrix_OR_H5D_OR_NULL",
        scale.data = "matrix_OR_H5D_OR_NULL",
        scale.unshared.data = "ANY",
        var.unshared.features = "character",
        H = "matrix_OR_NULL",
        V = "matrix_OR_NULL",
        A = "matrix_OR_NULL",
        B = "matrix_OR_NULL",
        U = "matrix_OR_NULL",
        agg.data = "list",
        h5file.info = "list",
        feature.meta = "DataFrame",
        colnames = "character",
        rownames = "character"
    )
)

# ------------------------------------------------------------------------------
# Dataset creatinfg function ####
# ------------------------------------------------------------------------------

#' Create ligerDataset object
#' @param raw.data matrix
#' @param norm.data matrix
#' @param annotation data.frame
#' @param V matrix
#' @param A matrix
#' @param B matrix
#' @param U matrix
#' @export
createLigerDataset <- function(
        raw.data = NULL,
        modal = c("default", "rna", "atac"),
        norm.data = NULL,
        scale.data = NULL,
        scale.unshared.data = NULL,
        var.unshared.features = NULL,
        H = NULL,
        V = NULL,
        A = NULL,
        B = NULL,
        U = NULL,
        h5file.info = NULL,
        feature.meta = NULL,
        ...
) {
    modal <- match.arg(modal)
    args <- as.list(environment())
    additional <- list(...)
    # TODO h5 file support
    # Necessary initialization of slots
    if (is.null(raw.data) & is.null(norm.data) & is.null(scale.data)) {
        stop("At least one type of expression data (raw.data, norm.data or ",
             "scale.data) has to be provided")
    }
    # Look for proper colnames and rownames
    cn <- NULL
    rn <- NULL
    for (i in c("raw.data", "norm.data", "scale.data")) {
        cn <- colnames(args[[i]])
        if (!is.null(cn)) break
    }
    if (!is.null(raw.data)) {
        rn <- rownames(raw.data)
        if (!inherits(raw.data, "dgCMatrix"))
            raw.data <- as(raw.data, "CsparseMatrix")
    }
    if (!is.null(norm.data)) {
        if (is.null(rn)) rn <- rownames(norm.data)
        if (!inherits(norm.data, "dgCMatrix"))
            raw.data <- as(norm.data, "CsparseMatrix")
    }
    if (!is.null(scale.data)) {
        if (is.null(rn)) rn <- rownames(scale.data)
    }
    if (is.null(h5file.info)) h5file.info <- list()
    if (is.null(feature.meta))
        feature.meta <- S4Vectors::DataFrame(row.names = rn)
    # Create ligerDataset
    allData <- list(.modalClassDict[[modal]],
                    raw.data = raw.data,
                    norm.data = norm.data,
                    scale.data = scale.data,
                    H =H, V = V, A = A, B = B, U = U,
                    h5file.info = h5file.info, feature.meta = feature.meta,
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

#' @export
createLigerDataset.h5 <- function(
        h5file,
        format.type = NULL,
        raw.data = NULL,
        norm.data = NULL,
        scale.data = NULL,
        barcodes.name = NULL,
        genes.name = NULL,
        indices.name = NULL,
        indptr.name = NULL,
        modal = c("default", "rna", "atac"),
        feature.meta = NULL,
        ...
) {
    if (!hdf5r::is_hdf5(h5file)) {
        stop("Please specify an HDF5 filename to argument `h5file`.")
    }
    modal <- match.arg(modal)
    additional <- list(...)
    h5file <- hdf5r::H5File$new(h5file, mode = "r+")
    if (!is.null(format.type) &&
        format.type %in% c("10X", "AnnData")) {
        if (format.type == "10X") {
            barcodes.name <- "matrix/barcodes"
            barcodes <- h5file[[barcodes.name]][]
            raw.data <- "matrix/data"
            indices.name <- "matrix/indices"
            indptr.name <- "matrix/indptr"
            genes.name <- "matrix/features/name"
            genes <- h5file[[genes.name]][]
        } else if (format.type == "AnnData") {
            barcodes.name <- "obs"
            barcodes <- h5file[[barcodes.name]][]$cell
            raw.data <- "raw.X/data"
            indices.name <- "raw.X/indices"
            indptr.name <- "raw.X/indptr"
            genes.name <- "raw.var"
            genes <- h5file[[genes.name]][]
        } else {
            stop("Specified `format.type` '", format.type,
                 "' is not supported for now.")
        }
    } else {
        barcodes <- h5file[[barcodes.name]][]
        genes <- h5file[[genes.name]][]
    }
    # The order of list elements matters. Put "paths" together so easier for
    # checking link existence.
    h5.meta <- list(
        H5File = h5file,
        filename = h5file$filename,
        format.type = format.type,
        indices.name = indices.name,
        indptr.name = indptr.name,
        barcodes.name = barcodes.name,
        genes.name = genes.name,
        raw.data = raw.data,
        norm.data = norm.data,
        scale.data = scale.data
    )
    if (!is.null(raw.data)) raw.data <- h5file[[raw.data]]
    if (!is.null(norm.data)) norm.data <- h5file[[norm.data]]
    if (!is.null(scale.data)) scale.data <- h5file[[scale.data]]
    if (is.null(feature.meta))
        feature.meta <- S4Vectors::DataFrame(row.names = genes)
    allData <- list(.modalClassDict[[modal]],
                    raw.data = raw.data,
                    norm.data = norm.data,
                    scale.data = scale.data,
                    #H = H, V = V, A = A, B = B, U = U,
                    h5file.info = h5.meta, feature.meta = feature.meta,
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
        if (length(object@h5file.info) == 0) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    } else if (inherits(object, "liger")) {
        if (is.null(dataset)) dataset <- names(object)
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
    for (slot in c("raw.data", "norm.data", "scale.data", "scale.unshared.data",
                   "H", "peak")) {
        if (!slot %in% methods::slotNames(x)) next
        data <- methods::slot(x, slot)
        if (!is.null(data)) {
            #if (slot == "annotation") {
            #    barcodes.slot <- rownames(data)
            #} else {
            barcodes.slot <- colnames(data)
            #}
            if (!identical(colnames(x), barcodes.slot)) {
                return(paste0("Inconsistant cell identifiers in ", slot, "."))
            }
        }
    }
    for (slot in c("scale.data", "scale.unshared.data", "V")) {
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
    if (!"H5File" %in% names(h5file.info(x))) {
        return(paste("`h5file.info` incomplete.", restoreGuide))
    }
    h5file <- getH5File(x)
    if (is.null(h5file)) {
        return(paste("`H5File` is NULL in `h5file.info` slot.", restoreGuide))
    }
    if (!h5file$is_valid) {
        return(paste("`H5File` is invalid in `h5file.info` slot.", restoreGuide))
    }
    if (!is.null(raw.data(x))) {
        if (!raw.data(x)$is_valid) {
            return(paste("`raw.data` slot is invalid.", restoreGuide))
        }
    }
    if (!is.null(norm.data(x))) {
        if (!norm.data(x)$is_valid) {
            return(paste("`norm.data` slot is invalid.", restoreGuide))
        }
    }
    if (!is.null(scale.data(x))) {
        if (!scale.data(x)$is_valid) {
            return(paste("`scale.data` slot is invalid.", restoreGuide))
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
    # TODO debating on whether to have check of the matching between scale.data
    # features and selected variable features.
}

setValidity("ligerDataset", .valid.ligerDataset)

# ------------------------------------------------------------------------------
# Generics ####
# ------------------------------------------------------------------------------

setGeneric("raw.data", function(x) standardGeneric("raw.data"))
setGeneric("raw.data<-", function(x, check = TRUE, value) standardGeneric("raw.data<-"))
setGeneric("norm.data", function(x) standardGeneric("norm.data"))
setGeneric("norm.data<-", function(x, check = TRUE, value) standardGeneric("norm.data<-"))
setGeneric("scale.data", function(x, name = NULL) standardGeneric("scale.data"))
setGeneric("scale.data<-", function(x, check = TRUE, value) standardGeneric("scale.data<-"))
setGeneric("getH5File", function(x, dataset = NULL) standardGeneric("getH5File"))
setGeneric("h5file.info", function(x, info = NULL) standardGeneric("h5file.info"))
setGeneric("h5file.info<-", function(x, info = NULL, check = TRUE, value) standardGeneric("h5file.info<-"))
setGeneric("feature.meta", function(x, check = NULL) standardGeneric("feature.meta"))
setGeneric("feature.meta<-", function(x, check = TRUE, value) standardGeneric("feature.meta<-"))
setGeneric("getMatrix", function(x, slot = "raw.data", dataset = NULL) standardGeneric("getMatrix"))
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
            !isTRUE(validObject(object, test = TRUE))) {
            warning("Link to HDF5 file fails. Please try running ",
                    "`restoreH5Liger(object)`.")
            return()
        }
        for (slot in c("raw.data", "norm.data", "scale.data")) {
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
        if ("peak" %in% methods::slotNames(object)) {
            if (!is.null(peak(object)))
                cat("peak:", nrow(peak(object)), "regions\n")
        }

        invisible(x = NULL)
    }
)

setMethod("dim", "ligerDataset", function(x) {
    nr <- length(x@rownames)
    nc <- length(x@colnames)
    c(nr, nc)
})

setMethod("dimnames", "ligerDataset", function(x) {
    rn <- x@rownames
    cn <- x@colnames
    list(rn, cn)
})

setReplaceMethod("dimnames", c("ligerDataset", "list"), function(x, value) {
    if (!isH5Liger(x)) {
        if (!is.null(raw.data(x))) {
            rownames(x@raw.data) <- value[[1L]]
            colnames(x@raw.data) <- value[[2L]]
        }
        if (!is.null(norm.data(x))) {
            rownames(x@norm.data) <- value[[1L]]
            colnames(x@norm.data) <- value[[2L]]
        }
        if (!is.null(scale.data(x))) {
            colnames(x@scale.data) <- value[[2L]]
            rownames(x@scale.data) <-
                value[[1L]][x@rownames %in% rownames(x@scale.data)]
        }
        if (!is.null(x@scale.unshared.data)) {
            colnames(x@scale.unshared.data) <- value[[2L]]
            rownames(x@scale.unshared.data) <-
                value[[1L]][x@rownames %in% rownames(x@scale.unshared.data)]
        }
    }
    # else {
        #h5file <- getH5File(x)
        #h5meta <- h5file.info(x)
        #h5file[[h5meta$genes.name]][1:length(value[[1L]])] <- value[[1L]]
        #h5file[[h5meta$barcodes.name]][1:length(value[[2L]])] <- value[[2L]]
    #}
    if (!is.null(x@H))
        rownames(x@H) <- value[[2L]]
    if (!is.null(x@V))
        colnames(x@V) <- value[[1L]][x@rownames %in% rownames(x@V)]
    if (!is.null(x@B))
        rownames(x@B) <- value[[1L]][x@rownames %in% rownames(x@V)]
    if (!is.null(x@U)) colnames(x@U) <- value[[2L]]
    if ("peak" %in% methods::slotNames(x)) {
        if (!is.null(peak(x))) colnames(peak(x)) <- value[[2L]]
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

#' Access data in a ligerDataset object
#' @param x \linkS4class{ligerDataset} object
#' @param value matrix like or \code{NULL}
#' @param check Whether to perform object validity check on setting new value.
#' @export
#' @rdname data-access
setMethod("raw.data", "ligerDataset",
          function(x) x@raw.data)

#' @export
#' @rdname data-access
setReplaceMethod(
    "raw.data",
    signature(x = "ligerDataset", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, check = TRUE, value) {
        if (isH5Liger(x))
            stop("Cannot replace slot with in-memory data for H5 based object.")
        x@raw.data <- value
        if (isTRUE(check)) validObject(x)
        x
    }
)

#' @export
#' @rdname data-access
setReplaceMethod(
    "raw.data",
    signature(x = "ligerDataset", check = "ANY", value = "H5D"),
    function(x, check = TRUE, value) {
        if (!isH5Liger(x))
            stop("Cannot replace slot with on-disk data for in-memory object.")
        x@raw.data <- value
        if (isTRUE(check)) validObject(x)
        x
    }
)

#' @export
#' @rdname data-access
setMethod("norm.data", "ligerDataset",
          function(x) x@norm.data)

#' @export
#' @rdname data-access
setReplaceMethod(
    "norm.data",
    signature(x = "ligerDataset", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, check = TRUE, value) {
        if (isH5Liger(x))
            stop("Cannot replace slot with in-memory data for H5 based object.")
        x@norm.data <- value
        if (isTRUE(check)) validObject(x)
        x
    }
)

#' @export
#' @rdname data-access
setReplaceMethod(
    "norm.data",
    signature(x = "ligerDataset", check = "ANY", value = "H5D"),
    function(x, check = TRUE, value) {
        if (!isH5Liger(x))
            stop("Cannot replace slot with on-disk data for in-memory object.")
        x@norm.data <- value
        if (isTRUE(check)) validObject(x)
        x
    }
)

#' @export
#' @rdname data-access
setMethod("scale.data", "ligerDataset",
          function(x, name = NULL) x@scale.data)

#' @export
#' @rdname data-access
setMethod(
    "scale.data",
    signature(x = "liger", name = "character"),
    function(x, name) {
        scale.data(dataset(x, name))
    }
)

#' @export
#' @rdname data-access
setMethod(
    "scale.data",
    signature(x = "liger", name = "numeric"),
    function(x, name) {
        scale.data(dataset(x, name))
    }
)

#' @export
#' @rdname data-access
setReplaceMethod(
    "scale.data",
    signature(x = "ligerDataset", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, check = TRUE, value) {
        if (isH5Liger(x))
            stop("Cannot replace slot with in-memory data for H5 based object.")
        x@scale.data <- value
        if (isTRUE(check)) validObject(x)
        x
    }
)

#' @export
#' @rdname data-access
setReplaceMethod(
    "scale.data",
    signature(x = "ligerDataset", check = "ANY", value = "H5D"),
    function(x, check = TRUE, value) {
        if (!isH5Liger(x))
            stop("Cannot replace slot with on-disk data for in-memory object.")
        x@scale.data <- value
        if (isTRUE(check)) validObject(x)
        x
    }
)

#' Access the H5File object
#' @param x \linkS4class{ligerDataset} object or \linkS4class{liger} object.
#' @param dataset Get H5File from specific dataset(s) if using a
#' \linkS4class{liger} object.
#' @export
#' @rdname getH5File
setMethod("getH5File",
          signature = signature(x = "ligerDataset", dataset = "missing"),
          function(x, dataset = NULL) h5file.info(x, "H5File"))

#' @rdname getH5File
#' @export
setMethod("getH5File",
          signature = signature(x = "liger", dataset = "character"),
          function(x, dataset = NULL) {
              if (is.null(dataset)) dataset <- names(x)
              if (any(!dataset %in% names(x))) {
                  stop("Specified dataset name(s) not found: ",
                       paste(dataset[!dataset %in% names(x)], collapse = ", "))
              }
              results <- lapply(datasets(x)[dataset],
                                function(ld) h5file.info(ld, "H5File"))
              if (length(results) == 1) results <- results[[1]]
              results
          })

#' Access the H5File information
#' @description H5File information generally stores the "paths" to all types of
#' data in the H5 file. Setter method changes not only information stored in
#' \code{h5file.info} slot, but also data slots (i.e. \code{raw.data},
#' \code{norm.data} or \code{scale.data}) that a related with the corresponding
#' information.
#' @param x \linkS4class{ligerDataset}
#' @param info Character vector of queried information. Default \code{NULL}
#' returns a full list of information.
#' @export
#' @rdname h5file.info
setMethod(
    "h5file.info",
    signature = signature(x = "ligerDataset", info = "ANY"),
    function(x, info = NULL) {
        if (is.null(info)) result <- x@h5file.info
        else {
            if (length(info) == 1) result <- x@h5file.info[[info]]
            else {
                if (any(!info %in% names(x@h5file.info))) {
                    stop("Specified h5file info not found: ",
                         paste(info[!info %in% names(x@h5file.info)],
                               collapse = ", "))
                }
                result <- x@h5file.info[info]
                names(result) <- info
            }
        }
        return(result)
    })

#' @export
#' @rdname h5file.info
setReplaceMethod(
    "h5file.info",
    signature = signature(
        x = "ligerDataset",
        info = "ANY",
        check = "ANY",
        value = "ANY"
    ),
    function(x, info = NULL, check = TRUE, value) {
        if (is.null(info)) {
            x@h5file.info <- value
        } else {
            if (!is.character(info) | length(info) != 1)
                stop('`info` has to be a single character.')
            if (info %in% c('indices.name', 'indptr.name', 'barcodes.name',
                            'genes.name', 'raw.data', 'norm.data',
                            'scale.data')) {
                if (!getH5File(x)$exists(value)) {
                    stop("Specified info is invalid, '", info,
                         "' does not exists in the HDF5 file.")
                }
            }
            x@h5file.info[[info]] <- value
            if (info %in% c('raw.data', 'norm.data', 'scale.data')) {
                x <- do.call(paste0(info, "<-"),
                             list(x = x,
                                  value = getH5File(x)[[h5file.info(x, info)]],
                                  check = check))
            }
        }
        if (isTRUE(check)) validObject(x)
        x
    }
)

#' Access the feature metadata in a ligerDataset object
#' @param x \linkS4class{ligerDataset} object
#' @rdname feature.meta
#' @export
setMethod("feature.meta", signature(x = "ligerDataset", check = "ANY"),
          function(x, check = NULL) {
    x@feature.meta
})

#' @rdname feature.meta
#' @export
setReplaceMethod(
    "feature.meta",
    signature(x = "ligerDataset", check = "ANY"),
    function(x, check = TRUE, value) {
        if (!inherits(value, "DFrame"))
            value <- S4Vectors::DataFrame(value)
        x@feature.meta <- value
        if (isTRUE(check)) validObject(x)
        x
    }
)

#' Get matrix from liger/ligerDataset object
#' @param x \linkS4class{liger} or \linkS4class{ligerDataset} object
#' @param slot The type of matrix to retrieve. Choose from \code{"raw.data"},
#' \code{"norm.data"}, \code{"scale.data"}, \code{"H"}, \code{"V"} for both
#' classes. For \linkS4class{liger} object, \code{"W"} can also be an option,
#' ignoring \code{dataset}.
#' @param dataset When using a \linkS4class{liger} object, from which dataset(s)
#' to retrieve the data. Single dataset returns the matrix, multiple datasets
#' return a named list of the type of matrices, default \code{NULL} return a
#' named list of the type of matrix from all datasets.
#' @return A matrix or a list, depending on \code{dataset}
#' @rdname getMatrix
#' @export
setMethod("getMatrix", signature(x = "ligerDataset", dataset = "missing"),
          function(x, slot = c("raw.data", "norm.data", "scale.data", "H", "V"),
                   dataset = NULL) {
              # TODO: Currently directly find the data with slot, but need to
              # think about maintainability when we need to change slot name.
              slot <- match.arg(slot)
              methods::slot(x, slot)
          })

#' @rdname getMatrix
#' @export
setMethod("getMatrix", signature(x = "liger"),
          function(x, slot = c("raw.data", "norm.data", "scale.data", "H", "V", "W", "H.norm"),
                   dataset = NULL) {
              slot <- match.arg(slot)
              if (slot == "W") return(x@W)
              if (slot == "H.norm") return(x@H.norm)
              if (is.null(dataset)) {
                  return(lapply(datasets(x), function(ld) getMatrix(ld, slot)))
              } else {
                  if (length(dataset) == 1)
                      return(getMatrix(dataset(x, dataset), slot))
                  else {
                      lds <- datasets(x)[dataset]
                      return(lapply(lds, function(ld) getMatrix(ld, slot)))
                  }
              }
          })

