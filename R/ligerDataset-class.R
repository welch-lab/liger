setClassUnion("dgCMatrix_OR_NULL", c("dgCMatrix", "NULL"))
setClassUnion("matrix_OR_NULL", c("matrix", "NULL"))
setClassUnion("matrixLike_OR_NULL", c(
    "matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix", "NULL"
))

setClassUnion("dgCMatrix_OR_H5D_OR_NULL", c("dgCMatrix_OR_NULL", "H5D"))
setClassUnion("matrix_OR_H5D_OR_NULL", c("matrix_OR_NULL", "H5D"))

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
#' @slot W matrix
#' @slot V matrix
#' @slot A matrix
#' @slot B matrix
#' @slot U matrix
#' @slot agg.data list
#' @slot h5file.info list
#' @exportClass ligerDataset
ligerDataset <- setClass(
    "ligerDataset",
    representation(
        raw.data = "dgCMatrix_OR_H5D_OR_NULL",
        norm.data = "dgCMatrix_OR_H5D_OR_NULL",
        scale.data = "matrix_OR_H5D_OR_NULL",
        scale.unshared.data = "matrix_OR_H5D_OR_NULL",
        var.unshared.features = "character",
        W = "matrix_OR_NULL",
        V = "matrix_OR_NULL",
        A = "matrix_OR_NULL",
        B = "matrix_OR_NULL",
        U = "matrix_OR_NULL",
        agg.data = "list",
        h5file.info = "list",
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
#' @param W matrix
#' @param V matrix
#' @param A matrix
#' @param B matrix
#' @param U matrix
#' @export
ligerDataset <- function(raw.data = NULL,
                         modal = c("default", "rna", "atac"),
                         norm.data = NULL,
                         scale.data = NULL,
                         scale.unshared.data = NULL,
                         var.unshared.features = NULL,
                         W = NULL,
                         V = NULL,
                         A = NULL,
                         B = NULL,
                         U = NULL,
                         h5file.info = NULL,
                         ...) {
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
    if (!is.null(raw.data) & !inherits(raw.data, "dgCMatrix")) {
        raw.data <- as(raw.data, "CsparseMatrix")
        rn <- rownames(raw.data)
    }
    if (!is.null(norm.data) & !inherits(norm.data, "dgCMatrix")) {
        norm.data <- as(norm.data, "CsparseMatrix")
    }
    if (is.null(h5file.info)) h5file.info <- list()
    # Create ligerDataset
    allData <- list(.modalClassDict[[modal]],
                    raw.data = raw.data,
                    norm.data = norm.data,
                    scale.data = scale.data,
                    W = W, V = V, A = A, B = B, U = U,
                    h5file.info = h5file.info,
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

.ligerDataset.h5 <- function(
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
    allData <- list(.modalClassDict[[modal]],
                    raw.data = raw.data,
                    norm.data = norm.data,
                    scale.data = scale.data,
                    #W = W, V = V, A = A, B = B, U = U,
                    h5file.info = h5.meta,
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
isOnlineLiger <- function(object, dataset = NULL) {
    if (inherits(object, "ligerDataset")) {
        if (length(object@h5file.info) == 0) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    } else if (inherits(object, "liger")) {
        if (!is.null(dataset)) {
            return(isOnlineLiger(dataset(object, dataset)))
        } else {
            allCheck <- unlist(lapply(datasets(object), isOnlineLiger))
            return(all(allCheck))
        }
    } else {
        warning("Given object is not of liger or ligerDataset class.")
        return(FALSE)
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
                   "W", "V", "A", "B", "U", "peak")) {
        if (!slot %in% methods::slotNames(x)) next
        data <- methods::slot(x, slot)
        if (!is.null(data)) {
            #if (slot == "annotation") {
            #    barcodes.slot <- rownames(data)
            #} else {
            barcodes.slot <- colnames(data)
            #}
            if (!identical(colnames(x), barcodes.slot)) {
                return(paste0("Inconsistant cell identifiers."))
            }
        }
    }
    TRUE
}

.checkOnlineLigerDatasetLink <- function(x) {
    restoreGuide <- "Please try running `restoreOnlineLiger(object)`."
    if (!"H5File" %in% names(x@h5file.info)) {
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
    if (isOnlineLiger(object)) {
        message("Pretend to have checked h5 validity")
        .checkOnlineLigerDatasetLink(object)
    } else {
        .checkLigerDatasetBarcodes(object)
    }
    # TODO more checks
}

setValidity("ligerDataset", .valid.ligerDataset)

# ------------------------------------------------------------------------------
# Generics ####
# ------------------------------------------------------------------------------

setGeneric("raw.data", function(x, ...) standardGeneric("raw.data"))
setGeneric("raw.data<-", function(x, check = TRUE, ...) standardGeneric("raw.data<-"))
setGeneric("norm.data", function(x, ...) standardGeneric("norm.data"))
setGeneric("norm.data<-", function(x, check = TRUE, ...) standardGeneric("norm.data<-"))
setGeneric("scale.data", function(x, ...) standardGeneric("scale.data"))
setGeneric("scale.data<-", function(x, check = TRUE, ...) standardGeneric("scale.data<-"))
setGeneric("getH5File", function(x, dataset = NULL) standardGeneric("getH5File"))

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
        if (isOnlineLiger(object) &
            !isTRUE(validObject(object, test = TRUE))) {
            warning("Link to HDF5 file fails. Please try running ",
                    "`restoreOnlineLiger(object)`.")
            return()
        }
        for (slot in c("raw.data", "norm.data", "scale.data")) {
            data <- methods::slot(object, slot)
            if (!is.null(data)) {
                if (inherits(data, "dgCMatrix")) {
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
    rownames(x@raw.data) <- value[[1L]]
    x@colnames <- value[[2L]]
    if (!is.null(raw.data(x))) colnames(raw.data(x)) <- value[[2L]]
    if (!is.null(norm.data(x))) colnames(norm.data(x)) <- value[[2L]]
    if (!is.null(scale.data(x))) colnames(scale.data(x)) <- value[[2L]]
    if (!is.null(x@scale.unshared.data))
        colnames(x@scale.unshared.data) <- value[[2L]]
    if (!is.null(x@W)) colnames(x@W) <- value[[2L]]
    if (!is.null(x@V)) colnames(x@V) <- value[[2L]]
    if (!is.null(x@A)) colnames(x@A) <- value[[2L]]
    if (!is.null(x@B)) colnames(x@B) <- value[[2L]]
    if (!is.null(x@U)) colnames(x@U) <- value[[2L]]
    if ("peak" %in% methods::slotNames(x)) {
        if (!is.null(peak(x))) colnames(peak(x)) <- value[[2L]]
    }
    return(x)
})

setMethod(
    "[",
    signature = "ligerDataset",
    function(x, i, ...) {
        modal <- .classModalDict[[class(x)]]
        subsetData <- list(
            modal = modal,
            raw.data = raw.data(x)[, i, drop = FALSE],
            norm.data = norm.data(x)[, i, drop = FALSE],
            scale.data = scale.data(x)[, i, drop = FALSE],
            scale.unshared.data = x@scale.unshared.data[, i, drop = FALSE],
            W = x@W[, i, drop = FALSE],
            V = x@V[, i, drop = FALSE],
            A = x@A[, i, drop = FALSE],
            B = x@B[, i, drop = FALSE],
            U = x@U[, i, drop = FALSE]
        )
        # Additional subsetting for sub-classes, if applicable
        if (modal == "atac") {
            subsetData$peak <- peak(x)[, i, drop = FALSE]
        }

        x <- do.call("ligerDataset", subsetData)
        validObject(x)
        x
    }
)

#' Access raw.data of ligerDataset object
#' @param x ligerDataset object
#' @export
setMethod("raw.data", "ligerDataset",
          function(x) x@raw.data)

#' Set raw.data to ligerDataset object
#' @param x ligerDataset object
#' @param value dgCMatrix
#' @export
setReplaceMethod("raw.data", "ligerDataset",
                 function(x, check = TRUE, value) {
                     x@raw.data <- value
                     if (isTRUE(check)) validObject(x)
                     x
                 })

#' Access norm.data of ligerDataset object
#' @param x ligerDataset object
#' @export
setMethod("norm.data", "ligerDataset",
          function(x) x@norm.data)

#' Set norm.data to ligerDataset object
#' @param x ligerDataset object
#' @param value dgCMatrix
#' @export
setReplaceMethod("norm.data", "ligerDataset",
                 function(x, check = TRUE, value) {
                     x@norm.data <- value
                     if (isTRUE(check)) validObject(x)
                     x
                 })

#' Access norm.data of ligerDataset object
#' @param x ligerDataset object
#' @export
setMethod("scale.data", "ligerDataset",
          function(x) x@scale.data)

#' Set norm.data to ligerDataset object
#' @param x ligerDataset object
#' @param value dgCMatrix
#' @export
setReplaceMethod("scale.data", "ligerDataset",
                 function(x, check = TRUE, value) {
                     x@scale.data <- value
                     if (isTRUE(check)) validObject(x)
                     x
                 })

#' Access the H5File object in a ligerDataset object
#' @param x ligerDataset object
#' @param dataset missing
#' @export
setMethod("getH5File",
          signature = signature(x = "ligerDataset", dataset = "missing"),
          function(x, dataset = NULL) x@h5file.info$H5File)

#' Access the H5File object in a specific dataset or all datasets of a liger
#' object
#' @param x liger object
#' @param dataset One or more names of datasets that can be found in \code{x}
#' @export
setMethod("getH5File",
          signature = signature(x = "liger", dataset = "character"),
          function(x, dataset = NULL) {
              if (is.null(dataset)) dataset <- names(x)
              if (any(!dataset %in% names(x))) {
                  stop("Specified dataset name(s) not found: ",
                       paste(dataset[!dataset %in% names(x)], collapse = ", "))
              }
              lapply(datasets(x)[dataset], function(ld) ld@h5file.info$H5File)
          })

