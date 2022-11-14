setClassUnion("dgCMatrix_OR_NULL", c("dgCMatrix", "NULL"))
setClassUnion("matrix_OR_NULL", c("matrix", "NULL"))
setClassUnion("matrixLike_OR_NULL", c(
    "matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix", "NULL"
))

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
        raw.data = "dgCMatrix_OR_NULL",
        norm.data = "dgCMatrix_OR_NULL",
        scale.data = "matrix_OR_NULL",
        scale.unshared.data = "matrix_OR_NULL",
        var.unshared.features = "character",
        W = "matrix_OR_NULL",
        V = "matrix_OR_NULL",
        A = "matrix_OR_NULL",
        B = "matrix_OR_NULL",
        U = "matrix_OR_NULL",
        agg.data = "list",
        h5file.info = "list",
        colnames = "character"
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
                         norm.data = NULL,
                         scale.data = NULL,
                         scale.unshared.data = NULL,
                         var.unshared.features = NULL,
                         W = NULL,
                         V = NULL,
                         A = NULL,
                         B = NULL,
                         U = NULL,
                         h5file.info = NULL) {
    args <- as.list(environment())
    # TODO h5 file support
    # Necessary initialization of slots
    if (is.null(raw.data) & is.null(norm.data) & is.null(scale.data)) {
        stop("At least one type of expression data (raw.data, norm.data or ",
             "scale.data) has to be provided")
    } else {
        cn <- NULL
        for (i in seq_along(args)) {
            cn <- colnames(args[[i]])
            if (!is.null(cn)) break
        }
    }
    if (!is.null(raw.data) & !inherits(raw.data, "dgCMatrix")) {
        raw.data <- as(raw.data, "CsparseMatrix")
    }
    if (!is.null(norm.data) & !inherits(norm.data, "dgCMatrix")) {
        norm.data <- as(norm.data, "CsparseMatrix")
    }
    if (is.null(h5file.info)) h5file.info <- list()
    # Create ligerDataset
    x <- methods::new("ligerDataset",
                      raw.data = raw.data,
                      norm.data = norm.data,
                      scale.data = scale.data,
                      W = W, V = V, A = A, B = B, U = U,
                      h5file.info = h5file.info,
                      colnames = cn)
    return(x)
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
                   "W", "V", "A", "B", "U")) {
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

.valid.ligerDataset <- function(object) {
    .checkLigerDatasetBarcodes(object)
    # TODO more checks
}

setValidity("ligerDataset", .valid.ligerDataset)

# ------------------------------------------------------------------------------
# Generics ####
# ------------------------------------------------------------------------------

setGeneric("raw.data", function(x, ...) standardGeneric("raw.data"))
setGeneric("raw.data<-", function(x, ...) standardGeneric("raw.data<-"))
setGeneric("norm.data", function(x, ...) standardGeneric("norm.data"))
setGeneric("norm.data<-", function(x, ...) standardGeneric("norm.data<-"))
setGeneric("scale.data", function(x, ...) standardGeneric("scale.data"))
setGeneric("scale.data<-", function(x, ...) standardGeneric("scale.data<-"))

# ------------------------------------------------------------------------------
# Methods ####
# ------------------------------------------------------------------------------

setMethod(
    f = "show",
    signature(object = "ligerDataset"),
    definition = function(object) {
        cat("An object of class ligerDataset with", ncol(object), "cells\n")
        for (slot in c("raw.data", "norm.data", "scale.data")) {
            data <- methods::slot(object, slot)
            if (!is.null(data)) {
                cat(paste0(slot, ":"), nrow(data), "features\n")
            }
        }
        invisible(x = NULL)
    }
)

setMethod("dim", "ligerDataset", function(x) {
    nr <- ifelse(is.null(x@raw.data), NA, nrow(x@raw.data))
    nc <- length(x@colnames)
    c(nr, nc)
})

setMethod("dimnames", "ligerDataset", function(x) {
    cn <- x@colnames
    rn <- rownames(raw.data(x))
    return(list(rn, cn))
})

setReplaceMethod("dimnames", c("ligerDataset", "list"), function(x, value) {
    rownames(x@raw.data) <- value[[1L]]
    x@colnames <- value[[2L]]
    if (!is.null(x@raw.data)) colnames(x@raw.data) <- value[[2L]]
    if (!is.null(x@norm.data)) colnames(x@norm.data) <- value[[2L]]
    if (!is.null(x@scale.data)) colnames(x@scale.data) <- value[[2L]]
    if (!is.null(x@scale.unshared.data))
        colnames(x@scale.unshared.data) <- value[[2L]]
    if (!is.null(x@W)) colnames(x@W) <- value[[2L]]
    if (!is.null(x@V)) colnames(x@V) <- value[[2L]]
    if (!is.null(x@A)) colnames(x@A) <- value[[2L]]
    if (!is.null(x@B)) colnames(x@B) <- value[[2L]]
    if (!is.null(x@U)) colnames(x@U) <- value[[2L]]
    return(x)
})

setMethod(
    "[",
    signature = "ligerDataset",
    function(x, i, ...) {
        x <- ligerDataset(
            raw.data = raw.data(x)[, i, drop = FALSE],
            norm.data = norm.data(x)[, i, drop = FALSE],
            scale.data = scale.data(x)[, i, drop = FALSE],
            #annotation = x@annotation[i, , drop = FALSE],
            scale.unshared.data = x@scale.unshared.data[, i, drop = FALSE],
            W = x@W[, i, drop = FALSE],
            V = x@V[, i, drop = FALSE],
            A = x@A[, i, drop = FALSE],
            B = x@B[, i, drop = FALSE],
            U = x@U[, i, drop = FALSE]
        )
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
                 function(x, value) {
                     x@raw.data <- value
                     validObject(x)
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
                 function(x, value) {
                     x@norm.data <- value
                     validObject(x)
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
                 function(x, value) {
                     x@scale.data <- value
                     validObject(x)
                     x
                 })
