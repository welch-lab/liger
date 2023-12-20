
################################################################################
# Developer guide for adding a new sub-class of `ligerDataset` for new modality
################################################################################
#
# Below is a check-list of the TODOs when new sub-classes need to be added.
# Please follow them carefully, and refer to existing code as examples.
#
# 1. Add `setClass` chunk for defining the new subclass, pay attention to:
#   a. Naming convention should be `liger{Modal}Dataset`, in camelCase
#   b. contains = "ligerDataset"
#   c. add new slots for modality specific information with `representation`
#   d. if the default new information could be empty, add `prototype`
# 2. In files `zzz.R` and `import.R`, search for
#    text "modal". When seeing a multi-option vector argument, add a unique
#    abbreviation of this new data type to the vector.
# 3. In file `ligerDataset-class.R`, find the list object `.modalClassDict`,
#    add the new entry for your sub-class. The name is the abbr you add in step
#    2., and the value should be the name you set for the new sub-class in step
#    1. Then the next list object `.classModalDict` should also be taken care
#    of in a similar way.
# 4. If the new slot(s) added is thought to be retrieved by future developers
#    or users, getter and setter methods MUST be implemented.
# 5. Please go through the implementation of the following functions in file
#    `ligerDataset-class.R`, and make sure data in the new slot(s) is properly
#    handled.
#   a. .checkLigerDatasetBarcodes()
#   b. `dimnames<-()` (search: `setReplaceMethod("dimnames"`)
#   c. `[` (search: "[")
#
################################################################################

setClassUnion("matrixLike_OR_NULL", c(
    "matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix", "NULL"
))

#' Subclass of ligerDataset for RNA modality
#' @description Inherits from \linkS4class{ligerDataset} class. Contained slots
#' can be referred with the link. This subclass does not have any different from
#' the default \code{ligerDataset} class except the class name.
#' @export
#' @exportClass ligerRNADataset
ligerRNADataset <- setClass(
    "ligerRNADataset", contains = "ligerDataset"
)

#-------------------------------------------------------------------------------
# Sub-class for ATAC data ####
#-------------------------------------------------------------------------------

#' Subclass of ligerDataset for ATAC modality
#'
#' @description Inherits from \linkS4class{ligerDataset} class. Contained slots
#' can be referred with the link.
#' @slot rawPeak sparse matrix
#' @slot normPeak sparse matrix
#' @exportClass ligerATACDataset
#' @export
ligerATACDataset <- setClass(
    "ligerATACDataset",
    contains = "ligerDataset",
    representation = representation(rawPeak = "matrixLike_OR_NULL",
                                    normPeak = "matrixLike_OR_NULL"),
    prototype = prototype(rawPeak = NULL, normPeak = NULL)
)

#' Access ligerATACDataset peak data
#' @description Similar as how default \linkS4class{ligerDataset} data is
#' accessed.
#' @param x \linkS4class{ligerATACDataset} object or a \linkS4class{liger}
#' object.
#' @param dataset Name or numeric index of an ATAC dataset.
#' @param check Logical, whether to perform object validity check on setting new
#' value.
#' @param value \code{\link[Matrix]{dgCMatrix-class}} matrix.
#' @return The retrieved peak count matrix or the updated \code{x} object.
#' @rdname peak
#' @export
setGeneric("rawPeak", function(x, dataset) standardGeneric("rawPeak"))

#' @rdname peak
#' @export
setGeneric("rawPeak<-", function(x, dataset, check = TRUE, value) standardGeneric("rawPeak<-"))

#' @rdname peak
#' @export
setGeneric("normPeak", function(x, dataset) standardGeneric("normPeak"))

#' @rdname peak
#' @export
setGeneric("normPeak<-", function(x, dataset, check = TRUE, value) standardGeneric("normPeak<-"))

#' @rdname peak
#' @export
setMethod("rawPeak", signature(x = "liger", dataset = "character"),
          function(x, dataset) {
              atac <- dataset(x, dataset)
              if (!inherits(atac, "ligerATACDataset")) {
                  stop("Specified dataset is not of ligerATACDataset class.")
              }
              atac@rawPeak
          })

#' @rdname peak
#' @export
setReplaceMethod(
    "rawPeak",
    signature(x = "liger", dataset = "character"),
    function(x, dataset, check = TRUE, value) {
        if (!inherits(dataset(x, dataset), "ligerATACDataset"))
            stop("Specified dataset is not of `ligerATACDataset` class.")
        x@datasets[[dataset]]@rawPeak <- value
        if (isTRUE(check)) methods::validObject(dataset(x, dataset))
        x
    })

#' @rdname peak
#' @export
setMethod("rawPeak", signature(x = "ligerATACDataset", dataset = "missing"),
          function(x, dataset = NULL) {
              x@rawPeak
          })

#' @rdname peak
#' @export
setReplaceMethod(
    "rawPeak",
    signature(x = "ligerATACDataset", dataset = "missing"),
    function(x, dataset = NULL, check = TRUE, value) {
        x@rawPeak <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    })

#' @rdname peak
#' @export
setMethod("normPeak", signature(x = "liger", dataset = "character"),
          function(x, dataset) {
              atac <- dataset(x, dataset)
              if (!inherits(atac, "ligerATACDataset")) {
                  stop("Specified dataset is not of ligerATACDataset class.")
              }
              atac@normPeak
          })

#' @rdname peak
#' @export
setReplaceMethod(
    "normPeak",
    signature(x = "liger", dataset = "character"),
    function(x, dataset, check = TRUE, value) {
        if (!inherits(dataset(x, dataset), "ligerATACDataset"))
            stop("Specified dataset is not of `ligerATACDataset` class.")
        x@datasets[[dataset]]@normPeak <- value
        if (isTRUE(check)) methods::validObject(dataset(x, dataset))
        x
    })

#' @rdname peak
#' @export
setMethod("normPeak", signature(x = "ligerATACDataset", dataset = "missing"),
          function(x, dataset = NULL) {
              x@normPeak
          })

#' @rdname peak
#' @export
setReplaceMethod(
    "normPeak",
    signature(x = "ligerATACDataset", dataset = "missing"),
    function(x, dataset = NULL, check = TRUE, value) {
        x@normPeak <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    })

.valid.ligerATACDataset <- function(object) {
    passSuperClassCheck <- .valid.ligerDataset(object)
    if (!isTRUE(passSuperClassCheck)) return(passSuperClassCheck)
    for (slot in c("rawPeak", "normPeak")) {
        data <- methods::slot(object, slot)
        if (!is.null(data)) {
            barcodes.slot <- colnames(data)
            if (!identical(object@colnames, barcodes.slot)) {
                return(paste0("Inconsistant cell identifiers in `", slot,
                              "` slot."))
            }
        }
    }
}

setValidity("ligerATACDataset", .valid.ligerATACDataset)

#-------------------------------------------------------------------------------
# Sub-class for Spatial data ####
#-------------------------------------------------------------------------------

#' Subclass of ligerDataset for Spatial modality
#'
#' @description Inherits from \linkS4class{ligerDataset} class. Contained slots
#' can be referred with the link.
#' @slot coordinate dense matrix
#' @exportClass ligerSpatialDataset
#' @export
ligerSpatialDataset <- setClass(
    "ligerSpatialDataset",
    contains = "ligerDataset",
    representation = representation(coordinate = "matrix_OR_NULL"),
    prototype = prototype(coordinate = NULL)
)

#' Access ligerSpatialDataset coordinate data
#' @description Similar as how default \linkS4class{ligerDataset} data is
#' accessed.
#' @param x \linkS4class{ligerSpatialDataset} object or a \linkS4class{liger}
#' object.
#' @param dataset Name or numeric index of an spatial dataset.
#' @param check Logical, whether to perform object validity check on setting new
#' value.
#' @param value \code{\link{matrix}}.
#' @return The retrieved coordinate matrix or the updated \code{x} object.
#' @rdname coordinate
#' @export
setGeneric("coordinate", function(x, dataset) standardGeneric("coordinate"))

#' @rdname coordinate
#' @export
setGeneric("coordinate<-", function(x, dataset, check = TRUE, value) standardGeneric("coordinate<-"))


#' @rdname coordinate
#' @export
setMethod("coordinate", signature(x = "liger", dataset = "character"),
          function(x, dataset) {
              spatial <- dataset(x, dataset)
              if (!inherits(spatial, "ligerSpatialDataset")) {
                  stop("Specified dataset is not of `ligerSpatialDataset` class.")
              }
              spatial@coordinate
          })

#' @rdname coordinate
#' @export
setReplaceMethod(
    "coordinate",
    signature(x = "liger", dataset = "character"),
    function(x, dataset, check = TRUE, value) {
        if (!inherits(dataset(x, dataset), "ligerSpatialDataset"))
            stop("Specified dataset is not of `ligerSpatialDataset` class.")
        value <- .checkCoords(ld = dataset(x, dataset), value = value)
        x@datasets[[dataset]]@coordinate <- value
        if (isTRUE(check)) methods::validObject(dataset(x, dataset))
        x
    })

#' @rdname coordinate
#' @export
setMethod("coordinate", signature(x = "ligerSpatialDataset", dataset = "missing"),
          function(x, dataset = NULL) {
              x@coordinate
          })

#' @rdname coordinate
#' @export
setReplaceMethod(
    "coordinate",
    signature(x = "ligerSpatialDataset", dataset = "missing"),
    function(x, dataset = NULL, check = TRUE, value) {
        value <- .checkCoords(ld = x, value = value)
        x@coordinate <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    })

.checkCoords <- function(ld, value) {
    if (is.null(rownames(value))) {
        warning("No rownames with given spatial coordinate, ",
                "assuming they match with the cells.")
        rownames(value) <- colnames(ld)
    }
    if (is.null(colnames(value))) {
        if (ncol(value) <= 3) {
            colnames(value) <- c("x", "y", "z")[seq(ncol(value))]
        } else {
            stop("More than 3 dimensions for the coordinates but no ",
                 "colnames are given.")
        }
        warning("No colnames with given spatial coordinate, ",
                "setting to ", paste0(colnames(value), collapse = ", "))
    }
    if (nrow(value) != ncol(ld)) {
        full <- matrix(NA, nrow = ncol(ld), ncol = ncol(value),
                       dimnames = list(colnames(ld), colnames(value)))
        full[rownames(value), colnames(value)] <- value
        value <- full
        warning("NA generated for missing cells.")
    }
    return(value)
}

.valid.ligerSpatialDataset <- function(object) {
    passSuperClassCheck <- .valid.ligerDataset(object)
    if (!isTRUE(passSuperClassCheck)) return(passSuperClassCheck)
    coord <- object@coordinate
    if (!is.null(coord)) {
        barcodes.slot <- rownames(coord)
        if (!identical(object@colnames, barcodes.slot)) {
            return(paste0("Inconsistant cell identifiers in `coordinate` slot."))
        }
    }
}

setValidity("ligerSpatialDataset", .valid.ligerSpatialDataset)

