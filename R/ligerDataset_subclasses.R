
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
# 2. In files `liger-class.R`, `ligerDataset-class.R` and `IO.R`, search for
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

setGeneric("rawPeak", function(x, dataset) standardGeneric("rawPeak"))
setGeneric("rawPeak<-", function(x, dataset, check = TRUE, value) standardGeneric("rawPeak<-"))
setGeneric("normPeak", function(x, dataset) standardGeneric("normPeak"))
setGeneric("normPeak<-", function(x, dataset, check = TRUE, value) standardGeneric("normPeak<-"))

#' @export
setMethod("rawPeak", signature(x = "liger", dataset = "character"),
          function(x, dataset) {
              atac <- dataset(x, dataset)
              if (!inherits(atac, "ligerATACDataset")) {
                  stop("Specified dataset is not of ligerATACDataset class.")
              }
              atac@rawPeak
          })

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

#' @export
setMethod("rawPeak", signature(x = "ligerATACDataset", dataset = "missing"),
          function(x, dataset = NULL) {
              x@rawPeak
          })

#' @export
setReplaceMethod(
    "rawPeak",
    signature(x = "ligerATACDataset", dataset = "missing"),
    function(x, dataset = NULL, check = TRUE, value) {
        x@rawPeak <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    })

#' @export
setMethod("normPeak", signature(x = "liger", dataset = "character"),
          function(x, dataset) {
              atac <- dataset(x, dataset)
              if (!inherits(atac, "ligerATACDataset")) {
                  stop("Specified dataset is not of ligerATACDataset class.")
              }
              atac@normPeak
          })

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

#' @export
setMethod("normPeak", signature(x = "ligerATACDataset", dataset = "missing"),
          function(x, dataset = NULL) {
              x@normPeak
          })

#' @export
setReplaceMethod(
    "normPeak",
    signature(x = "ligerATACDataset", dataset = "missing"),
    function(x, dataset = NULL, check = TRUE, value) {
        x@normPeak <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    })

