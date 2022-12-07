
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

#' @exportClass ligerATACDataset
ligerATACDataset <- setClass(
    "ligerATACDataset",
    contains = "ligerDataset",
    representation = representation(peak = "matrixLike_OR_NULL"),
    prototype = prototype(peak = NULL)
)

#' @export
setGeneric("peak", function(x, dataset) standardGeneric("peak"))
setGeneric("peak<-", function(x, dataset, value) standardGeneric("peak<-"))

#' @export
setMethod("peak", signature(x = "liger", dataset = "character"),
          function(x, dataset) {
              atac <- dataset(x, dataset)
              if (!inherits(atac, "ligerATACDataset")) {
                  stop("Specified dataset is not of ligerATACDataset class.")
              }
              atac@peak
          })

#' @export
setReplaceMethod(
    "peak",
    signature(x = "liger", dataset = "character"),
    function(x, dataset, value) {
        peak(dataset(x, dataset)) <- value
        x
    })

#' @export
setMethod("peak", signature(x = "ligerATACDataset", dataset = "missing"),
          function(x, dataset = NULL) {
              x@peak
          })

#' @export
setReplaceMethod(
    "peak",
    signature(x = "ligerATACDataset", dataset = "missing"),
    function(x, dataset = NULL, value) {
        if (!all(colnames(x) %in% colnames(value))) {
            stop("Cells in the dataset not all found in given matrix")
        }
        if (!identical(colnames(value), colnames(x))) {
            warning("colnames of given matrix not identical to the dataset. ",
                    "Reordering for the existing cells in the dataset.")
            value <- value[, colnames(x)]
        }
        x@peak <- value
        validObject(x)
        x
    })

