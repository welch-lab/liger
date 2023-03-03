setClass("ligerDataset")
setClassUnion("matrixLike", c("matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"))

##############################################
# From other things to ligerDataset class ####
##############################################

#' Convert an object of various classes to a ligerDataset object
#' @rdname as.ligerDataset
#' @description S4 method to convert objects of various types to a
#' \linkS4class{ligerDataset} object or a modality specific sub-class of
#' \linkS4class{ligerDataset} object. Supported classes include a matrix like
#' object, a \linkS4class{SingleCellExperiment} object, a \linkS4class{Seurat}
#' object, and \code{AnnData} object. This method also supports modality
#' setting to a \linkS4class{ligerDataset} inherited object.
#' @param x An object to be converted
#' @param modal The modality of this dataset. Default \code{"default"} for RNA.
#' Can choose from \code{"rna"}, \code{"atac"}.
#' @return A ligerDataset object by default, or a modality specific sub-class of
#' ligerDataset object according to \code{modal}.
#' @export
#' @author Yichen Wang
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
        if (!requireNamespace("SingleCellExperiment",quietly = "TRUE"))
            stop("Package \"SingleCellExperiment\" needed for this function to work. ",
                 "Please install it by command:\n",
                 "BiocManager::install('SingleCellExperiment')",
                 call. = FALSE)
        if ("counts" %in% SummarizedExperiment::assayNames(x))
            raw.data <- SingleCellExperiment::counts(x)
        else raw.data <- NULL
        if ("logcounts" %in% SummarizedExperiment::assayNames(x))
            norm.data <- SingleCellExperiment::logcounts(x)
        else norm.data <- NULL
        if ("counts" %in% SummarizedExperiment::assayNames(x))
            raw.data <- SingleCellExperiment::counts(x)
        createLigerDataset(raw.data = raw.data, norm.data = norm.data,
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
    function(x, modal= c("default", "rna", "atac")) {
        if (!requireNamespace("Seurat",quietly = "TRUE"))
            stop("Package \"Seurat\" needed for this function to work. ",
                 "Please install it by command:\n",
                 "BiocManager::install('Seurat')",
                 call. = FALSE)
        counts <- Seurat::GetAssayData(x, "counts")
        norm.data <- Seurat::GetAssayData(x, "data")
        if (identical(counts, norm.data)) norm.data <- NULL
        scale.data <- Seurat::GetAssayData(x, "scale.data")
        if (sum(dim(scale.data)) == 0) scale.data <- NULL
        createLigerDataset(raw.data = counts, norm.data = norm.data,
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

##############################################
# From ligerDataset class to other things ####
##############################################

#setAs("ligerDataset", "SingleCellExperiment", function(from) {
#    requireNamespace("SingleCellExperiment")
#    assays <- list()
#    if (!is.null(raw.data(from)))
#        assays <- c(assays, list(counts = raw.data(from)))
#    if (!is.null(norm.data(from)))
#        assays <- c(assays, list(normcounts = norm.data(from)))
#    SingleCellExperiment::SingleCellExperiment(
#        assays = assays
#    )
#})
