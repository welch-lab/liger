setClassUnion("dataframe", c("data.frame", "DFrame", "NULL", "missing"))
setClassUnion("matrixLike", c("matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"))

setGeneric("createLiger", function(raw.data = NULL,
                                   cell.meta = NULL,
                                   indices.name = NULL)
    standardGeneric("createLiger")
)
setGeneric("createLigerDataset",
           function(raw.data, format, modal = c("default", "rna", "atac")) {
               standardGeneric("createLigerDataset")
           })

setMethod(
    "createLigerDataset",
    signature(raw.data = "ligerDataset",
              format = "ANY",
              modal = "ANY"),
    function(raw.data, modal) {
        raw.data
    }
)

setMethod(
    "createLigerDataset",
    signature(raw.data = "matrixLike",
              format = "ANY",
              modal = "ANY"),
    function(raw.data, modal) {
        ligerDataset(raw.data, modal)
    }
)

setMethod(
    "createLigerDataset",
    signature(raw.data = "character",
              format = "ANY",
              modal = "ANY"),
    function(raw.data, format, modal) {
        message("Character input. Trying to read H5 file")
    }
)

#'
setMethod(
    "createLigerDataset",
    signature(raw.data = "Seurat",
              modal = "ANY"),
    function(raw.data, modal= c("default", "rna", "atac")) {
        counts <- Seurat::GetAssayData(raw.data, "counts")
        norm.data <- Seurat::GetAssayData(raw.data, "data")
        if (identical(counts, norm.data)) norm.data <- NULL
        scale.data <- Seurat::GetAssayData(raw.data, "scale.data")
        if (sum(dim(scale.data)) == 0) scale.data <- NULL
        ligerDataset(raw.data = counts, norm.data = norm.data,
                     scale.data = scale.data, modal = modal)
    }
)

setClass("anndata._core.anndata.AnnData")
setMethod(
    "createLigerDataset",
    signature(raw.data = "anndata._core.anndata.AnnData",
              modal = "ANY"),
    function(raw.data, modal) {
        message("Python object AnnData input. ")

    }
)
