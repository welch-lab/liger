setClassUnion("dataframe", c("data.frame", "DFrame", "NULL", "missing"))
setClassUnion("matrixLike", c("matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"))

setGeneric("createLiger", function(raw.data = NULL,
                                   cell.meta = NULL,
                                   indices.name = NULL)
    standardGeneric("createLiger")
)
setGeneric("createLigerDataset", function(raw.data) standardGeneric("createLigerDataset"))

setMethod(
    "createLiger",
    signature(raw.data = "list",
              cell.meta = "dataframe",
              indices.name = "missing"),
    function(raw.data, cell.meta = NULL) {
        print("yes1")
    }
)

setMethod(
    "createLigerDataset",
    signature(raw.data = "matrixLike"),
    function(raw.data) {
        ligerDataset(raw.data)
    }
)

setMethod(
    "createLigerDataset",
    signature(raw.data = "character"),
    function(raw.data) {
        message("Character input. Trying to read H5 file")
    }
)
