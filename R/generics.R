#' @section Matrix access:
#' For \code{ligerDataset} object, \code{rawData()}, \code{normData},
#' \code{scaleData()} and \code{scaleUnsharedData()} methods are exported for
#' users to access the corresponding feature expression matrix. Replacement
#' methods are also available to modify the slots.
#'
#' For other matrices, such as the \eqn{H} and \eqn{V}, which are dataset
#' specific, please use \code{getMatrix()} method with specifying slot name.
#' Directly accessing slot with \code{@} is generally not recommended.
#' @export
#' @rdname ligerDataset-class
setGeneric("rawData", function(x, dataset = NULL) standardGeneric("rawData"))

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "rawData<-",
    function(x, dataset = NULL, check = TRUE, value) standardGeneric("rawData<-")
)

#' @export
#' @rdname ligerDataset-class
setGeneric("normData", function(x, dataset = NULL) standardGeneric("normData"))

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "normData<-",
    function(x, dataset = NULL, check = TRUE, value) standardGeneric("normData<-")
)

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "scaleData",
    function(x, dataset = NULL) standardGeneric("scaleData")
)

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "scaleData<-",
    function(x, dataset = NULL, check = TRUE, value) standardGeneric("scaleData<-")
)

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "scaleUnsharedData",
    function(x, dataset = NULL) standardGeneric("scaleUnsharedData")
)

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "scaleUnsharedData<-",
    function(x, dataset = NULL, check = TRUE, value) standardGeneric("scaleUnsharedData<-")
)

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "getMatrix",
    function(x, slot = "rawData", dataset = NULL, returnList = FALSE) {
        standardGeneric("getMatrix")
    }
)

#' @section H5 file and information access:
#' A \code{ligerDataset} object has a slot called \code{h5fileInfo}, which is a
#' list object. The first element is called \code{$H5File}, which is an
#' \code{H5File} class object and is the connection to the input file. The
#' second element is \code{$filename} which stores the absolute path of the H5
#' file in the current machine. The third element \code{$formatType} stores the
#' name of preset being used, if applicable. The other following keys pair with
#' paths in the H5 file that point to specific data for constructing a feature
#' expression matrix.
#'
#' \code{h5fileInfo()} method access the list described above and simply
#' retrieves the corresponding value. When \code{info = NULL}, returns the whole
#' list. When \code{length(info) == 1}, returns the requested list value. When
#' more info requested, returns a subset list.
#'
#' The replacement method modifies the list elements and corresponding slot
#' value (if applicable) at the same time. For example, running
#' \code{h5fileInfo(obj, "rawData") <- newPath} not only updates the list, but
#' also updates the \code{rawData} slot with the \code{H5D} class data at
#' "newPath" in the \code{H5File} object.
#'
#' \code{getH5File()} is a wrapper and is equivalent to
#' \code{h5fileInfo(obj, "H5File")}.
#' @export
#' @rdname ligerDataset-class
setGeneric("h5fileInfo", function(x, info = NULL) standardGeneric("h5fileInfo"))

#' @export
#' @rdname ligerDataset-class
setGeneric(
    "h5fileInfo<-",
    function(x, info = NULL, check = TRUE, value) {
        standardGeneric("h5fileInfo<-")
    }
)


#' @export
#' @rdname ligerDataset-class
setGeneric("getH5File", function(x, dataset = NULL) standardGeneric("getH5File"))

#' @export
#' @rdname ligerDataset-class
setMethod("getH5File",
          signature = signature(x = "ligerDataset", dataset = "missing"),
          function(x, dataset = NULL) h5fileInfo(x, "H5File"))


#' @section Feature metadata access:
#' A slot \code{featureMeta} is included for each \code{ligerDataset} object.
#' This slot requires a \code{\link[S4Vectors]{DataFrame-class}} object, which
#' is the same as \code{cellMeta} slot of a \linkS4class{liger} object. However,
#' the associated S4 methods only include access to the whole table for now.
#' Internal information access follows the same way as data.frame operation.
#' For example, \code{featureMeta(ligerD)$nCell} or
#' \code{featureMeta(ligerD)[varFeatures(ligerObj), "gene_var"]}.
#' @export
#' @rdname ligerDataset-class
setGeneric("featureMeta", function(x, check = NULL) {
    standardGeneric("featureMeta")
})

#' @export
#' @rdname ligerDataset-class
setGeneric("featureMeta<-", function(x, check = TRUE, value) {
    standardGeneric("featureMeta<-")
})







#' @section Dataset access:
#' \code{datasets()} method only accesses the \code{datasets} slot, the list of
#' \linkS4class{ligerDataset} objects. \code{dataset()} method accesses a single
#' dataset, with subsequent cell metadata updates and checks bonded when adding
#' or modifying a dataset. Therefore, when users want to modify something inside
#' a \code{ligerDataset} while no cell metadata change should happen, it is
#' recommended to use: \code{datasets(x)[[name]] <- ligerD} for efficiency,
#' though the result would be the same as \code{dataset(x, name) <- ligerD}.
#'
#' \code{length()} and \code{names()} methods are implemented to access the
#' number and names of datasets. \code{names<-} method is supported for
#' modifying dataset names, with taking care of the "dataset" variable in cell
#' metadata.
#' @section Matrix access:
#' For \code{liger} object, \code{rawData()}, \code{normData},
#' \code{scaleData()} and \code{scaleUnsharedData()} methods are exported for
#' users to access the corresponding feature expression matrix with
#' specification of one dataset. For retrieving a type of matrix from multiple
#' datasets, please use \code{getMatrix()} method.
#'
#' When only one matrix is expected to be retrieved by \code{getMatrix()}, the
#' matrix itself will be returned. A list will be returned if multiple matrices
#' is requested (by querying multiple datasets) or \code{returnList} is set to
#' \code{TRUE}.
#' @export
#' @rdname liger-class
setGeneric("datasets", function(x, check = NULL) standardGeneric("datasets"))

#' @export
#' @rdname liger-class
setGeneric(
    "datasets<-",
    function(x, check = TRUE, value) standardGeneric("datasets<-")
)

#' @export
#' @rdname liger-class
setGeneric("dataset", function(x, dataset = NULL) standardGeneric("dataset"))

#' @export
#' @rdname liger-class
setGeneric("dataset<-", function(x, dataset, type = NULL, qc = TRUE, value) {
    standardGeneric("dataset<-")
})


#' @export
#' @rdname liger-class
#' @section Cell metadata access:
#' Three approaches are provided for access of cell metadata. A generic function
#' \code{cellMeta} is implemented with plenty of options and multi-variable
#' accessibility. Besides, users can use double-bracket (e.g.
#' \code{ligerObj[[varName]]}) or dollor-sign (e.g. \code{ligerObj$nUMI}) to
#' access or modify single variables.
#'
#' For users' convenience of generating a customized ggplot with available cell
#' metadata, the S3 method \code{fortify.liger} is implemented. With this under
#' the hook, users can create simple ggplots by directly starting with
#' \code{ggplot(ligerObj, aes(...))} where cell metadata variables can be
#' directly thrown into \code{aes()}.
#'
#' Special partial metadata insertion is implemented specifically for mapping
#' categorical annotation from sub-population (subset object) back to original
#' experiment (full-size object). For example, when sub-clustering and
#' annotation is done for a specific cell-type of cells (stored in
#' \code{subobj}) subset from an experiment (stored as \code{obj}), users can do
#' \code{cellMeta(obj, "sub_ann", cellIdx = colnames(subobj)) <- subobj$sub_ann}
#' to map the value back, leaving other cells non-annotated with NAs. Plotting
#' with this variable will then also show NA cells with default grey color.
#' Furthermore, sub-clustering labels for other cell types can also be mapped
#' to the same variable. For example, \code{cellMeta(obj, "sub_ann",
#' cellIdx = colnames(subobj2)) <- subobj2$sub_ann}. As long as the labeling
#' variables are stored as factor class (categorical), the levels (category
#' names) will be properly handled and merged. Other situations follow the R
#' default behavior (e.g. categories might be converted to integer numbers if
#' mapped to numerical variable in the original object). Note that this feature
#' is only available with using the generic function \code{cellMeta} but not
#' with the \code{`[[`} or \code{`$`} accessing methods due to syntax reasons.
#'
#' The generic \code{defaultCluster} works as both getter and setter. As a
#' setter, users can do \code{defaultCluster(obj) <- "existingVariableName"} to
#' set a categorical variable as default cluster used for visualization or
#' downstream analysis. Users can also do \code{defaultCluster(obj,
#' "newVarName") <- factorOfLabels} to push new labeling into the object and set
#' as default. For getter method, the function returns a factor object of the
#' default cluster labeling. Argument \code{useDatasets} can be used for
#' requiring that given or retrieved labeling should match with cells in
#' specified datasets. We generally don't recommend setting \code{"dataset"} as
#' a default cluster because it is a preserved (always existing) field in
#' metadata and can lead to meaningless result when running analysis that
#' utilizes both clustering information and the dataset source information.
setGeneric(
    "cellMeta",
    function(x, columns = NULL, useDatasets = NULL, cellIdx = NULL, as.data.frame = FALSE, ...) {
        standardGeneric("cellMeta")
    }
)

#' @export
#' @rdname liger-class
setGeneric(
    "cellMeta<-",
    function(x, columns = NULL, useDatasets = NULL, cellIdx = NULL, inplace = FALSE, check = FALSE, value) {
        standardGeneric("cellMeta<-")
    }
)

#' @export
#' @rdname liger-class
setGeneric(
    "defaultCluster",
    function(x, useDatasets = NULL, ...) {
        standardGeneric("defaultCluster")
    }
)

#' @export
#' @rdname liger-class
setGeneric(
    "defaultCluster<-",
    function(x, name = NULL, useDatasets = NULL, ..., value) {
        standardGeneric("defaultCluster<-")
    }
)

#' @export
#' @rdname liger-class
#' @section Dimension reduction access:
#' Currently, low-dimensional representaion of cells, presented as dense
#' matrices, are all stored in \code{dimReds} slot, and can totally be accessed
#' with generics \code{dimRed} and \code{dimRed<-}. Adding a dimRed to the
#' object looks as simple as \code{dimRed(obj, "name") <- matrixLike}. It can
#' be retrieved back with \code{dimRed(obj, "name")}. Similar to having a
#' default cluster labeling, we also constructed the feature of default dimRed.
#' It can be set with \code{defaultDimRed(obj) <- "existingMatLikeVar"} and the
#' matrix can be retrieved with \code{defaultDimRed(obj)}.
setGeneric(
    "dimReds",
    function(x) standardGeneric("dimReds")
)

#' @export
#' @rdname liger-class
setGeneric(
    "dimReds<-",
    function(x, value) standardGeneric("dimReds<-")
)

#' @export
#' @rdname liger-class
setGeneric(
    "dimRed",
    function(x, name = NULL, useDatasets = NULL, cellIdx = NULL, ...) {
        standardGeneric("dimRed")
    }
)

#' @export
#' @rdname liger-class
setGeneric(
    "dimRed<-",
    function(x, name = NULL, useDatasets = NULL, cellIdx = NULL, ..., value) {
        standardGeneric("dimRed<-")
    }
)

#' @export
#' @rdname liger-class
setGeneric(
    "defaultDimRed",
    function(x, useDatasets = NULL, cellIdx = NULL) {
        standardGeneric("defaultDimRed")
    }
)

#' @export
#' @rdname liger-class
setGeneric(
    "defaultDimRed<-",
    function(x, value) {
        standardGeneric("defaultDimRed<-")
    }
)

#' @export
#' @rdname liger-class
#' @section Variable feature access:
#' The \code{varFeatures} slot allows for character vectors of gene names.
#' \code{varFeatures(x)} returns this vector and \code{value} for
#' \code{varFeatures<-} method has to be a character vector or \code{NULL}.
#' The replacement method, when \code{check = TRUE} performs checks on gene
#' name consistency check across the \code{scaleData}, \code{H}, \code{V} slots
#' of inner \code{ligerDataset} objects as well as the \code{W} and
#' \code{H.norm} slots of the input \code{liger} object.
setGeneric("varFeatures", function(x) standardGeneric("varFeatures"))

#' @export
#' @rdname liger-class
setGeneric(
    "varFeatures<-",
    function(x, check = TRUE, value) standardGeneric("varFeatures<-")
)



#' @export
#' @rdname liger-class
setGeneric("varUnsharedFeatures", function(x, dataset = NULL) {
    standardGeneric("varUnsharedFeatures")
})

#' @export
#' @rdname liger-class
setGeneric(
    "varUnsharedFeatures<-",
    function(x, dataset, check = TRUE, value) {
        standardGeneric("varUnsharedFeatures<-")
    }
)

#' @section Command records:
#' rliger functions, that perform calculation and update the \code{liger}
#' object, will be recorded in a \code{ligerCommand} object and stored in the
#' \code{commands} slot, a list, of \code{liger} object. Method
#' \code{commands()} is implemented to retrieve or show the log history.
#' Running with \code{funcName = NULL} (default) returns all command labels.
#' Specifying \code{funcName} allows partial matching to all command labels
#' and returns a subset list (of \code{ligerCommand} object) of matches (or
#' the \code{ligerCommand} object if only one match found). If \code{arg} is
#' further specified, a subset list of parameters from the matches will be
#' returned. For example, requesting a list of resolution values used in
#' all louvain cluster attempts: \code{commands(ligerObj, "louvainCluster",
#' "resolution")}
#' @export
#' @rdname liger-class
setGeneric(
    "commands",
    function(x, funcName = NULL, arg = NULL) standardGeneric("commands")
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






#' Converting other classes of data to a liger object
#' @description
#' This function converts data stored in SingleCellExperiment (SCE), Seurat
#' object or a merged sparse matrix (dgCMatrix) into a liger object. This is
#' designed for a container object or matrix that already contains multiple
#' datasets to be integerated with LIGER. For individual datasets, please use
#' \code{\link{createLiger}} instead.
#' @export
#' @param object Object.
#' @param datasetVar Specify the dataset belonging by: 1. Select a variable from
#' existing metadata in the object (e.g. colData column); 2. Specify a
#' vector/factor that assign the dataset belonging. 3. Give a single character
#' string which means that all data is from one dataset (must not be a metadata
#' variable, otherwise it is understood as 1.). Default \code{NULL} gathers
#' things into one dataset and names it "sample" for dgCMatrix, attempts
#' to find variable "sample" from SCE or "orig.ident" from Seurat.
#' @param modal Modality setting for each dataset. See
#' \code{\link{createLiger}}.
#' @param ... Additional arguments passed to \code{\link{createLiger}}
#' @details
#' For Seurat V5 structure, it is highly recommended that users make use of its
#' split layer feature, where things like "counts", "data", and "scale.data"
#' can be held for each dataset in the same Seurat object, e.g. with
#' "count.ctrl", "count.stim", not merged. If a Seurat object with split layers
#' is given, \code{datasetVar} will be ignored and the layers will be directly
#' used.
#' @return a \linkS4class{liger} object.
#' @rdname as.liger
#' @examples
#' # dgCMatrix (common sparse matrix class), usually obtained from other
#' # container object, and contains multiple samples merged in one.
#' matList <- rawData(pbmc)
#' multiSampleMatrix <- mergeSparseAll(matList)
#' # The `datasetVar` argument expects the variable assigning the sample source
#' pbmc2 <- as.liger(multiSampleMatrix, datasetVar = pbmc$dataset)
#' pbmc2
#'
#' \donttest{
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'     sce <- SingleCellExperiment::SingleCellExperiment(
#'         assays = list(counts = multiSampleMatrix)
#'     )
#'     sce$sample <- pbmc$dataset
#'     pbmc3 <- as.liger(sce, datasetVar = "sample")
#'     pbmc3
#' }
#'
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'     seu <- SeuratObject::CreateSeuratObject(multiSampleMatrix)
#'     # Seurat creates variable "orig.ident" by identifying the cell barcode
#'     # prefixes, which is indeed what we need in this case. Users might need
#'     # to be careful and have it confirmed first.
#'     pbmc4 <- as.liger(seu, datasetVar = "orig.ident")
#'     pbmc4
#'
#'     # As per Seurat V5 updates with layered data, specifically helpful udner the
#'     # scenario of dataset integration. "counts" and etc for each datasets can be
#'     # split into layers.
#'     seu5 <- seu
#'     seu5[["RNA"]] <- split(seu5[["RNA"]], pbmc$dataset)
#'     print(SeuratObject::Layers(seu5))
#'     pbmc5 <- as.liger(seu5)
#'     pbmc5
#' }
#' }
as.liger <- function(object, ...) UseMethod("as.liger", object)

#' Converting other classes of data to a ligerDataset object
#' @description
#' Works for converting a matrix or container object to a single ligerDataset,
#' and can also convert the modality preset of a ligerDataset. When used with
#' a dense matrix object, it automatically converts the matrix to sparse form
#' (\code{\link[Matrix]{dgCMatrix-class}}). When used with container objects
#' such as Seurat or SingleCellExperiment, it is highly recommended that the
#' object contains only one dataset/sample which is going to be integrated with
#' LIGER. For multi-sample objects, please use \code{\link{as.liger}} with
#' dataset source variable specified.
#' @export
#' @param object Object.
#' @param modal Modality setting for each dataset. Choose from \code{"default"},
#' \code{"rna"}, \code{"atac"}, \code{"spatial"}, \code{"meth"}.
#' @param ... Additional arguments passed to \code{\link{createLigerDataset}}
#' @return a \linkS4class{liger} object.
#' @rdname as.ligerDataset
#' @examples
#' ctrl <- dataset(pbmc, "ctrl")
#' ctrl
#' # Convert the modality preset
#' as.ligerDataset(ctrl, modal = "atac")
#' rawCounts <- rawData(ctrl)
#' class(rawCounts)
#' as.ligerDataset(rawCounts)
as.ligerDataset <- function(object, ...) UseMethod("as.ligerDataset", object)
