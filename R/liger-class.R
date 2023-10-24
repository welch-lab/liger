setClassUnion("dgCMatrix_OR_NULL", c("dgCMatrix", "NULL"))
setClassUnion("matrix_OR_NULL", c("matrix", "NULL"))
setClassUnion("matrixLike_OR_NULL", c(
    "matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix", "NULL"
))
setClassUnion("character_OR_NULL", c("character", "NULL"))
setClassUnion("matrixLike", c(
    "matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"
))
setClassUnion("Number_or_NULL", c("integer", "numeric", "NULL"))

setClassUnion("dataframe", c("data.frame", "DataFrame", "NULL", "missing"))

#' @title liger class
#' @rdname liger-class
#' @docType class
#' @description \code{liger} object is the main data container for LIGER
#' analysis in R. The slot \code{datasets} is a list where each element should
#' be a \linkS4class{ligerDataset} object containing dataset specific
#' information, such as the expression matrices. The other parts of liger object
#' stores information that can be shared across the analysis, such as the cell
#' metadata and factorization result matrices.
#'
#' This manual provides explanation to the \code{liger} object structure as well
#' as usage of class-specific methods. Please see detail sections for more
#' information.
#'
#' For \code{liger} objects created with older versions of rliger package,
#' please try updating the objects individually with
#' \code{\link{convertOldLiger}}.
#' @slot datasets list of \linkS4class{ligerDataset} objects. Use generic
#' \code{dataset}, \code{dataset<-}, \code{datasets} or \code{datasets<-} to
#' interact with. See detailed section accordingly.
#' @slot cellMeta \linkS4class{DFrame} object for cell metadata. Pre-existing
#' metadata, QC metrics, cluster labeling, low-dimensional embedding and etc.
#' are all stored here. Use generic \code{cellMeta}, \code{cellMeta<-},
#' \code{$}, \code{[[]]} or \code{[[]]<-} to interact with. See detailed section
#' accordingly.
#' @slot varFeatures Character vector of feature names. Use generic
#' \code{varFeatures} or \code{varFeatures<-} to interact with. See detailed
#' section accordingly.
#' @slot W Matrix of gene loading for each factor. See
#' \code{\link{runIntegration}}.
#' @slot H.norm Matrix of aligned factor loading for each cell. See
#' \code{\link{quantileNorm}} and \code{\link{runIntegration}}.
#' @slot commands List of \linkS4class{ligerCommand} objects. Record of
#' analysis. Use \code{commands} to retrieve information. See detailed section
#' accordingly.
#' @slot uns List for unstructured meta-info of analyses or presets.
#' @slot version Record of version of rliger2 package
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom ggplot2 fortify
liger <- setClass(
    "liger",
    representation(
        datasets = "list",
        cellMeta = "DataFrame",
        varFeatures = "character_OR_NULL",
        W = "matrix_OR_NULL",
        H.norm = "matrix_OR_NULL",
        uns = "list",
        commands = "list",
        version = "ANY"
    ),
    methods::prototype(
        cellMeta = methods::new("DFrame"),
        version = utils::packageVersion("rliger2")
    )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Validity ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.checkLigerBarcodes <- function(x) {
    bcFromDatasets <- unlist(lapply(datasets(x), colnames), use.names = FALSE)
    if (!identical(colnames(x), bcFromDatasets)) {
        return("liger object barcodes do not match to barcodes in datasets")
    }
    if (!is.null(x@H.norm)) {
        if (!identical(rownames(x@H.norm), bcFromDatasets)) {
            return("H.norm barcodes do not match to barcodes in datasets.")
        }
    }
    if (!"dataset" %in% names(cellMeta(x))) {
        return("`datasets` variable missing in cellMeta(x)")
    }
    datasetNamesFromDatasets <- rep(names(x), lapply(datasets(x), ncol))
    names(datasetNamesFromDatasets) <- NULL

    if (!identical(datasetNamesFromDatasets, as.character(x$dataset))) {
        return("names of datasets do not match
               `datasets` variable in cellMeta")
    }
    return(NULL)
}

.checkLigerVarFeature <- function(x) {
    if (!is.null(varFeatures(x)) &&
        length(varFeatures(x)) > 0) {
        if (!is.null(x@W))
            if (!identical(rownames(x@W), varFeatures(x)))
                return("Variable features do not match dimension of W matrix")
        for (d in names(x)) {
            ld <- dataset(x, d)
            if (!is.null(ld@V)) {
                if (!identical(rownames(ld@V), varFeatures(x)))
                    return(paste("Variable features do not match dimension",
                                 "of V matrix in dataset", d))
            }

            if (!is.null(scaleData(ld))) {
                if (!isH5Liger(ld)) {
                    if (!identical(rownames(scaleData(ld)), varFeatures(x)))
                        return(paste("Variable features do not match dimension",
                                     "of scaleData in dataset", d))
                } else {
                    if (inherits(scaleData(ld), "H5D")) {
                        if (scaleData(ld)$dims[1] != length(varFeatures(x)))
                            return(paste("Variable features do not match ",
                                         "dimension of scaleData in dataset ",
                                         "(H5)", d))
                    } else if (inherits(scaleData(ld), "H5Group")) {
                        if (scaleData(ld)[["featureIdx"]]$dims != length(varFeatures(x))) {
                            return(paste("Variable features do not match ",
                                         "dimension of scaleData in dataset ",
                                         "(H5)", d))
                        }
                    }
                }
            }
        }
    }
    return(NULL)
}

.valid.liger <- function(object) {
    # message("Checking liger object validity")
    res <- .checkLigerBarcodes(object)
    if (!is.null(res)) return(res)
    res <- .checkLigerVarFeature(object)
    if (!is.null(res)) return(res)
    # TODO more checks
}

setValidity("liger", .valid.liger)

#' Check if given liger object if under new implementation
#' @param object A liger object
#' @return \code{TRUE} if the version of \code{object} is later than or equal to
#' 1.99.0. Otherwise \code{FALSE}
#' @noRd
is.newLiger <- function(object) {
    v <- object@version
    v1990 <- package_version("1.99.0")
    if (v >= v1990) TRUE
    else FALSE
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Base Generic Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @param x,object,model A \linkS4class{liger} object
#' @param dataset Name or numeric index of a dataset
#' @param value Check detail sections for requirements.
#' @param type When using \code{dataset<-} with a matrix like \code{value},
#' specify what type the matrix is. Choose from \code{"rawData"},
#' \code{"normData"} or \code{"scaleData"}.
#' @param qc Logical, whether to perform general qc on added new dataset.
#' @param check Logical, whether to perform object validity check on setting new
#' value. Users are not supposed to set \code{FALSE} here.
#' @param columns The names of available variables in \code{cellMeta} slot. When
#' \code{as.data.frame = TRUE}, please use variable names after coercion.
#' @param name The name of available variables in \code{cellMeta} slot or the
#' name of a new variable to store.
#' @param cellIdx Valid cell subscription to subset retrieved variables. Default
#' \code{NULL} uses all cells.
#' @param as.data.frame Logical, whether to apply
#' \code{\link[base]{as.data.frame}} on the subscription. Default \code{FALSE}.
#' @param i,j Feature and cell index for \code{`[`} method. For \code{`[[`}
#' method, use a single variable name with \code{i} while \code{j} is not
#' applicable.
#' @param drop Not applicable.
#' @param slot Name of slot to retrieve matrix from. Options shown in Usage.
#' @param returnList Logical, whether to force return a list even when only one
#' dataset-specific matrix (i.e. expression matrices, H, V or U) is requested.
#' Default \code{FALSE}.
#' @param useDatasets Setter or getter method should only apply on cells in
#' specified datasets. Any valid character, numeric or logical subscriber is
#' acceptable. Default \code{NULL} works with all datasets.
#' @param funcName,arg See Command records section.
#' @param data fortify method required argument. Not used.
#' @param ... See detailed sections for explanation.
#' @return See detailed sections for explanetion.
#' @export
#' @rdname liger-class
#' @examples
#' # Methods for base generics
#' pbmcPlot
#' print(pbmcPlot)
#' dim(pbmcPlot)
#' ncol(pbmcPlot)
#' colnames(pbmcPlot)[1:5]
#' pbmcPlot[varFeatures(pbmcPlot)[1:10], 1:10]
#' names(pbmcPlot)
#' length(pbmcPlot)
#'
#' # rliger2 generics
#' ## Retrieving dataset(s), replacement methods available
#' datasets(pbmcPlot)
#' dataset(pbmcPlot, "ctrl")
#' dataset(pbmcPlot, 2)
#'
#' ## Retrieving cell metadata, replacement methods available
#' cellMeta(pbmcPlot)
#' head(pbmcPlot[["nUMI"]])
#' head(pbmcPlot$UMAP)
#'
#' ## Retrieving variable features, replacement methods available
#' varFeatures(pbmcPlot)
#'
#' ## Command record/history
#' pbmcPlot <- scaleNotCenter(pbmcPlot)
#' commands(pbmcPlot)
#' commands(pbmcPlot, funcName = "scaleNotCenter")
#'
#' # S3 methods
#' c(pbmcPlot, pbmcPlot)
#'
#' library(ggplot2)
#' ggplot(pbmcPlot, aes(x = UMAP.1, y = UMAP.2)) + geom_point()
setMethod(
    f = "show",
    signature(object = "liger"),
    definition = function(object) {
        .checkObjVersion(object)
        cat("An object of class liger with", ncol(object), "cells\n")
        cat(paste0("datasets(", length(object), "): "))
        datasetInfos <- paste0(names(object), " (",
                               unlist(lapply(datasets(object), ncol)),
                               " cells)")
        cat(.collapseLongNames(datasetInfos), "\n")
        cat(paste0("cellMeta(", ncol(cellMeta(object)), "): "))
        cat(.collapseLongNames(colnames(cellMeta(object))), "\n")
        cat(paste0("varFeatures(", length(varFeatures(object)), "): "))
        cat(.collapseLongNames(varFeatures(object)), "\n")
        invisible(x = NULL)
    }
)

#' @section Dimensionality:
#' For a \code{liger} object, the column orientation is assigned for
#' cells. Due to the data structure, it is hard to define a row index for the
#' \code{liger} object, which might contain datasets that vary in number of
#' genes.
#'
#' Therefore, for \code{liger} objects, \code{dim} and \code{dimnames} returns
#' \code{NA}/\code{NULL} for rows and total cell counts/barcodes for the
#' columns.
#'
#' For direct call of \code{dimnames<-} method, \code{value} should be a list
#' with \code{NULL} as the first element and valid cell identifiers as the
#' second element. For \code{colnames<-} method, the character vector of cell
#' identifiers. \code{rownames<-} method is not applicable.
#' @rdname liger-class
#' @export
setMethod("dim", "liger", function(x) {
    c(NA, nrow(cellMeta(x)))
})

#' @rdname liger-class
#' @export
setMethod("dimnames", "liger", function(x) {
    list(NULL, rownames(cellMeta(x)))
})

#' @rdname liger-class
#' @export
setReplaceMethod("dimnames", c("liger", "list"), function(x, value) {
    dataset.var <- x$dataset
    dataset.uniq <- unique(dataset.var)
    for (d in dataset.uniq) {
        dataset.idx <- dataset.var == d
        # Not using `dataset()` accessor because it does object validity check
        colnames(x@datasets[[d]]) <- value[[2L]][dataset.idx]
    }
    rownames(cellMeta(x)) <- value[[2L]]
    if (!is.null(x@H.norm)) rownames(x@H.norm) <- value[[2L]]
    x
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subsetting ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @section Subsetting:
#' For more detail of subsetting a \code{liger} object or a
#' \linkS4class{ligerDataset} object, please check out \code{\link{subsetLiger}}
#' and \code{\link{subsetLigerDataset}}. Here, we set the S4 method
#' "single-bracket" \code{[} as a quick wrapper to subset a \code{liger} object.
#' Note that \code{j} serves as cell subscriptor which can be any valid index
#' refering the collection of all cells (i.e. \code{rownames(cellMeta(obj))}).
#' While \code{i}, the feature subscriptor can only be character vector because
#' the features for each dataset can vary. \code{...} arugments are passed to
#' \code{subsetLiger} so that advanced options are allowed.
#' @export
#' @rdname liger-class
setMethod(
    "[",
    signature(x = "liger", i = "character", j = "missing"),
    function(x, i, j, ...) subsetLiger(x, featureIdx = i, cellIdx = NULL, ...)
)

#' @export
#' @rdname liger-class
setMethod(
    "[",
    signature(x = "liger", i = "missing", j = "index"),
    function(x, i, j, ...) subsetLiger(x, featureIdx = NULL, cellIdx = j, ...)
)

#' @export
#' @rdname liger-class
setMethod(
    "[",
    signature(x = "liger", i = "character", j = "index"),
    function(x, i, j, ...) subsetLiger(x, featureIdx = i, cellIdx = j, ...)
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Datasets ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
setMethod("datasets", signature(x = "liger", check = "ANY"),
          function(x, check = NULL) {
              x@datasets
    })

#' @export
#' @rdname liger-class
setReplaceMethod("datasets", signature(x = "liger", check = "logical"),
                 function(x, check = TRUE, value) {
                     x@datasets <- value
                     if (isTRUE(check)) methods::validObject(x)
                     x
                 })

#' @export
#' @rdname liger-class
setReplaceMethod("datasets", signature(x = "liger", check = "missing"),
                 function(x, check = TRUE, value) {
                     x@datasets <- value
                     methods::validObject(x)
                     x
                 })

#' @export
#' @rdname liger-class
setGeneric("dataset", function(x, dataset = NULL) standardGeneric("dataset"))

#' @export
#' @rdname liger-class
setGeneric("dataset<-", function(x, dataset, type = NULL, qc = FALSE, value) {
    standardGeneric("dataset<-")
})

#' @export
#' @rdname liger-class
setMethod("dataset", signature(x = "liger", dataset = "character_OR_NULL"),
          function(x, dataset = NULL) {
              if (is.null(dataset)) return(datasets(x)[[1]])
              else {
                  if (!dataset %in% names(x)) {
                      stop('Specified dataset name "', dataset,
                           '" not found in liger object.')
                  }
                  return(datasets(x)[[dataset]])
              }
          })

#' @export
#' @rdname liger-class
setMethod("dataset", signature(x = "liger", dataset = "missing"),
          function(x, dataset = NULL) {
              datasets(x)[[1]]
          })

#' @export
#' @rdname liger-class
setMethod("dataset", signature(x = "liger", dataset = "numeric"),
          function(x, dataset = NULL) {
              datasets(x)[[dataset]]
          })

.mergeCellMeta <- function(cm1, cm2) {
    newDF <- S4Vectors::DataFrame(row.names = c(rownames(cm1), rownames(cm2)))
    for (var in names(cm1)) {
        value <- cm1[[var]]
        if (var %in% names(cm2)) {
            # TODO: check type
            tryCatch(
                expr = {
                    if (is.null(dim(value))) {
                        value <- c(value, cm2[[var]])
                    } else {
                        value <- rbind(value, cm2[[var]])
                    }
                },
                finally = {
                    cm2Idx <- seq(nrow(cm1) + 1, nrow(cm1) + nrow(cm2))
                    if (is.null(dim(value))) value[cm2Idx] <- NA
                    else {
                        empty <- matrix(NA, nrow = nrow(cm2), ncol = ncol(value))
                        value <- rbind(value, empty)
                    }
                }
            )
        } else {
            cm2Idx <- seq(nrow(cm1) + 1, nrow(cm1) + nrow(cm2))
            if (is.null(dim(value))) value[cm2Idx] <- NA
            else {
                empty <- matrix(NA, nrow = nrow(cm2), ncol = ncol(value))
                value <- rbind(value, empty)
            }
        }
        newDF[[var]] <- value
    }
    for (var in names(cm2)[!names(cm2) %in% names(cm1)]) {
        value <- cm2[[var]]
        if (is.null(dim(value))) {
            if (is.factor(value)) {
                value <- factor(c(rep(NA, nrow(cm1)), value),
                                levels = levels(value))
            } else {
                value <- c(rep(NA, nrow(cm1)), value)
            }
        } else {
            empty <- matrix(NA, nrow = nrow(cm1), ncol = ncol(value))
            value <- rbind(empty, value)
        }
        newDF[[var]] <- value
    }
    return(newDF)
}
.expandDataFrame <- function(df, idx) {
    dfList <- as.list(df)
    dfList <- lapply(dfList, function(x, idx) {
        if (is.null(dim(x))) {
            x[idx] <- NA
        } else {
            empty <- matrix(NA, nrow = length(idx), ncol = ncol(x))
            x <- rbind(x, empty)
        }
        x
    }, idx = idx)
    include <- sapply(dfList, function(x) {
        is.vector(x) | is.factor(x)
    })
    newdf <- S4Vectors::DataFrame(dfList[include])
    for (i in which(!include)) newdf[[names(dfList)[i]]] <- dfList[[i]]
    newdf
}

#' @export
#' @rdname liger-class
setReplaceMethod("dataset", signature(x = "liger", dataset = "character",
                                      type = "missing", qc = "ANY",
                                      value = "ligerDataset"),
                 function(x, dataset, type = NULL, qc = FALSE, value) {
                     if (dataset %in% names(x)) {
                         dataset(x, dataset) <- NULL
                     }
                     new.idx <- seq(ncol(x) + 1, ncol(x) + ncol(value))
                     x@datasets[[dataset]] <- value
                     cm <- x@cellMeta
                     remainingRowname <- rownames(cm)
                     cm <- .expandDataFrame(cm, new.idx)
                     rownames(cm) <- c(remainingRowname, colnames(value))

                     levels(cm$dataset) <- c(levels(cm$dataset), dataset)
                     cm$dataset[new.idx] <- dataset
                     cm$barcode[new.idx] <- colnames(value)
                     x@cellMeta <- cm
                     # x@W is genes x k, no need to worry
                     if (!is.null(x@H.norm)) {
                         message("Filling in NAs to H.norm matrix")
                         H.normNew <- matrix(
                             NA, ncol(value), ncol(x@H.norm),
                             dimnames = list(colnames(value), NULL))
                         x@H.norm <- rbind(x@H.norm, H.normNew)
                     }
                     methods::validObject(x)
                     if (qc) x <- runGeneralQC(x, useDatasets = dataset)
                     x
                 })

#' @export
#' @rdname liger-class
setReplaceMethod("dataset", signature(x = "liger", dataset = "character",
                                      type = "ANY", qc = "ANY",
                                      value = "matrixLike"),
                 function(x, dataset,
                          type = c("rawData", "normData", "scaleData"),
                          qc = FALSE,
                          value) {
                     type <- match.arg(type)
                     if (type == "rawData") {
                         ld <- createLigerDataset(rawData = value)
                         colnames(ld) <- paste0(dataset, "_", colnames(ld))
                     } else if (type == "normData") {
                         ld <- createLigerDataset(normData = value)
                     } else if (type == "scaleData") {
                         ld <- createLigerDataset(scaleData = value)
                     }
                     dataset(x, dataset, qc = qc) <- ld
                     x
                 })

#' @export
#' @rdname liger-class
setReplaceMethod("dataset", signature(x = "liger", dataset = "character",
                                      type = "missing", qc = "ANY",
                                      value = "NULL"),
                 function(x, dataset, type = NULL, qc = TRUE, value) {
                     if (!dataset %in% names(x)) {
                         warning("Specified dataset name not found in ",
                                 "liger object. Nothing would happen.")
                     } else {
                         idxToRemove <- x$dataset == dataset
                         x@datasets[[dataset]] <- NULL
                         x@cellMeta <- x@cellMeta[!idxToRemove, , drop = FALSE]
                         x@H.norm <- x@H.norm[!idxToRemove, , drop = FALSE]
                         x@cellMeta$dataset <- droplevels(x@cellMeta$dataset)
                     }
                     x
                 })

#' @rdname liger-class
#' @export
setMethod("names", signature(x = "liger"), function(x) {
    names(datasets(x))
})

#' @rdname liger-class
#' @export
setReplaceMethod(
    "names",
    signature(x = "liger", value = "character"),
    function(x, value) {
        originalNames <- names(x)
        if (!identical(value, originalNames)) {
            dataset.idx <- lapply(originalNames, function(n) {
                x$dataset == n
            })
            x@cellMeta$dataset <- as.character(x@cellMeta$dataset)
            for (i in seq_along(value)) {
                x@cellMeta$dataset[dataset.idx[[i]]] <- value[i]
            }
            x@cellMeta$dataset <- factor(x@cellMeta$dataset, levels = value)
            names(x@datasets) <- value
        }
        x
    })

#' @rdname liger-class
#' @export
setMethod("length", signature(x = "liger"), function(x) {
    .checkObjVersion(x)
    length(datasets(x))
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Cell metadata ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    function(x, columns = NULL, cellIdx = NULL, as.data.frame = FALSE, ...) {
        standardGeneric("cellMeta")
    }
)

#' @export
#' @rdname liger-class
setGeneric(
    "cellMeta<-",
    function(x, columns = NULL, check = FALSE, value) {
        standardGeneric("cellMeta<-")
    }
)

.subsetCellMeta <- function(
        object,
        columns = NULL,
        cellIdx = NULL,
        as.data.frame = FALSE,
        ...
) {
    res <- object@cellMeta
    if (isTRUE(as.data.frame)) res <- as.data.frame(res)
    if (!is.null(columns)) {
        notFound <- !columns %in% colnames(res)
        if (any(notFound)) {
            warning("Specified variables from cellMeta not found: ",
                    paste(columns[notFound], collapse = ", "))
            columns <- columns[!notFound]
        }
        res <- res[, columns, ...]
    }
    if (!is.null(cellIdx)) {
        cellIdx <- .idxCheck(object, idx = cellIdx, orient = "cell")
        if (is.vector(res) || is.factor(res)) res <- res[cellIdx]
        else if (!is.null(dim(res))) res <- res[cellIdx, , ...]
        else stop("Result before idx subscription corrupted")
    }
    return(res)
}

#' @export
#' @rdname liger-class
setMethod(
    "cellMeta",
    signature(x = "liger", columns = "NULL"),
    function(x, columns = NULL, cellIdx = NULL, as.data.frame = FALSE, ...) {
        NULL
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "cellMeta",
    signature(x = "liger", columns = "character"),
    function(x, columns = NULL, cellIdx = NULL, as.data.frame = FALSE, ...) {
        .subsetCellMeta(x, columns = columns, cellIdx = cellIdx,
                        as.data.frame = as.data.frame, ...)
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "cellMeta",
    signature(x = "liger", columns = "missing"),
    function(x, columns = NULL, cellIdx = NULL, as.data.frame = FALSE, ...) {
        .subsetCellMeta(x, columns = NULL, cellIdx = cellIdx,
                        as.data.frame = as.data.frame, ...)
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "cellMeta",
    signature(x = "liger", columns = "missing"),
    function(x, columns = NULL, check = FALSE, value) {
        if (!inherits(value, "DFrame"))
            value <- S4Vectors::DataFrame(value)
        x@cellMeta <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "cellMeta",
    signature(x = "liger", columns = "character"),
    function(x, columns = NULL, check = FALSE, value) {
        x@cellMeta[[columns]] <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @rdname liger-class
#' @export
setMethod(
    "[[",
    c("liger", "ANY", "missing"),
    function(x, i, j, ...) {
        if (is.null(i)) NULL
        else if (!i %in% names(cellMeta(x))) NULL
        else cellMeta(x, columns = i, ...)
    }
)

#' @rdname liger-class
#' @export
setReplaceMethod(
    "[[",
    c("liger", "ANY", "missing"),
    function(x, i, j, ..., value) {
        cellMeta(x, columns = i, ...) <- value
        x
    }
)

#' @export
#' @method .DollarNames liger
.DollarNames.liger <- function(x, pattern = "") {
    grep(pattern, colnames(x@cellMeta), value = TRUE)
}


#' @export
#' @rdname liger-class
setMethod(
    "$",
    signature(x = "liger"),
    function(x, name) {
        if (!name %in% colnames(cellMeta(x))) NULL
        else cellMeta(x, columns = name)
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod("$", signature(x = "liger"),
                 function(x, name, value) {
                     cellMeta(x, columns = name) <- value
                     return(x)
                 })
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
setMethod(
    "defaultCluster",
    signature = c(x = "liger", useDatasets = "ANY"),
    function(x, useDatasets = NULL, droplevels = FALSE, ...) {
        # No name given, retrieve default
        useDatasets <- .checkUseDatasets(x, useDatasets)
        name <- x@uns$defaultCluster
        if (is.null(name)) return(NULL)
        cluster <- cellMeta(x, name, x$dataset %in% useDatasets)
        names(cluster) <- colnames(x)[x$dataset %in% useDatasets]
        if (isTRUE(droplevels)) cluster <- droplevels(cluster)
        return(cluster)
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "defaultCluster",
    signature(x = "liger", value = "character"),
    function(x, name = NULL, useDatasets = NULL, ..., value) {
        useDatasets <- .checkUseDatasets(x, useDatasets)
        cellIdx <- x$dataset %in% useDatasets
        if (length(value) == 1) {
            # If doing defaultCluster(obj) <- "someName"
            if (!is.null(name)) {
                warning("Cannot have `name` when selecting a name with ",
                        "`value`.")
            }
            name <- value
            if (!name %in% colnames(cellMeta(x))) {
                stop("Selected name does not exist in `cellMeta(x)`")
            }
            x@uns$defaultCluster <- name
        } else {
            # If doing defaultCluster(obj) <- c(blah, blah, blah, ...)
            defaultCluster(x, name, useDatasets) <- factor(value)
        }
        return(x)
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "defaultCluster",
    signature(x = "liger", value = "factor"),
    function(x, name = NULL, useDatasets = NULL, droplevels = TRUE, ..., value) {
        if (isTRUE(droplevels)) value <- droplevels(value)
        useDatasets <- .checkUseDatasets(x, useDatasets)
        cellIdx <- x$dataset %in% useDatasets
        if (length(value) != sum(cellIdx)) {
            stop("Length of `value` does not match with the number of cells")
        }
        if (is.null(name)) {
            .log("Storing given cluster labels to cellMeta(x) field: ",
                 "\"defaultCluster\"")
            name <- "defaultCluster"
        }
        if (is.null(names(value))) names(value) <- colnames(x)[cellIdx]
        cellMeta(x, name, cellIdx) <- value
        x@uns$defaultCluster <- name
        return(x)
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "defaultCluster",
    signature(x = "liger", value = "NULL"),
    function(x, name = NULL, useDatasets = NULL, ..., value) {
        x@uns$defaultCluster <- NULL
        if ("defaultCluster" %in% colnames(cellMeta(x))) {
            cellMeta(x, "defaultCluster") <- NULL
        }
        return(x)
    }
)

#' @export
#' @rdname liger-class
#' @section Dimension reduction access:
#' Currently, low-dimensional representaion of cells, presented as dense
#' matrices, are all stored in \code{cellMeta} slot, and can totally be accessed
#' with generics \code{cellMeta} and \code{cellMeta<-}. In addition to that,
#' we provide specific generics \code{dimRed} and \code{dimRed<-} for getting
#' and setting matrix like cell metadata, respectively. Adding a matrix to the
#' object looks as simple as \code{dimRed(obj, "name") <- matrixLike}. It can
#' be retrived back with \code{dimRed(obj, "name")}. Similar to having a default
#' cluster labeling, we also constructed the feature of default dimRed. It can
#' be set with \code{defaultDimRed(obj) <- "existingMatLikeVar"} and the matrix
#' can be retrieved with \code{defaultDimRed(obj)}.
setGeneric(
    "dimRed",
    function(x, name = NULL, useDatasets = NULL, ...) {
        standardGeneric("dimRed")
    }
)

#' @export
#' @rdname liger-class
setGeneric(
    "dimRed<-",
    function(x, name = NULL, useDatasets = NULL, ..., value) {
        standardGeneric("dimRed<-")
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "dimRed",
    signature = c(x = "liger", name = "missing", useDatasets = "ANY"),
    function(x, name = NULL, useDatasets = NULL, ...) {
        # No name given, retrieve default
        useDatasets <- .checkUseDatasets(x, useDatasets)
        name <- x@uns$defaultDimRed
        dimred <- NULL
        if (is.null(name)) {
            for (i in seq_along(cellMeta(x))) {
                if (!is.null(dim(cellMeta(x)[[i]]))) {
                    warning("No default dimRed recorded. Returning the first ",
                            "matrix like object in cellMeta(object)")
                    dimred <- cellMeta(x)[[i]]
                    break
                }
            }
            if (is.null(dimred)) {
                stop("No possible dimRed can be found in this liger object.")
            }
        } else {
            dimred <- cellMeta(x, name, x$dataset %in% useDatasets)
        }
        dimred <- as.matrix(dimred)
        rownames(dimred) <- colnames(x)[x$dataset %in% useDatasets]
        colnames(dimred) <- paste0(name, "_", seq_len(ncol(dimred)))
        return(dimred)
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "dimRed",
    signature = c(x = "liger", name = "character", useDatasets = "ANY"),
    function(x, name, useDatasets = NULL, ...) {
        # No name given, retrieve default
        useDatasets <- .checkUseDatasets(x, useDatasets)
        dimred <- cellMeta(x, name, x$dataset %in% useDatasets)
        if (is.null(dim(dimred))) {
            stop("Retrieved data for \"", name, "\" is not a matrix.")
        }
        dimred <- as.matrix(dimred)
        rownames(dimred) <- colnames(x)[x$dataset %in% useDatasets]
        colnames(dimred) <- paste0(name, "_", seq_len(ncol(dimred)))
        return(dimred)
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "dimRed",
    signature(x = "liger", name = "character", value = "matrixLike"),
    function(x, name = NULL, useDatasets = NULL, asDefault = NULL, ..., value) {
        useDatasets <- .checkUseDatasets(x, useDatasets)
        cellIdx <- x$dataset %in% useDatasets
        value <- as.matrix(value)
        if (is.null(asDefault)) {
            if (!is.null(x@uns$defaultDimRed)) asDefault <- FALSE
            else asDefault <- TRUE
        }
        colnames(value) <- seq_len(ncol(value))
        rownames(value) <- colnames(x)[cellIdx]
        cellMeta(x, name, cellIdx) <- value
        if (isTRUE(asDefault)) defaultDimRed(x) <- name
        return(x)
    }
)

#' @export
#' @rdname liger-class
setGeneric(
    "defaultDimRed",
    function(x, useDatasets = NULL) {
        standardGeneric("defaultDimRed")
    }
)

#' @export
#' @rdname liger-class
setGeneric(
    "defaultDimRed<-",
    function(x, name, useDatasets = NULL, value) {
        standardGeneric("defaultDimRed<-")
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "defaultDimRed",
    signature(x = "liger", useDatasets = "ANY"),
    function(x, useDatasets = NULL) {
        name <- x@uns$defaultDimRed
        if (is.null(name)) return(NULL)
        else dimRed(x, name = name, useDatasets = useDatasets)
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "defaultDimRed",
    signature(x = "liger", name = "missing", value = "character"),
    function(x, name = NULL, useDatasets = NULL, value) {
        value <- value[1]
        dimred <- cellMeta(x, value)
        if (is.null(dim(dimred))) {
            stop("Specified variable is not matrix like.")
        }
        if (ncol(dimred) == 0) {
            stop("Cannot set unexisting variable as default dimRed.")
        }
        x@uns$defaultDimRed <- value
        return(x)
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "defaultDimRed",
    signature(x = "liger", name = "character", value = "matrixLike"),
    function(x, name, useDatasets = NULL, value) {
        useDatasets <- .checkUseDatasets(x, useDatasets)
        cellIdx <- x$dataset %in% useDatasets
        colnames(value) <- seq_len(ncol(value))
        rownames(value) <- colnames(x)[cellIdx]
        cellMeta(x, name, cellIdx) <- value
        x@uns$defaultDimRed <- name
        return(x)
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
setMethod("varFeatures", signature(x = "liger"),
          function(x) x@varFeatures)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "varFeatures",
    signature(x = "liger", check = "ANY", value = "character"),
    function(x, check = TRUE, value) {
        x@varFeatures <- value
        for (d in names(x)) {
            ld <- dataset(x, d)
            featureMeta(ld, check = FALSE)$selected <- rownames(ld) %in% value
            datasets(x, check = FALSE)[[d]] <- ld
        }
        if (isTRUE(check)) {
            checkResult <- unlist(lapply(datasets(x), function(ld) {
                all(value %in% rownames(ld))
            }), use.names = FALSE)
            if (!all(checkResult)) {
                problem <- names(x)[!checkResult]
                warning("Not all variable features passed are ",
                        "found in datasets: ",
                        paste(problem, collapse = ", "))
            }
        }
        x
    }
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# S3 methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname liger-class
#' @export
#' @method fortify liger
fortify.liger <- function(model, data, ...) {
    df <- cellMeta(model, as.data.frame = TRUE)
    if (!is.null(model@H.norm)) df <- cbind(df, model@H.norm)
    df
}

#' @section Combining multiple liger object: The list of \code{datasets} slot,
#' the rows of \code{cellMeta} slot and the list of \code{commands} slot will
#' be simply concatenated. Variable features in \code{varFeatures} slot will be
#' taken a union. The \eqn{W} and \eqn{H.norm} matrices are not taken into
#' account for now.
#' @rdname liger-class
#' @export
#' @method c liger
c.liger <- function(...) {
    objList <- list(...)
    if (any(sapply(objList, function(obj) !inherits(obj, "liger"))))
        stop("Can only combine `liger` objects with `c(...)` method for now.")
    objList[[length(objList)]] <- recordCommand(objList[[length(objList)]])
    allDatasets <- list()
    allCellMeta <- NULL
    varFeatures <- character()
    allCommands <- list()
    for (i in seq_along(objList)) {
        obj <- objList[[i]]
        allDatasets <- c(allDatasets, datasets(obj))
        allCellMeta <- rbind(allCellMeta, obj@cellMeta)
        varFeatures <- union(varFeatures, varFeatures(obj))
        allCommands <- c(allCommands, obj@commands)
    }
    methods::new("liger", datasets = allDatasets, cellMeta = allCellMeta,
                 varFeatures = varFeatures, commands = allCommands,
                 version = utils::packageVersion("rliger2"))
}
