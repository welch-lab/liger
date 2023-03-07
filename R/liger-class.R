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
#' @description Data container for LIGER analysis. The slot \code{datasets} is a
#' list where each element should be \linkS4class{ligerDataset} containing
#' dataset specific information, such as the expression matrices. The other
#' parts of liger object stores information that can be shared across the
#' analysis, such as the cell metadata and factorization result matrices.
#' @details
#' Accessor methods are defined for getting and setting values from or
#' to the object, respectively.
#'
#' A liger object with version number >= 1.2.0 is fortify-able which means users
#' can produce ggplots with the cell metadata simply by calling
#' \code{ggplot(liger_object)}.
#' @slot datasets list of \linkS4class{ligerDataset} objects.
#' @slot cellMeta Cell metadata, \linkS4class{DFrame}
#' @slot varFeatures vector of feature names
#' @slot H.norm matrix
#' @slot H matrix
#' @slot commands Record of analysis
#' @slot uns Unstructured meta-info of analyses
#' @slot version Record of version of rliger package
#' @importClassesFrom S4Vectors DataFrame
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
        version = utils::packageVersion("rliger")
    )
)

# ------------------------------------------------------------------------------
# Object constructor ####
# ------------------------------------------------------------------------------

#' Create liger object
#' @description This function allows creating \linkS4class{liger} object from
#' multiple datasets of various forms (See \code{rawData}).
#' @param rawData Named list of datasets. Required. Elements allowed include a
#' matrix, a \code{Seurat} object, a \code{SingleCellExperiment} object, an
#' \code{AnnData} object, a \linkS4class{ligerDataset} object or a filename to
#' an HDF5 file. See detail for HDF5 reading.
#' @param modal Character vector for modality setting. Currently options of
#' \code{"default"}, \code{"rna"}, and \code{"atac"} are supported.
#' @param cellMeta data.frame of metadata at single-cell level. Default
#' \code{NULL}.
#' @param removeMissing Logical. Whether to remove cells that do not have any
#' counts and features not expressed in any cells from each dataset. Default
#' \code{TRUE}. H5 based dataset with less than 8000 cells will be subset into
#' memory.
#' @param formatType Select preset of H5 file structure. Current available
#' options are \code{"10X"} and \code{"AnnData"}. Can be either a single
#' specification for all datasets or a character vector that match with each
#' dataset.
#' @param dataName,indicesName,indptrName The path in a H5 file for the raw
#' sparse matrix data. These three types of data stands for the \code{x},
#' \code{i}, and \code{p} slots of a \code{\link[Matrix]{dgCMatrix-class}}
#' object. Default \code{NULL} uses \code{formatType} preset.
#' @param genesName,barcodesName The path in a H5 file for the gene names and
#' cell barcodes. Default \code{NULL} uses \code{formatType} preset.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{TRUE}.
#' @param remove.missing,format.type,data.name,indices.name,indptr.name,genes.name,barcodes.name
#' \bold{Deprecated.} See Usage section for replacement.
#' @export
#' @seealso \code{\link{createLigerDataset}}, \code{\link{createH5LigerDataset}}
createLiger <- function(
        rawData,
        modal = NULL,
        cellMeta = NULL,
        removeMissing = TRUE,
        formatType = "10X",
        dataName = NULL,
        indicesName = NULL,
        indptrName = NULL,
        genesName = NULL,
        barcodesName = NULL,
        verbose = TRUE,
        # Deprecated coding style
        remove.missing = removeMissing,
        format.type = formatType,
        data.name = dataName,
        indices.name = indicesName,
        indptr.name = indptrName,
        genes.name = genesName,
        barcodes.name = barcodesName
) {
    .deprecateArgs(list(remove.missing = "removeMissing",
                        format.type = "formatType", data.name = "dataName",
                        indices.name = "indicesName",
                        indptr.name = "indptrName", genes.name = "genesName",
                        barcodes.name = "barcodesName"))
    if (!is.list(rawData)) stop("`rawData` has to be a named list.")

    nData <- length(rawData)
    if (missing(modal) || is.null(modal)) modal <- "default"
    modal <- tolower(modal)
    if (length(modal) == 1) modal <- rep(modal, nData)
    else if (length(modal) != nData)
        stop("Wrong length of `modal`. ",
             "Specify only 1 or match the length of `datasets`. ",
             "See ?createLiger for valid options.")
    # TODO handle h5 specific argument for hybrid of H5 and in memory stuff.
    datasets <- list()
    for (i in seq_along(rawData)) {
        dname <- names(rawData)[i]
        data <- rawData[[i]]
        if (is.character(data)) {
            # Assuming character input is a filename
            datasets[[dname]] <- createH5LigerDataset(
                h5file = data,
                formatType = formatType,
                rawData = dataName,
                barcodesName = barcodesName,
                genesName = genesName,
                indicesName = indicesName,
                indptrName = indptrName,
                modal = modal[i]
            )
        } else {
            datasets[[dname]] <- as.ligerDataset(data, modal = modal[i])
        }
    }

    datasets <- .dedupLigerDatasets(datasets)
    barcodes <- unlist(lapply(datasets, colnames), use.names = FALSE)
    if (is.null(cellMeta)) {
        cellMeta <- S4Vectors::DataFrame(
            dataset = factor(rep(names(datasets), lapply(datasets, ncol))),
            row.names = barcodes)
    } else {
        cellMeta <- S4Vectors::DataFrame(cellMeta)
        cellMeta <- cellMeta[barcodes,,drop = FALSE]
        # Force writing `dataset` variable as named by @datasets
        cellMeta$dataset <- factor(rep(names(datasets),
                                        lapply(datasets, ncol)))
    }
    obj <- methods::new("liger",
                        datasets = datasets,
                        cellMeta = cellMeta)
    obj <- runGeneralQC(obj, verbose = verbose)
    if (isTRUE(removeMissing)) {
        obj <- removeMissing(obj, "both", filenameSuffix = "qc",
                             verbose = verbose)
    }

    return(obj)
}

#' Deduplicate barcodes from all datasets
#' @param  datasets list of ligerDataset object
#' @noRd
.dedupLigerDatasets <- function(datasets) {
    allBarcodes <- unlist(lapply(datasets, colnames), use.names = FALSE)
    if (any(duplicated(allBarcodes))) {
        datasetVar <- rep(names(datasets), lapply(datasets, ncol))
        # Main deduplication process. Wondering if make stand-alone dot-func
        dups <- unique(allBarcodes[duplicated(allBarcodes)])
        for (bc in dups) {
            idx <- which(allBarcodes == bc)
            allBarcodes[idx] <- paste0(bc, "-", seq_along(idx))
        }
        # Done. Then assign them back to ligerDataset objects
        for (d in names(datasets)) {
            colnames(datasets[[d]]) <- allBarcodes[datasetVar == d]
        }
    }
    return(datasets)
}

# Expand rows of raw matrix in each dataset to the union gene list of all
# datasets
.takeGeneUnion.rawData <- function(rawData) {

}
# ------------------------------------------------------------------------------
# Validity ####
# ------------------------------------------------------------------------------

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
            if (!is.null(ld@V))
                if (!identical(rownames(ld@V), varFeatures(x)))
                    return(paste("Variable features do not match dimension",
                                 "of V matrix in dataset", d))
            if (!is.null(scaleData(ld))) {
                if (!isH5Liger(ld)) {
                    if (!identical(rownames(scaleData(ld)), varFeatures(x)))
                        return(paste("Variable features do not match dimension",
                                     "of scaleData in dataset", d))
                } else {
                    if (scaleData(ld)$dims[1] != length(varFeatures(x)))
                        return(paste("Variable features do not match dimension",
                                     "of scaleData in dataset (H5)", d))
                }
            }

        }
    }
    return(NULL)
}

.valid.liger <- function(object) {
    message("Checking liger object validity")
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
#' 1.2.0. Otherwise \code{FALSE}
#' @noRd
is.newLiger <- function(object) {
    v <- as.character(object@version)
    result <- utils::compareVersion("1.2.0", v)
    if (result == 0 | result == -1) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# ------------------------------------------------------------------------------
# Generics ####
# ------------------------------------------------------------------------------
#' Access data in a liger object
#' @description The methods listed in this manual only apply to
#' \linkS4class{liger} objects, while there are other generics that have methods
#' for both \code{liger} and \linkS4class{ligerDataset} objects. Please see
#' \code{\link{rawData}} page.
#' @param x \linkS4class{liger} object
#' @param name Name or numeric index of a dataset
#' @param value Check \linkS4class{liger} class specification for required value
#' type.
#' @param type When using \code{dataset<-} with a matrix like \code{value},
#' specify what type the matrix is. Choose from \code{"rawData"},
#' \code{"normData"} or \code{"scaleData"}.
#' @param qc Logical, whether to perform general qc on added new dataset.
#' @param check Logical, whether to perform object validity check on setting new
#' value.
#' @param columns The names of available variables in \code{cellMeta} slot. When
#' \code{as.data.frame = TRUE}, please use variable names after coercion.
#' @param cellIdx Valid cell subscription to subset retrieved variables. Default
#' \code{NULL} uses all cells.
#' @param as.data.frame Logical, whether to apply
#' \code{\link[base]{as.data.frame}} on the subscription. Default \code{FALSE}.
#' @param i,j Feature or cell index for \code{`[`} method. For \code{`[[`}
#' method, use a single variable name with \code{i} and \code{j} is not
#' applicable.
#' @param ... For \code{cellMeta}, arguments passed to \code{[]} subscriptor,
#' e.g. \code{drop = FALSE}.
#' @export
#' @rdname liger-access
setGeneric("datasets", function(x, check = NULL) standardGeneric("datasets"))

#' @export
#' @rdname liger-access
setGeneric("datasets<-", function(x, check = TRUE, value) standardGeneric("datasets<-"))

#' @export
#' @rdname liger-access
setGeneric("dataset", function(x, name = NULL) standardGeneric("dataset"))

#' @export
#' @rdname liger-access
setGeneric("dataset<-", function(x, name, type = NULL, qc = TRUE, value) {
    standardGeneric("dataset<-")
})

#' @export
#' @rdname liger-access
setGeneric("cellMeta", function(x, columns = NULL, cellIdx = NULL, as.data.frame = FALSE, ...) standardGeneric("cellMeta"))

#' @export
#' @rdname liger-access
setGeneric("cellMeta<-", function(x, columns = NULL, check = FALSE, value) standardGeneric("cellMeta<-"))

#' @export
#' @rdname liger-access
setGeneric("varFeatures", function(x) standardGeneric("varFeatures"))

#' @export
#' @rdname liger-access
setGeneric("varFeatures<-", function(x, check = TRUE, value) standardGeneric("varFeatures<-"))

# ------------------------------------------------------------------------------
# S4 methods ####
# ------------------------------------------------------------------------------

#' show method for liger
#'
#' @param object liger object
#' @name show
#' @aliases show,liger-method
#' @docType methods
#' @rdname show-methods
setMethod(
    f = "show",
    signature(object = "liger"),
    definition = function(object) {
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

#' Access dataset names or length
#' @param x \linkS4class{liger} object
#' @param value Character vector of new names for datasets.
#' @return character vector of dataset names or number of datasets
#' @rdname names-liger
#' @export
setMethod("names", signature(x = "liger"), function(x) {
    names(datasets(x))
})

#' @rdname names-liger
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
            x@cellMeta$dataset <- factor(x@cellMeta$dataset)
            names(x@datasets) <- value
        }
        x
    })

#' @rdname names-liger
#' @export
setMethod("length", signature(x = "liger"), function(x) {
    length(datasets(x))
})

#' Access dimensionality of liger object
#' @description For a \linkS4class{liger} object or a
#' \linkS4class{ligerDataset} object, the column orientation is
#' assigned for cells while rows should be for features/genes. However, due to
#' the data structure, it is hard to define a row index for the outer container
#' -- \code{liger} object, which might contain datasets that vary in genes.
#'
#' Therefore, for \code{liger} objects, \code{dim} and \code{dimnames} returns
#' \code{NA}/\code{NULL} for rows and total cell counts/barcodes for the
#' columns. For \code{ligerDataset} objects, gene index will be available for
#' rows.
#' @param x \linkS4class{liger} object
#' @param value For direct call of \code{dimnames<-} method, should be a list
#' with \code{NULL} as the first element and valid cell identifiers as the
#' second element. For \code{colnames<-} method, the character vector of cell
#' identifiers. \code{rownames<-} method is not applicable.
#' @return See Description
#' @rdname dim-liger
#' @export
setMethod("dim", "liger", function(x) {
    c(NA, nrow(cellMeta(x)))
})

# `dimnames()` method set so that `rownames()` returns `NULL`, `colnames()`
# returns all barcodes.
#' @rdname dim-liger
#' @export
setMethod("dimnames", "liger", function(x) {
    list(NULL, rownames(cellMeta(x)))
})

#' @rdname dim-liger
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
    rownames(x@W) <- value[[1L]]
    if (!is.null(x@H.norm)) colnames(x@H.norm) <- value[[2L]]
    x
})

#' @rdname liger-access
#' @export
setMethod(
    "[[",
    c("liger", "ANY", "missing"),
    function(x, i, j, ...) {
        if (is.null(i)) return(NULL)
            cellMeta(x, columns = i, ...)
    }
)

#' @rdname liger-access
#' @export
setReplaceMethod(
    "[[",
    c("liger", "ANY", "missing"),
    function(x, i, j, ..., value) {
        cellMeta(x, columns = i, ...) <- value
        x
    }
)

.DollarNames.liger <- function(x, pattern = "")
    grep(pattern, colnames(x@cellMeta), value = TRUE)

#' @export
#' @rdname liger-access
setMethod("$", signature(x = "liger"),
          function(x, name) {
              cellMeta(x, columns = name)
          })

#' @export
#' @rdname liger-access
setReplaceMethod("$", signature(x = "liger"),
                 function(x, name, value) {
                     cellMeta(x, columns = name) <- value
                     return(x)
                 })

#' @export
#' @rdname subsetLiger
setMethod(
    "[",
    signature(x = "liger", i = "character", j = "missing"),
    function(x, i, j, ...) subsetLiger(x, featureIdx = i, cellIdx = NULL, ...)
)

#' @export
#' @rdname subsetLiger
setMethod(
    "[",
    signature(x = "liger", i = "missing", j = "index"),
    function(x, i, j, ...) subsetLiger(x, featureIdx = NULL, cellIdx = j, ...)
)

#' @export
#' @rdname subsetLiger
setMethod(
    "[",
    signature(x = "liger", i = "character", j = "index"),
    function(x, i, j, ...) subsetLiger(x, featureIdx = i, cellIdx = j, ...)
)

#' @section Whether to use datasets() or dataset()?:
#' \code{datasets()} method only accesses the \code{datasets} slot, the list of
#' \linkS4class{ligerDataset} objects. \code{dataset()} method accesses a single
#' dataset, with subsequent cell metadata updates and checks bonded when adding
#' or modifying a dataset. Therefore, when users want to modify something inside
#' a \code{ligerDataset} while no cell metadata change should happen, it is
#' recommended to use: \code{datasets(x)[[name]] <- ligerD} for efficiency,
#' though the result would be the same as \code{dataset(x, name) <- ligerD}.
#' @export
#' @rdname liger-access
setMethod("datasets", signature(x = "liger", check = "ANY"),
          function(x, check = NULL) x@datasets)

#' @export
#' @rdname liger-access
setReplaceMethod("datasets", signature(x = "liger", check = "logical"),
                 function(x, check = TRUE, value) {
                     x@datasets <- value
                     if (isTRUE(check)) methods::validObject(x)
                     x
                 })

#' @export
#' @rdname liger-access
setReplaceMethod("datasets", signature(x = "liger", check = "missing"),
                 function(x, check = TRUE, value) {
                     x@datasets <- value
                     methods::validObject(x)
                     x
                 })

#' @export
#' @rdname liger-access
setMethod("dataset", signature(x = "liger", name = "character_OR_NULL"),
          function(x, name = NULL) {
              if (is.null(name)) return(datasets(x)[[1]])
              else {
                  if (!name %in% names(x)) {
                      stop('Specified dataset name "', name,
                           '" not found in liger object.')
                  }
                  return(datasets(x)[[name]])
              }
          })

#' @export
#' @rdname liger-access
setMethod("dataset", signature(x = "liger", name = "missing"),
          function(x, name = NULL) {
              datasets(x)[[1]]
          })

#' @export
#' @rdname liger-access
setMethod("dataset", signature(x = "liger", name = "numeric"),
          function(x, name = NULL) {
              datasets(x)[[name]]
          })

#' @export
#' @rdname liger-access
setReplaceMethod("dataset", signature(x = "liger", name = "character",
                                      type = "missing", qc = "ANY",
                                      value = "ligerDataset"),
                 function(x, name, type = NULL, qc = TRUE, value) {
                     if (name %in% names(x)) {
                         dataset(x, name) <- NULL
                     }
                     new.idx <- seq(ncol(x) + 1, ncol(x) + ncol(value))
                     x@datasets[[name]] <- value
                     # TODO also add rows to cellMeta and H.norm
                     x@cellMeta[new.idx, ] <- NA
                     rownames(x@cellMeta)[new.idx] <- colnames(value)
                     levels(x@cellMeta$dataset) <-
                         c(levels(x@cellMeta$dataset), name)
                     x@cellMeta$dataset[new.idx] <- name
                     # x@W is genes x k, no need to worry
                     if (!is.null(x@H.norm)) {
                         message("Filling in NAs to H.norm matrix")
                         H.normNew <- matrix(
                             NA, ncol(value), ncol(x@H.norm),
                             dimnames = list(colnames(value), NULL))
                         x@H.norm <- rbind(x@H.norm, H.normNew)
                     }
                     methods::validObject(x)
                     if (qc) x <- runGeneralQC(x, useDatasets = name)
                     x
                 })

#' @export
#' @rdname liger-access
setReplaceMethod("dataset", signature(x = "liger", name = "character",
                                      type = "ANY", qc = "ANY",
                                      value = "matrixLike"),
                 function(x, name,
                          type = c("rawData", "normData", "scaleData"),
                          qc = TRUE,
                          value) {
                     type <- match.arg(type)
                     if (type == "rawData") {
                         ld <- createLigerDataset(rawData = value)
                     } else if (type == "normData") {
                         ld <- createLigerDataset(normData = value)
                     } else if (type == "scaleData") {
                         ld <- createLigerDataset(scaleData = value)
                     }
                     dataset(x, name, qc = qc) <- ld
                     x
                 })

#' @export
#' @rdname liger-access
setReplaceMethod("dataset", signature(x = "liger", name = "character",
                                      type = "missing", qc = "ANY",
                                      value = "NULL"),
                 function(x, name, type = NULL, qc = TRUE, value) {
                     if (!name %in% names(x)) {
                         warning("Specified dataset name not found in ",
                                 "liger object. Nothing would happen.")
                     } else {
                         idxToRemove <- x$dataset == name
                         x@datasets[[name]] <- NULL
                         x@cellMeta <- x@cellMeta[!idxToRemove, , drop = FALSE]
                         x@H.norm <- x@H.norm[!idxToRemove, , drop = FALSE]
                     }
                     x
                 })

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
        if (is.vector(res)) res <- res[cellIdx]
        else if (!is.null(dim(res))) res <- res[cellIdx, , ...]
        else stop("Result before idx subscription corrupted")
    }
    return(res)
}

#' @export
#' @rdname liger-access
#' @section Cell metadata access:
#' Three approaches are provided for access of cell metadata. A generic function
#' \code{cellMeta} is implemented with plenty of options and multi-variable
#' accessibility. Besides, users can use double-bracket (e.g.
#' \code{ligerObj[[Var]]}) or dollor-sign (e.g. \code{ligerObj$nUMI}) to access
#' or modify single variables.
setMethod(
    "cellMeta",
    signature(x = "liger", columns = "NULL"),
    function(x, columns = NULL, cellIdx = NULL, as.data.frame = FALSE, ...)
        NULL
)

#' @export
#' @rdname liger-access
setMethod(
    "cellMeta",
    signature(x = "liger", columns = "character"),
    function(x, columns = NULL, cellIdx = NULL, as.data.frame = FALSE, ...)
        .subsetCellMeta(x, columns = columns, cellIdx = cellIdx,
                        as.data.frame = as.data.frame, ...)
)

#' @export
#' @rdname liger-access
setMethod(
    "cellMeta",
    signature(x = "liger", columns = "missing"),
    function(x, columns = NULL, cellIdx = NULL, as.data.frame = FALSE, ...)
        .subsetCellMeta(x, columns = NULL, cellIdx = cellIdx,
                        as.data.frame = as.data.frame, ...)
)

#' @export
#' @rdname liger-access
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
#' @rdname liger-access
setReplaceMethod(
    "cellMeta",
    signature(x = "liger", columns = "character"),
    function(x, columns = NULL, check = FALSE, value) {
        x@cellMeta[[columns]] <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-access
setMethod("varFeatures", signature(x = "liger"),
          function(x) x@varFeatures)

#' @export
#' @rdname liger-access
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

# ------------------------------------------------------------------------------
# S3 methods ####
# ------------------------------------------------------------------------------

#' Fortify method for classes from the liger package
#'
#' @param x \code{liger} object to convert into a dataframe.
#' @return A \code{data.frame} combining \code{cellMeta(x)} with \code{H.norm}
#' matrix.
#' @name fortify.liger
#' @importFrom ggplot2 fortify
NULL

#' @rdname fortify.liger
#' @export
#' @method fortify liger
fortify.liger <- function(x) {
    df <- cellMeta(x, as.data.frame = TRUE)
    if (!is.null(x@H.norm)) df <- cbind(df, x@H.norm)
    df
}

#' Combine multiple liger object in to one
#' @description The list of \code{datasets} slot, the rows of
#' \code{cellMeta} slot and the list of \code{commands} slot will be simply
#' concatenated. Variable features in \code{varFeatures} slot will be taken a
#' union. The \eqn{W} and \eqn{H.norm} matrices are not taken into account for
#' now.
#' @param ... \linkS4class{liger} objects
#' @return A new \linkS4class{liger} object containing all datasets in the order
#' of specification.
#' @name c.liger
NULL

#' @rdname c.liger
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
                 version = utils::packageVersion("rliger"))
}
