setClassUnion("dgCMatrix_OR_NULL", c("dgCMatrix", "NULL"))
setClassUnion("matrix_OR_NULL", c("matrix", "NULL"))
setClassUnion("matrixLike_OR_NULL", c(
    "matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix", "NULL"
))
setClassUnion("character_OR_NULL", c("character", "NULL"))
setClassUnion("matrixLike", c(
    "matrix", "dgCMatrix", "dgTMatrix", "dgeMatrix"
))

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
#' @slot cell.meta Cell metadata, \linkS4class{DFrame}
#' @slot var.features vector of feature names
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
        cell.meta = "DataFrame",
        var.features = "character_OR_NULL",
        W = "matrix_OR_NULL",
        H.norm = "matrix_OR_NULL",
        commands = "list",
        uns = "list",
        version = "ANY"
    ),
    methods::prototype(
        cell.meta = new("DFrame"),
        version = packageVersion("rliger")
    )
)

# ------------------------------------------------------------------------------
# Object constructor ####
# ------------------------------------------------------------------------------

#' Create liger object
#' @description This function allows creating \linkS4class{liger} object from
#' multiple datasets of various forms (See \code{raw.data}).
#' @param raw.data Named list of datasets. Required. Elements allowed include a
#' matrix, a \linkS4class{Seurat} object, a \linkS4class{SingleCellExperiment}
#' object, an \code{AnnData} object, a \linkS4class{ligerDataset} object or a
#' filename to an HDF5 file. See detail for HDF5 reading.
#' @param modal Character vector for modality setting. See detail.
#' @param cell.meta data.frame of metadata at single-cell level. Default
#' @param remove.missing Whether to remove cells that do not have any counts and
#' features not expressed in any cells from each dataset. Default \code{TRUE}.
#' @param format.type Select preset of H5 file structure. Current available
#' options are \code{"10X"} and \code{"AnnData"}. Can be either a single
#' specification for all datasets or a character vector that match with each
#' dataset.
#' @export
#' @seealso \code{\link{createLigerDataset}},
#' \code{\link{createLigerDataset.h5}}
createLiger <- function(
        raw.data,
        modal = NULL,
        cell.meta = NULL,
        remove.missing = TRUE,
        format.type = "10X",
        data.name = NULL,
        indices.name = NULL,
        indptr.name = NULL,
        genes.name = NULL,
        barcodes.name = NULL,
        verbose = TRUE
) {
    if (!is.list(raw.data)) stop("`raw.data` has to be a named list.")

    nData <- length(raw.data)
    if (missing(modal) || is.null(modal)) modal <- "default"
    modal <- tolower(modal)
    if (length(modal) == 1) modal <- rep(modal, nData)
    else if (length(modal) != nData)
        stop("Wrong length of `modal`. ",
             "Specify only 1 or match the length of `datasets`. ",
             "See ?createLiger for valid options.")
    # TODO handle h5 specific argument for hybrid of H5 and in memory stuff.
    datasets <- list()
    for (i in seq_along(raw.data)) {
        dname <- names(raw.data)[i]
        data <- raw.data[[i]]
        if (is.character(data)) {
            # Assuming character input is a filename
            datasets[[dname]] <- createLigerDataset.h5(
                h5file = data,
                format.type = format.type,
                raw.data = data.name,
                barcodes.name = barcodes.name,
                genes.name = genes.name,
                indices.name = indices.name,
                indptr.name = indptr.name,
                modal = modal[i]
            )
        } else {
            datasets[[dname]] <- as.ligerDataset(data, modal = modal[i])
        }
    }

    datasets <- .dedupLigerDatasets(datasets)
    barcodes <- unlist(lapply(datasets, colnames), use.names = FALSE)
    if (is.null(cell.meta)) {
        cell.meta <- S4Vectors::DataFrame(
            dataset = factor(rep(names(datasets), lapply(datasets, ncol))),
            row.names = barcodes)
    } else {
        cell.meta <- S4Vectors::DataFrame(cell.meta)
        cell.meta <- cell.meta[barcodes,]
        # Force writing `dataset` variable as named by @datasets
        cell.meta$dataset <- factor(rep(names(datasets),
                                        lapply(datasets, ncol)))
    }
    obj <- methods::new("liger",
                        datasets = datasets,
                        cell.meta = cell.meta)
    obj <- runGeneralQC(obj, verbose = verbose, )
    if (isTRUE(remove.missing)) {
        obj <- removeMissing(obj, "both", filenameSuffix = "qc",
                             verbose = verbose)
    }

    return(obj)
}

#' Deduplicate barcodes from all datasets
#'
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
.takeGeneUnion.rawData <- function(raw.data) {

}
# ------------------------------------------------------------------------------
# Validity ####
# ------------------------------------------------------------------------------

.checkLigerBarcodes <- function(x) {
    #colnames <- get("colnames", .GlobalEnv)
    #ncol <- get("ncol", .GlobalEnv)
    #cell.meta <- get("cell.meta", .GlobalEnv)
    #datasets <- get("datasets", .GlobalEnv)
    bcFromDatasets <- unlist(lapply(datasets(x), colnames), use.names = FALSE)
    if (!identical(colnames(x), bcFromDatasets)) {
        return("liger object barcodes do not match to barcodes in datasets")
    }
    if (!"dataset" %in% names(cell.meta(x))) {
        return("`datasets` variable missing in cell.meta(x)")
    }
    datasetNamesFromDatasets <- rep(names(x), lapply(datasets(x), ncol))
    names(datasetNamesFromDatasets) <- NULL

    if (!identical(datasetNamesFromDatasets, as.character(x$dataset))) {
        return("names of datasets do not match
               `datasets` variable in cell.meta")
    }
    NULL
}

.valid.liger <- function(object) {
    message("Checking liger object validity")
    .checkLigerBarcodes(object)
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

setGeneric("datasets", function(x, check = NULL, ...) standardGeneric("datasets"))
setGeneric("datasets<-", function(x, check = TRUE, ...) standardGeneric("datasets<-"))
setGeneric("dataset", function(x, name = NULL) standardGeneric("dataset"))
setGeneric("dataset<-", function(x, name, type = NULL, qc = TRUE, value) {
    standardGeneric("dataset<-")
})
setGeneric("cell.meta", function(x, ...) standardGeneric("cell.meta"))
setGeneric("cell.meta<-", function(x, ...) standardGeneric("cell.meta<-"))
setGeneric("var.features", function(x) standardGeneric("var.features"))
setGeneric("var.features<-", function(x, check = TRUE, value) standardGeneric("var.features<-"))

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
        cat(paste0("cell.meta(", ncol(cell.meta(object)), "): "))
        cat(.collapseLongNames(colnames(cell.meta(object))), "\n")
        invisible(x = NULL)
    }
)

setMethod("names", signature(x = "liger"), function(x) {
    names(datasets(x))
})

setReplaceMethod(
    "names",
    signature(x = "liger", value = "character"),
    function(x, value) {
        originalNames <- names(x)
        if (!identical(value, originalNames)) {
            dataset.idx <- lapply(originalNames, function(n) {
                x$dataset == n
            })
            x@cell.meta$dataset <- as.character(x@cell.meta$dataset)
            for (i in seq_along(value)) {
                x@cell.meta$dataset[dataset.idx[[i]]] <- value[i]
            }
            x@cell.meta$dataset <- factor(x@cell.meta$dataset)
            names(x@datasets) <- value
        }
        x
    })

setMethod("length", signature(x = "liger"), function(x) {
    length(datasets(x))
})

# `dim()` method set so that `ncol()` returns total number of cells,
# `nrow()` returns `NA`
setMethod("dim", "liger", function(x) {
    c(NA, nrow(cell.meta(x)))
})

# `dimnames()` method set so that `rownames()` returns `NULL`, `colnames()`
# returns all barcodes.
setMethod("dimnames", "liger", function(x) {
    list(NULL, rownames(cell.meta(x)))
})

setReplaceMethod("dimnames", c("liger", "list"), function(x, value) {
    dataset.var <- x$dataset
    dataset.uniq <- unique(dataset.var)
    for (d in dataset.uniq) {
        dataset.idx <- dataset.var == d
        # Not using `dataset()` accessor because it does object validity check
        colnames(x@datasets[[d]]) <- value[[2L]][dataset.idx]
    }
    rownames(cell.meta(x)) <- value[[2L]]
    colnames(x@W) <- value[[1L]]
    if (!is.null(x@H.norm)) colnames(x@H.norm) <- value[[2L]]
    x
})

setMethod("[[", signature(x = "liger", i = "ANY", j = "missing"),
          function(x, i, j, ...)
          {
              cell.meta(x)[[i, ...]]
          })

setReplaceMethod("[[", signature(x = "liger", i = "ANY", j = "missing"),
                 function(x, i, j, ..., value)
                 {
                     x@cell.meta[[i, ...]] <- value
                     x
                 })

.DollarNames.liger <- function(x, pattern = "")
    grep(pattern, names(cell.meta(x)), value = TRUE)

setMethod("$", signature(x = "liger"),
          function(x, name) {
              cell.meta(x)[[name]]
          })

setReplaceMethod("$", signature(x = "liger"),
                 function(x, name, value) {
                     cell.meta(x)[[name]] <- value
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

#' Access datasets of liger object
#' @description This method access the whole slot \code{datasets}, a list, of
#' the \linkS4class{liger} object. So that operations like
#' \code{datasets(x)[[1]] <- y} modifies the inner \code{ligerDataset} object
#' without triggering cell metadata update, which will be done by
#' \code{dataset(x) <- y}.
#' @param x liger object
#' @param check Whether to perform object validity check when using setter
#' method. Default \code{TRUE}.
#' @param value list of \linkS4class{ligerDataset} object.
#' @seealso \code{\link{dataset}}
#' @export
#' @rdname datasets
setMethod("datasets", signature(x = "liger", check = "ANY"),
          function(x, check = NULL) x@datasets)

#' @export
#' @rdname datasets
setReplaceMethod("datasets", signature(x = "liger", check = "logical"),
                 function(x, check = TRUE, value) {
                     x@datasets <- value
                     if (isTRUE(check)) validObject(x)
                     x
                 })

#' @export
#' @rdname datasets
setReplaceMethod("datasets", signature(x = "liger", check = "missing"),
                 function(x, check = TRUE, value) {
                     x@datasets <- value
                     validObject(x)
                     x
                 })

#' @title Access single dataset of liger object
#' @description TODO: Debating on whether to have a separate (differently named)
#' method to directly access matrices in embedded ligerDataset objects, or to
#' allow accessing matrices with \code{type}.
#'
#' This method access a single dataset of the \linkS4class{liger} object. Using
#' the setter method to add/modify a dataset will subsequentially perform cell
#' metadata updates and checks.
#' @rdname dataset
#' @param x \linkS4class{liger} object
#' @param name name of dataset
#' @param type Type of specified matrix. Only works when \code{value} is a
#' matrix like object. choose from \code{"raw.data"}, \code{"norm.data"} or
#' \code{"scale.data"}. Default \code{"raw.data"}.
#' @param value A matrix of a dataset to be added, or a constructed
#' \linkS4class{ligerDataset} object. \code{NULL} to remove a dataset by
#' \code{name}.
#' @seealso \code{\link{datasets}}
#' @return Getter methods returns a \linkS4class{ligerDataset} object, while
#' setter methods update the input \code{x} (the \linkS4class{liger} object).
#' @export
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

#' @rdname dataset
#' @export
setMethod("dataset", signature(x = "liger", name = "missing"),
          function(x, name = NULL) {
              datasets(x)[[1]]
          })

#' @rdname dataset
#' @export
setMethod("dataset", signature(x = "liger", name = "numeric"),
          function(x, name = NULL) {
              datasets(x)[[name]]
          })

#' @rdname dataset
#' @export
setReplaceMethod("dataset", signature(x = "liger", name = "character",
                                      type = "missing", qc = "ANY",
                                      value = "ligerDataset"),
                 function(x, name, type = NULL, qc = TRUE, value) {
                     if (name %in% names(x)) {
                         dataset(x, name) <- NULL
                     }
                     new.idx <- seq(ncol(x) + 1, ncol(x) + ncol(value))
                     x@datasets[[name]] <- value
                     # TODO also add cols to H/H.norm and rows to cell.meta
                     x@cell.meta[new.idx, ] <- NA
                     rownames(x@cell.meta)[new.idx] <- colnames(value)
                     levels(x@cell.meta$dataset) <-
                         c(levels(x@cell.meta$dataset), name)
                     x@cell.meta$dataset[new.idx] <- name
                     # x@W is k x genes, no need to worry
                     if (!is.null(x@H.norm)) {
                         message("Filling in NAs to H.norm matrix")
                         x@H.norm[, new.idx] <- NA
                         colnames(x@H.norm)[new.idx] <- colnames(value)
                     }
                     validObject(x)
                     if (qc) x <- runGeneralQC(x, useDatasets = name)
                     x
                 })

#' @rdname dataset
#' @export
setReplaceMethod("dataset", signature(x = "liger", name = "character",
                                      type = "ANY", qc = "ANY",
                                      value = "matrixLike"),
                 function(x, name,
                          type = c("raw.data", "norm.data", "scale.data"),
                          qc = TRUE,
                          value) {
                     type <- match.arg(type)
                     if (type == "raw.data") {
                         ld <- ligerDataset(raw.data = value)
                     } else if (type == "norm.data") {
                         ld <- ligerDataset(norm.data = value)
                     } else if (type == "scale.data") {
                         ld <- ligerDataset(scale.data = value)
                     }
                     dataset(x, name, qc = qc) <- ld
                     x
                 })

#' @rdname dataset
#' @export
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
                         x@cell.meta <- x@cell.meta[!idxToRemove, , drop = FALSE]
                         x@H.norm <- x@H.norm[!idxToRemove,]
                     }
                     x
                 })

#' Access cell.meta of liger object
#' @param x \linkS4class{liger} object
#' @rdname cell.meta
#' @export
setMethod("cell.meta", signature(x = "liger"),
          function(x) x@cell.meta)

#' @rdname cell.meta
#' @export
setReplaceMethod("cell.meta", signature(x = "liger"),
                 function(x, value) {
                     if (!inherits(value, "DFrame"))
                         value <- S4Vectors::DataFrame(value)
                     x@cell.meta <- value
                     validObject(x)
                     x
                 })

#' Access variable features of liger object
#' @param x \linkS4class{liger} object
#' @rdname var.features
#' @export
setMethod("var.features", signature(x = "liger"),
          function(x) x@var.features)

#' @rdname var.features
#' @export
setReplaceMethod(
    "var.features",
    signature(x = "liger", check = "ANY", value = "character"),
    function(x, check = TRUE, value) {
        x@var.features <- value
        for (d in names(x)) {
            ld <- dataset(x, d)
            feature.meta(ld, check = FALSE)$selected <- rownames(ld) %in% value
            datasets(x)[[d]] <- ld
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
#' @param ... not used by this method
#' @return A \code{data.frame} combining \code{cell.meta(x)} with
#' \code{H} and \code{H.norm} matrices, if any of the latters exist.
#' @name fortify.liger
#' @importFrom ggplot2 fortify
#' @examples
#' \dontrun{
#' library(ggplot2)
#' # Create liger object, ligerex
#' ggplot(ligerex,  aes(x = dataset, y = nUMI)) + geom_violin()
#' }
NULL

#' @rdname fortify.liger
#' @export
#' @method fortify liger
fortify.liger <- function(x) {
    df <- as.data.frame(cell.meta(x))
    if (!is.null(x@H)) df <- cbind(df, x@H)
    if (!is.null(x@H.norm)) df <- cbind(df, x@H.norm)
    df
}
