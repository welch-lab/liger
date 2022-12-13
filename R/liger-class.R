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
        H = "matrix_OR_NULL",
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
#' @param datasets Named list of datasets. List elements can be one of in memory
#' expression matrix, file name (.h5), constructed ligerDataset object.
#' @param cell.meta data.frame
#' @param var.features vector of feature names
#' @param H matrix
#' @param H.norm matrix
#' @export
createLiger <- function(
        raw.data = NULL,
        format.type = NULL,
        modal = NULL,
        cell.meta = NULL,
        take.gene.union = FALSE,
        remove.missing = FALSE,
        h5.format = "10X",
        data.name = NULL,
        indices.name = NULL,
        indptr.name = NULL,
        genes.name = NULL,
        barcodes.name = NULL,
        verbose = TRUE
        )
{
    if (!is.list(raw.data)) stop("`raw.data` has to be a named list.")
    nData <- length(raw.data)
    if (length(format.type) == 1) format.type <- rep(format.type, nData)
    else if (length(format.type) != nData)
        stop("Wrong length of `format.type`. ",
             "Specify only 1 or match the length of `datasets`. ",
             "See ?createLiger for valid options.")
    if (length(modal) == 1) modal <- rep(modal, nData)
    else if (length(modal) != nData)
        stop("Wrong length of `modal`. ",
             "Specify only 1 or match the length of `datasets`. ",
             "See ?createLiger for valid options.")
    datasets <- mapply(function(data, type, modal) {
        createLigerDataset(raw.data = data,
                           format = format.type,
                           modal = modal)
    }, raw.data, format.type, modal)

    if (isTRUE(take.gene.union)) {
        merged.data <- MergeSparseDataAll(raw.data)
        if (isTRUE(remove.missing)) {
            missing_genes <- rowSums(merged.data) == 0
            if (length(missing_genes) > 0) {
                if (isTRUE(verbose)) {
                    message(
                        "Removing ",
                        sum(missing_genes),
                        " genes not expressed in any cells across all datasets."
                    )
                }
                if (length(missing_genes) < 25) {
                    if (isTRUE(verbose)) {
                        message(rownames(merged.data)[missing_genes])
                    }
                }
                merged.data <- merged.data[!missing_genes,]
            }
        }
        raw.data <- lapply(raw.data, function(x) {
            merged.data[, colnames(x)]
        })
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
    #if (is.null(uns)) uns <- list()
    obj <- methods::new("liger",
                        datasets = datasets,
                        cell.meta = cell.meta)
    if (isTRUE(remove.missing)) {
        obj <- removeMissingObs(obj, use.cols = TRUE, verbose = verbose)
        # remove missing genes if not already merged
        if (!isTRUE(take.gene.union)) {
            obj <- removeMissingObs(obj, use.cols = FALSE, verbose = verbose)
        }
    }
    obj <- runGeneralQC(obj)
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
        for (n in names(datasets)) {
            colnames(datasets[[n]]) <- allBarcodes[datasetVar == n]
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

setGeneric("datasets", function(x, ...) standardGeneric("datasets"))
setGeneric("datasets<-", function(x, ...) standardGeneric("datasets<-"))
setGeneric("dataset", function(x, name = NULL) standardGeneric("dataset"))
setGeneric("dataset<-", function(x, name, type = NULL, qc = TRUE, value) {
    standardGeneric("dataset<-")
})
setGeneric("cell.meta", function(x, ...) standardGeneric("cell.meta"))
setGeneric("cell.meta<-", function(x, ...) standardGeneric("cell.meta<-"))

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

.collapseLongNames <- function(x) {
    if (length(x) < 5) {
        return(paste(x, collapse = ", "))
    } else {
        head <- paste(x[seq(3)], collapse = ", ")
        tail <- x[length(x)]
        return(paste0(head, ", ..., ", tail))
    }
}

setMethod("names", signature(x = "liger"), function(x) {
    names(datasets(x))
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
    if (!is.null(x@H)) colnames(x@H) <- value[[2L]]
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
                     cell.meta(x)[[i, ...]] <- value
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

#' Subset liger object by cell
#' @param x liger object
#' @param i subscriber for cell, based on cell.meta(x)
setMethod(
    "[",
    signature(x = "liger", i = "ANY"),
    function(x, i, ...) {
        validObject(x)
        bc <- rownames(cell.meta(x))[i]
        datasetVar <- x$dataset[i]
        datasetUniq <- unique(as.character(datasetVar))
        datasets.new <- lapply(datasetUniq, function(d) {
            # Use character colname subscription
            dataset(x, d)[bc[datasetVar == d]]
        })
        names(datasets.new) <- datasetUniq
        methods::new(
            "liger",
            datasets = datasets.new,
            cell.meta = cell.meta(x)[i, , drop = FALSE],
            var.features = "character_OR_NULL",
            H = x@H[, i, drop = FALSE],
            H.norm = x@H.norm[, i, drop = FALSE],
            version = packageVersion("rliger")
        )
    }
)

#' Access datasets of liger object
#' @param x liger object
#' @export
setMethod("datasets", signature(x = "liger"),
          function(x) x@datasets)

#' Set datasets to liger object
#' @param x liger object
#' @param value ligerDataset object
#' @export
setReplaceMethod("datasets", signature(x = "liger"),
                 function(x, value) {
                     x@datasets <- value
                     validObject(x)
                     x
                 })

#' @title Access single dataset of liger object
#' @description TODO: Debating on whether to have a separate (differently named)
#' method to directly access matrices in embedded ligerDataset objects, or to
#' allow accessing matrices with \code{type}.
#' @rdname dataset
#' @param x \linkS4class{liger} object
#' @param name name of dataset
#' @param type Type of specified matrix. Only works when \code{value} is a
#' matrix like object. choose from \code{"raw.data"}, \code{"norm.data"} or
#' \code{"scale.data"}. Default \code{"raw.data"}.
#' @param value A matrix of a dataset to be added, or a constructed
#' \linkS4class{ligerDataset} object. \code{NULL} to remove a dataset by
#' \code{name}.
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
                     if (!is.null(x@H)) {
                         message("Filling in NAs to H matrix")
                         x@H[, new.idx] <- NA
                         colnames(x@H)[new.idx] <- colnames(value)
                     }
                     if (!is.null(x@H.norm)) {
                         message("Filling in NAs to H.norm matrix")
                         x@H.norm[, new.idx] <- NA
                         colnames(x@H.norm)[new.idx] <- colnames(value)
                     }
                     validObject(x)
                     if (qc) x <- runGeneralQC(x)
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
                         ldToRemove <- dataset(x, name)
                         bcToRemove <- colnames(ldToRemove)
                         bc.idx <- colnames(x) %in% bcToRemove
                         x@datasets[[name]] <- NULL
                         x@cell.meta <- x@cell.meta[!bc.idx, , drop = FALSE]
                         x@H <- x@H[!bc.idx, ]
                         x@H.norm <- x@H.norm[!bc.idx, ]
                     }
                     x
                 })

#' Access cell.meta of liger object
#' @param x liger object
#' @export
setMethod("cell.meta", signature(x = "liger"),
          function(x) x@cell.meta)

#' Set cell.meta to liger object
#' @param x liger object
#' @param value ligerDataset object
#' @export
setReplaceMethod("cell.meta", signature(x = "liger"),
                 function(x, value) {
                     x@cell.meta <- value
                     validObject(x)
                     x
                 })

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
