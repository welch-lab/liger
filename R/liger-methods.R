#' Check if given liger object if under new implementation
#' @param object A liger object
#' @return \code{TRUE} if the version of \code{object} is later than or equal to
#' 1.99.0. Otherwise \code{FALSE}. It raises an error if input object is not of
#' \linkS4class{liger} class.
#' @export
#' @examples
#' is.newLiger(pbmc) # TRUE
is.newLiger <- function(object) {
    if (!inherits(object, "liger"))
        cli::cli_abort("{.var object} is not even of {.cls liger} class.")
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
#' # rliger generics
#' ## Retrieving dataset(s), replacement methods available
#' datasets(pbmcPlot)
#' dataset(pbmcPlot, "ctrl")
#' dataset(pbmcPlot, 2)
#'
#' ## Retrieving cell metadata, replacement methods available
#' cellMeta(pbmcPlot)
#' head(pbmcPlot[["nUMI"]])
#'
#' ## Retrieving dimemtion reduction matrix
#' head(dimRed(pbmcPlot, "UMAP"))
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
#' pbmcPlot2 <- pbmcPlot
#' names(pbmcPlot2) <- paste0(names(pbmcPlot), 2)
#' c(pbmcPlot, pbmcPlot2)
#'
#' library(ggplot2)
#' ggplot(pbmcPlot, aes(x = UMAP_1, y = UMAP_2)) + geom_point()
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
        cat(paste0("dimReds(", length(object@dimReds), "): "))
        cat(.collapseLongNames(names(object@dimReds)), "\n")
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

#' Subset liger with brackets
#' @name sub-liger
#' @param x A \linkS4class{liger} object
#' @param i Feature subscriptor, passed to \code{featureIdx} of
#' \code{\link{subsetLiger}}.
#' @param j Cell subscriptor, passed to \code{cellIdx} of
#' \code{\link{subsetLiger}}.
#' @param ... Additional arguments passed to \code{\link{subsetLiger}}.
#' @export
#' @return Subset of \code{x} with specified features and cells.
#' @seealso \code{\link{subsetLiger}}
#' @method [ liger
#' @examples
#' pbmcPlot[varFeatures(pbmcPlot)[1:10], 1:10]
`[.liger` <- function(x, i, j, ...) {
    if (missing(i)) i <- NULL
    if (missing(j)) j <- NULL
    subsetLiger(x, featureIdx = i, cellIdx = j, ...)
}
# setMethod(
#     "[",
#     signature(x = "liger", i = "character", j = "missing"),
#     function(x, i, j, ...) subsetLiger(x, featureIdx = i, cellIdx = NULL, ...)
# )
#
# #' @export
# #' @rdname liger-class
# setMethod(
#     "[",
#     signature(x = "liger", i = "missing", j = "index"),
#     function(x, i, j, ...) subsetLiger(x, featureIdx = NULL, cellIdx = j, ...)
# )
#
# #' @export
# #' @rdname liger-class
# setMethod(
#     "[",
#     signature(x = "liger", i = "character", j = "index"),
#     function(x, i, j, ...) subsetLiger(x, featureIdx = i, cellIdx = j, ...)
# )

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Datasets ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
setMethod("dataset", signature(x = "liger", dataset = "character_OR_NULL"),
          function(x, dataset = NULL) {
              if (is.null(dataset)) return(datasets(x)[[1]])
              else {
                  if (!dataset %in% names(x)) {
                      cli::cli_abort("Specified dataset name {.val {dataset}} not found in {.cls liger} object")
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
                 function(x, dataset, type = NULL, qc = TRUE, value) {
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
                         cli::cli_alert_info("Filling in NAs to {.field H.norm} matrix")
                         H.normNew <- matrix(
                             NA, ncol(value), ncol(x@H.norm),
                             dimnames = list(colnames(value), NULL))
                         x@H.norm <- rbind(x@H.norm, H.normNew)
                     }
                     if (length(x@dimReds) != 0) {
                         cli::cli_alert_info("Filling in NAs to {.field dimReds}")
                         for (dr in names(x@dimReds)) {
                             x@dimReds[[dr]] <- rbind(
                                 x@dimReds[[dr]],
                                 matrix(NA, ncol(value), ncol(x@dimReds[[dr]]),
                                        dimnames = list(colnames(value), NULL))
                             )
                         }
                     }
                     methods::validObject(x)
                     if (qc) x <- runGeneralQC(x, useDatasets = dataset,
                                               verbose = FALSE)
                     x
                 })

#' @export
#' @rdname liger-class
setReplaceMethod("dataset", signature(x = "liger", dataset = "character",
                                      type = "ANY", qc = "ANY",
                                      value = "matrixLike"),
                 function(x, dataset,
                          type = c("rawData", "normData"),
                          qc = FALSE,
                          value) {
                     type <- match.arg(type)
                     if (!all(startsWith(colnames(value), paste0(dataset, "_")))) {
                         cli::cli_alert_warning(
                             "Colnames of {.var value} do not all start with {.val {dataset}_}.
                                  Prefix added."
                         )
                         colnames(value) <- paste0(dataset, "_", colnames(value))
                     }
                     if (type == "rawData") {
                         ld <- createLigerDataset(rawData = value)
                     } else {
                         ld <- createLigerDataset(normData = value)
                     }
                     dataset(x, dataset, qc = qc) <- ld
                     x
                 })

#' @export
#' @rdname liger-class
setReplaceMethod(
    "dataset",
    signature(x = "liger", dataset = "character", type = "missing", qc = "ANY",
              value = "NULL"),
    function(x, dataset, type = NULL, qc = TRUE, value) {
        if (!dataset %in% names(x)) {
            cli::cli_alert_warning(
                "Specified dataset name {.val {dataset}} not found in {.cls liger} object. Nothing would happen.")
        } else {
            idxToRemove <- x$dataset == dataset
            x@datasets[[dataset]] <- NULL
            x@cellMeta <- x@cellMeta[!idxToRemove, , drop = FALSE]
            x@H.norm <- x@H.norm[!idxToRemove, , drop = FALSE]
            x@cellMeta$dataset <- droplevels(x@cellMeta$dataset)
            for (i in seq_along(x@dimReds)) {
                x@dimReds[[i]] <- x@dimReds[[i]][!idxToRemove, , drop = FALSE]
            }
        }
        x
    }
)

#' @rdname liger-class
#' @export
#' @method names liger
names.liger <- function(x) {
    names(x@datasets)
}

#' @rdname liger-class
#' @export
#' @method names<- liger
`names<-.liger` <- function(x, value) {
    originalNames <- names(x)
    if (identical(value, originalNames)) return(x)
    if (any(duplicated(value)))
        cli::cli_abort("Duplicated dataset names are not allowed.")

    dataset.idx <- lapply(originalNames, function(n) {
        x$dataset == n
    })
    x@cellMeta$dataset <- as.character(x@cellMeta$dataset)
    for (i in seq_along(value)) {
        x@cellMeta$dataset[dataset.idx[[i]]] <- value[i]
    }
    x@cellMeta$dataset <- factor(x@cellMeta$dataset, levels = value)
    names(x@datasets) <- value
    x
}

#' @rdname liger-class
#' @export
#' @method length liger
length.liger <- function(x) {
    .checkObjVersion(x)
    length(x@datasets)
}

#' @rdname liger-class
#' @export
#' @param use.names Whether returned vector should be named with dataset names.
#' @method lengths liger
lengths.liger <- function(x, use.names = TRUE) {
    .checkObjVersion(x)
    len <- sapply(x@datasets, ncol)
    if (isTRUE(use.names)) names(len) <- names(x@datasets)
    else names(len) <- NULL
    return(len)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Cell metadata ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.subsetCellMeta <- function(
        object,
        columns = NULL,
        cellIdx = NULL,
        as.data.frame = FALSE,
        ...
) {
    full <- object@cellMeta
    if (isTRUE(as.data.frame)) res <- .DataFrame.as.data.frame(full)
    else res <- full
    if (!is.null(columns)) {
        notFound <- !columns %in% colnames(res)
        if (any(notFound)) {
            cli::cli_alert_danger(
                "Specified variables from cellMeta not found: {.val {columns[notFound]}}")
            columns <- columns[!notFound]
        }
        res <- res[, columns, ...]
    }
    if (length(columns) == 1) {
        if (!is.null(dim(res))) {
            rownames(res) <- rownames(full)
        } else {
            names(res) <- rownames(full)
        }
    }
    if (!is.null(cellIdx)) {
        cellIdx <- .idxCheck(object, idx = cellIdx, orient = "cell")
        if (is.vector(res) || is.factor(res)) res <- res[cellIdx]
        else if (!is.null(dim(res))) res <- res[cellIdx, , ...]
        else cli::cli_abort("Result before idx subscription corrupted")
    }
    return(res)
}

#' @export
#' @rdname liger-class
setMethod(
    "cellMeta",
    signature(x = "liger", columns = "NULL"),
    function(x, columns = NULL, useDatasets = NULL, cellIdx = NULL, as.data.frame = FALSE, ...) {
        NULL
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "cellMeta",
    signature(x = "liger", columns = "character"),
    function(x, columns = NULL, useDatasets = NULL, cellIdx = NULL, as.data.frame = FALSE, ...) {
        if (is.null(cellIdx) && !is.null(useDatasets)) {
            if (!is.character(useDatasets)) useDatasets <- names(x)[useDatasets]
            cellIdx <- x@cellMeta$dataset %in% useDatasets
        }
        .subsetCellMeta(x, columns = columns, cellIdx = cellIdx,
                        as.data.frame = as.data.frame, ...)
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "cellMeta",
    signature(x = "liger", columns = "missing"),
    function(x, columns = NULL, useDatasets = NULL, cellIdx = NULL, as.data.frame = FALSE, ...) {
        if (is.null(cellIdx) && !is.null(useDatasets)) {
            if (!is.character(useDatasets)) useDatasets <- names(x)[useDatasets]
            cellIdx <- x@cellMeta$dataset %in% useDatasets
        }
        .subsetCellMeta(x, columns = NULL, cellIdx = cellIdx,
                        as.data.frame = as.data.frame, ...)
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "cellMeta",
    signature(x = "liger", columns = "missing"),
    function(x, columns = NULL, useDatasets = NULL, cellIdx = NULL, check = FALSE, value) {
        if (!inherits(value, "DFrame"))
            value <- S4Vectors::DataFrame(value)
        x@cellMeta <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-class
#' @param inplace For \code{cellMeta<-} method, when \code{columns} is for
#' existing variable and \code{useDatasets} or \code{cellIdx} indicate partial
#' insertion to the object, whether to by default (\code{TRUE}) in-place insert
#' \code{value} into the variable for selected cells or to replace the whole
#' variable with non-selected part left as NA.
setReplaceMethod(
    "cellMeta",
    signature(x = "liger", columns = "character"),
    function(x, columns = NULL, useDatasets = NULL, cellIdx = NULL,
             inplace = TRUE, check = FALSE, value) {
        # 1 - check cell selection
        if (is.null(cellIdx) && !is.null(useDatasets)) {
            if (!is.character(useDatasets)) useDatasets <- names(x)[useDatasets]
            cellIdx <- which(x@cellMeta$dataset %in% useDatasets)
        } else {
            cellIdx <- .idxCheck(x, cellIdx, "cell")
        }
        if (length(cellIdx) == 0)
            cli::cli_abort("No cell selected with either {.val cellIdx} or {.var useDatasets}.")

        # 2 - check value matching or length/dimension
        barcodes <- colnames(x)[cellIdx]
        if (is.null(dim(value))) {
            # Vector/factor like
            if (is.null(names(value))) {
                # No name matching, require exact length
                value <- .checkArgLen(value, n = length(cellIdx), class = c("vector", "factor"))
            } else {
                if (!all(endsWith(barcodes, names(value)))) {
                    cli::cli_abort(
                        c("x" = "{.code names(value)} do match to cells selected. ",
                          "i" = "The first three given names: {.val {names(value)[1:3]}}",
                          "i" = "The first three selected names: {.val {barocdes[1:3]}}")
                    )
                }
                names(value) <- barcodes
            }
        } else {
            # matrix like
            if (is.null(rownames(value))) {
                # no rowname matching, require extact nrow
                if (nrow(value) != length(cellIdx)) {
                    cli::cli_abort(
                        "{.code nrow(value)} ({nrow(value)}) does not match with cells selected ({length(cellIdx)}).")
                }
            } else {
                if (!all(barcodes %in% rownames(value))) {
                    cli::cli_abort(
                        c("{.code rownames(value)} do not contain all cells selected. ",
                          "These are not involved: ",
                          "{.val {barcodes[!barcodes %in% rownames(value)]}}")
                    )
                }
                value <- value[barcodes, , drop = FALSE]
            }
        }

        # 3 - Insert value
        if (length(cellIdx) == ncol(x)) {
            x@cellMeta[[columns]] <- value
        } else if (length(cellIdx) < ncol(x)) {
            # Partial replacement/addition
            if (!columns %in% names(x@cellMeta)) {
                # Add new variable
                if (is.null(dim(value))) {
                    # vector/factor
                    x@cellMeta[[columns]] <- NA
                    if (is.factor(value)) {
                        charValue <- as.character(value)
                        x@cellMeta[[columns]][cellIdx] <- charValue
                        x@cellMeta[[columns]] <- factor(x@cellMeta[[columns]],
                                                        levels = levels(value))
                    } else {
                        x@cellMeta[[columns]][cellIdx] <- value
                    }
                } else {
                    # matrix like
                    x@cellMeta[[columns]] <- matrix(
                        NA, ncol(x), ncol(value),
                        dimnames = list(colnames(x), colnames(value))
                    )
                    x@cellMeta[[columns]][cellIdx,] <- value
                }
            } else {
                if (isTRUE(inplace)) {
                    # Modifying existing variable
                    if (is.null(dim(value)) && is.null(dim(x@cellMeta[[columns]]))) {
                        # Both are 1-D
                        if (is.factor(value) && is.factor(x@cellMeta[[columns]])) {
                            charVar <- as.character(x@cellMeta[[columns]])
                            charVar[cellIdx] <- as.character(value)
                            x@cellMeta[[columns]] <-
                                factor(
                                    charVar,
                                    levels = unique(c(levels(x@cellMeta[[columns]]),
                                                      levels(value)))
                                )
                        } else {
                            x@cellMeta[[columns]][cellIdx] <- value
                        }
                    } else if (!is.null(dim(value)) && !is.null(dim(x@cellMeta[[columns]]))) {
                        # Both are dimensional
                        if (ncol(value) != ncol(x@cellMeta[[columns]])) {
                            cli::cli_abort("Cannot insert value to a variable of different dimensionality")
                        }
                        x@cellMeta[[columns]][cellIdx,] <- value
                    } else {
                        cli::cli_abort("Cannot insert value to a variable of different dimensionality")
                    }
                } else {
                    x@cellMeta[[columns]] <- NULL
                    # Remove and go to "Add new variable" case above
                    cellMeta(x, columns = columns, cellIdx = cellIdx, check = check) <- value
                }
            }
        } else {
            cli::cli_abort("{.var cellIdx} pointing to more cells than available")
        }
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)


#' @export
#' @rdname liger-class
setMethod("rawData", c("liger", "ANY"),
          function(x, dataset = NULL) {
              if (is.null(dataset)) {
                  getMatrix(x, slot = "rawData", returnList = TRUE)
              } else {
                  getMatrix(x, "rawData", dataset = dataset)
              }
          }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "rawData",
    signature(x = "liger", dataset = "ANY", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, dataset = NULL, check = TRUE, value) {
        dataset <- .checkUseDatasets(x, dataset)
        if (length(dataset) != 1) cli::cli_abort("Need to specify one dataset to insert.")
        if (isH5Liger(x, dataset))
            cli::cli_abort("Cannot replace slot with in-memory data for H5 based object.")
        x@datasets[[dataset]]@rawData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "rawData",
    signature(x = "liger", dataset = "ANY", check = "ANY", value = "H5D"),
    function(x, dataset = NULL, check = TRUE, value) {
        dataset <- .checkUseDatasets(x, dataset)
        if (length(dataset) != 1) cli::cli_abort("Need to specify one dataset to insert.")
        if (!isH5Liger(x, dataset))
            cli::cli_abort("Cannot replace slot with on-disk data for in-memory object.")
        x@datasets[[dataset]]@rawData <- value
        if (isTRUE(check)) methods::validObject(x@datasets[[dataset]])
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-class
setMethod("normData", c("liger", "ANY"),
          function(x, dataset = NULL) {
              if (is.null(dataset)) {
                  getMatrix(x, slot = "normData", returnList = TRUE)
              } else {
                  getMatrix(x, "normData", dataset = dataset)
              }
          }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "normData",
    signature(x = "liger", dataset = "ANY", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, dataset = NULL, check = TRUE, value) {
        dataset <- .checkUseDatasets(x, dataset)
        if (length(dataset) != 1) cli::cli_abort("Need to specify one dataset to insert.")
        if (isH5Liger(x, dataset))
            cli::cli_abort("Cannot replace slot with in-memory data for H5 based object.")
        x@datasets[[dataset]]@normData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "normData",
    signature(x = "liger", dataset = "ANY", check = "ANY", value = "H5D"),
    function(x, dataset = NULL, check = TRUE, value) {
        dataset <- .checkUseDatasets(x, dataset)
        if (length(dataset) != 1) cli::cli_abort("Need to specify one dataset to insert.")
        if (!isH5Liger(x, dataset))
            cli::cli_abort("Cannot replace slot with on-disk data for in-memory object.")
        x@datasets[[dataset]]@normData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "scaleData",
    signature(x = "liger", dataset = "ANY"),
    function(x, dataset = NULL) {
        if (is.null(dataset)) {
            getMatrix(x, slot = "scaleData", returnList = TRUE)
        } else {
            getMatrix(x, "scaleData", dataset = dataset)
        }
    }
)


#' @export
#' @rdname liger-class
setReplaceMethod(
    "scaleData",
    signature(x = "liger", dataset = "ANY", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, dataset = NULL, check = TRUE, value) {
        dataset <- .checkUseDatasets(x, dataset)
        if (length(dataset) != 1) cli::cli_abort("Need to specify one dataset to insert.")
        if (isH5Liger(x, dataset))
            cli::cli_abort("Cannot replace slot with in-memory data for H5 based object.")
        x@datasets[[dataset]]@scaleData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "scaleData",
    signature(x = "liger", dataset = "ANY", check = "ANY", value = "H5D"),
    function(x, dataset = NULL, check = TRUE, value) {
        dataset <- .checkUseDatasets(x, dataset)
        if (length(dataset) != 1) cli::cli_abort("Need to specify one dataset to insert.")
        if (!isH5Liger(x, dataset))
            cli::cli_abort("Cannot replace slot with on-disk data for in-memory object.")
        x@datasets[[dataset]]@scaleData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "scaleData",
    signature(x = "liger", dataset = "ANY", check = "ANY", value = "H5Group"),
    function(x, dataset = NULL, check = TRUE, value) {
        dataset <- .checkUseDatasets(x, dataset)
        if (length(dataset) != 1) cli::cli_abort("Need to specify one dataset to insert.")
        if (!isH5Liger(x, dataset))
            cli::cli_abort("Cannot replace slot with on-disk data for in-memory object.")
        x@datasets[[dataset]]@scaleData <- value
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "scaleUnsharedData",
    signature(x = "liger", dataset = "character"),
    function(x, dataset) {
        scaleUnsharedData(dataset(x, dataset))
    }
)


#' @export
#' @rdname liger-class
setMethod(
    "scaleUnsharedData",
    signature(x = "liger", dataset = "numeric"),
    function(x, dataset) {
        scaleUnsharedData(dataset(x, dataset))
    }
)


#' @export
#' @rdname liger-class
setReplaceMethod(
    "scaleUnsharedData",
    signature(x = "liger", dataset = "ANY", check = "ANY", value = "matrixLike_OR_NULL"),
    function(x, dataset = NULL, check = TRUE, value) {
        dataset <- .checkUseDatasets(x, dataset)
        if (length(dataset) != 1) cli::cli_abort("Need to specify one dataset to insert.")
        ld <- dataset(x, dataset)
        scaleUnsharedData(ld) <- value
        x@datasets[[dataset]] <- ld
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "scaleUnsharedData",
    signature(x = "liger", dataset = "ANY", check = "ANY", value = "H5D"),
    function(x, dataset = NULL, check = TRUE, value) {
        dataset <- .checkUseDatasets(x, dataset)
        if (length(dataset) != 1) cli::cli_abort("Need to specify one dataset to insert.")
        ld <- dataset(x, dataset)
        scaleUnsharedData(ld) <- value
        x@datasets[[dataset]] <- ld
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "scaleUnsharedData",
    signature(x = "liger", dataset = "ANY", check = "ANY", value = "H5Group"),
    function(x, dataset = NULL, check = TRUE, value) {
        dataset <- .checkUseDatasets(x, dataset)
        if (length(dataset) != 1) cli::cli_abort("Need to specify one dataset to insert.")
        if (!isH5Liger(x, dataset))
            cli::cli_abort("Cannot replace slot with on-disk data for in-memory object.")
        ld <- dataset(x, dataset)
        scaleUnsharedData(ld) <- value
        x@datasets[[dataset]] <- ld
        if (isTRUE(check)) methods::validObject(x)
        x
    }
)


#' @export
#' @rdname liger-class
setMethod(
    "getMatrix", signature(x = "liger"),
    function(x,
             slot = c("rawData", "normData", "scaleData",
                      "scaleUnsharedData", "H", "V", "U", "A", "B",
                      "W", "H.norm"),
             dataset = NULL,
             returnList = FALSE) {
        slot <- match.arg(slot)
        if (slot == "W") return(x@W)
        if (slot == "H.norm") return(x@H.norm)
        if (is.null(dataset)) {
            return(lapply(datasets(x), function(ld) getMatrix(ld, slot)))
        } else {
            if (length(dataset) == 1) {
                if (isTRUE(returnList)) {
                    result <- list(getMatrix(dataset(x, dataset), slot))
                    names(result) <- dataset
                    return(result)
                } else if (isFALSE(returnList))
                    return(getMatrix(dataset(x, dataset), slot))
            } else {
                lds <- datasets(x)[dataset]
                return(lapply(lds, function(ld) getMatrix(ld, slot)))
            }
        }
    })

#' @export
#' @rdname liger-class
setMethod("getH5File",
          signature = signature(x = "liger", dataset = "ANY"),
          function(x, dataset = NULL) {
              if (is.null(dataset)) dataset <- names(x)
              dataset <- .checkUseDatasets(x, dataset)
              results <- lapply(datasets(x)[dataset],
                                function(ld) h5fileInfo(ld, "H5File"))
              if (length(results) == 1) results <- results[[1]]
              results
          })


#' Get cell metadata variable
#' @name sub-sub-liger
#' @param x A \linkS4class{liger} object
#' @param i Name or numeric index of cell meta data to fetch
#' @param ... Anything that \code{S4Vectors::\link[S4Vectors]{DataFrame}}
#' method allows.
#' @export
#' @method [[ liger
#' @return If \code{i} is given, the selected metadata will be returned; if it
#' is missing, the whole cell metadata table in
#' \code{S4Vectors::\link[S4Vectors]{DataFrame}} class will be returned.
#' @examples
#' # Retrieve whole cellMeta
#' pbmc[[]]
#' # Retrieve a variable
#' pbmc[["dataset"]]
`[[.liger` <- function(x, i, ...) {
    if (missing(i)) return(cellMeta(x))
    else return(x@cellMeta[[i, ...]])
}

#' Set cell metadata variable
#' @rdname liger-class
#' @param i Name or numeric index of cell meta variable to be replaced
#' @param value Metadata value to be inserted
#' @return Input liger object updated with replaced/new variable in
#' \code{cellMeta(x)}.
#' @export
#' @method [[<- liger
#' @examples
#' cellMeta(pbmc)
#' # Add new variable
#' pbmc[["newVar"]] <- 1
#' cellMeta(pbmc)
#' # Change existing variable
#' pbmc[["newVar"]][1:3] <- 1:3
#' cellMeta(pbmc)
`[[<-.liger` <- function(x, i, value) {
    name <- if (is.character(i)) i else colnames(x@cellMeta)[i]
    if (name == "dataset") {
        cli::cli_abort(
            c("x" = "Cannot directly modify {.var dataset} variable in {.cls liger} object.",
              "i" = "Please use {.code names(x) <- value} instead.")
        )
    }
    x@cellMeta[[i]] <- value
    methods::validObject(x)
    return(x)
}

#' @export
#' @method .DollarNames liger
.DollarNames.liger <- function(x, pattern = "") {
    grep(pattern, colnames(x@cellMeta), value = TRUE)
}


#' @export
#' @rdname liger-class
#' @method $ liger
`$.liger` <- function(x, name) {
    if (!name %in% colnames(cellMeta(x))) NULL
    else cellMeta(x, columns = name)
}

#' @export
#' @rdname liger-class
#' @method $<- liger
`$<-.liger` <- function(x, name, value) {
    cellMeta(x, columns = name) <- value
    return(x)
}

#' @export
#' @rdname liger-class
#' @param droplevels Whether to remove unused cluster levels from the factor
#' object fetched by \code{defaultCluster()}. Default \code{FALSE}.
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
        if (length(value) == 1) {
            # If doing defaultCluster(obj) <- "someName"
            if (!is.null(name)) {
                cli::cli_alert_danger("Cannot have {.code name} when selecting a variable with {.code value}.")
            }
            name <- value
            if (!name %in% colnames(cellMeta(x))) {
                cli::cli_abort("Selected variable does not exist in {.code cellMeta(x)}.")
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
            cli::cli_abort("Length of {.code value} does not match with the number of cells.")
        }
        if (is.null(name)) {
            cli::cli_alert_info(
                c("Storing given cluster labels to {.code cellMeta(x)} field: ",
                  "{.val defaultCluster}.")
            )
            name <- "defaultCluster"
        }
        if (is.null(names(value))) names(value) <- colnames(x)[cellIdx]
        else {
            if (all(names(value) %in% colnames(x)[cellIdx])) {
                value <- value[colnames(x)[cellIdx]]
            } else {
                cli::cli_abort(
                    c(x = "Not all {.code names(value)} match with target cells: ",
                      "{.val {names(value)[!names(value) %in% colnames(x)[cellIdx]]}}",
                      i = "Please have a check or try {.code unname(value)}.")
                )
            }
        }
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
setMethod(
    "dimReds",
    signature = c(x = "liger"),
    function(x) x@dimReds
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "dimReds",
    signature(x = "liger", value = "list"),
    function(x, value) {
        x@dimReds <- value
        methods::validObject(x)
        return(x)
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "dimRed",
    signature = c(x = "liger", name = "missing_OR_NULL"),
    function(x, name = NULL, useDatasets = NULL, cellIdx = NULL, ...) {
        name <- x@uns$defaultDimRed
        dimred <- NULL
        if (is.null(name)) {
            if (length(x@dimReds) > 0) {
                cli::cli_alert_warning(
                    "No default {.field dimRed} recorded. Returning the first available."
                )
                dimred <- dimRed(x, name = 1, useDatasets = useDatasets,
                                 cellIdx = cellIdx, ...)
            } else {
                cli::cli_abort(
                    "No {.field dimRed} available in this {.cls liger} object."
                )
            }
        } else {
            dimred <- dimRed(x, name = name, useDatasets = useDatasets,
                             cellIdx = cellIdx, ...)
        }
        return(dimred)
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "dimRed",
    signature = c(x = "liger", name = "index"),
    function(x, name, useDatasets = NULL, cellIdx = NULL, ...) {
        if (is.null(useDatasets) && is.null(cellIdx)) {
            cellIdx <- seq_len(ncol(x))
        } else if (!is.null(cellIdx)) {
            cellIdx <- .idxCheck(x, cellIdx, "cell")
        } else if (!is.null(useDatasets)) {
            useDatasets <- .checkUseDatasets(x, useDatasets)
            cellIdx <- which(x$dataset %in% useDatasets)
        }

        name <- .findDimRedName(x, name, stopOnNull = TRUE)
        dimred <- x@dimReds[[name]]
        dimred <- dimred[cellIdx, , drop = FALSE]
        rownames(dimred) <- colnames(x)[cellIdx]
        colnames(dimred) <- paste0(name, "_", seq_len(ncol(dimred)))
        return(dimred)
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "dimRed",
    signature(x = "liger", name = "index", value = "NULL"),
    function(x, name = NULL, useDatasets = NULL, cellIdx = NULL, ..., value = NULL) {
        name <- .findDimRedName(x, name, stopOnNull = TRUE, returnFirst = FALSE)
        x@dimReds[[name]] <- NULL
        if (name %in% x@uns$defaultDimRed) x@uns$defaultDimRed <- NULL
        return(x)
    }
)

#' @export
#' @rdname liger-class
#' @param asDefault Whether to set the inserted dimension reduction matrix as
#' default for visualization methods. Default \code{NULL} sets it when no
#' default has been set yet, otherwise does not change current default.
setReplaceMethod(
    "dimRed",
    signature(x = "liger", name = "character", value = "matrixLike"),
    function(x, name = NULL, useDatasets = NULL, cellIdx = NULL, asDefault = NULL, inplace = FALSE, ..., value) {
        if (is.null(useDatasets) && is.null(cellIdx)) {
            cellIdx <- seq_len(ncol(x))
        } else if (!is.null(cellIdx)) {
            cellIdx <- .idxCheck(x, cellIdx, "cell")
        } else if (!is.null(useDatasets)) {
            useDatasets <- .checkUseDatasets(x, useDatasets)
            cellIdx <- which(x$dataset %in% useDatasets)
        }

        if (!name %in% names(x@dimReds) || isFALSE(inplace)) {
            # Totally new thing or just replace
            init <- matrix(
                data = NA,
                nrow = ncol(x), ncol = ncol(value),
                dimnames = list(
                    colnames(x),
                    paste0(name, "_", seq_len(ncol(value)))
                )
            )
        } else {
            # Partial insertion
            init <- dimRed(x, name = name)
            if (ncol(init) != ncol(value)) {
                cli::cli_abort(
                    "Cannot partially insert {ncol(value)} columns to {ncol(init)} columns inplace, at {.field dimReds}: {.val {name}}"
                )
            }
        }

        value <- as.matrix(value)
        if (nrow(value) != length(cellIdx)) {
            cli::cli_abort(
                "{.code nrow(value)} does not match with the number of cells selected."
            )
        }
        if (is.null(rownames(value))) {
            cli::cli_alert_warning(
                "No rownames detected. Assume cells match to the same order as in the object."
            )
        } else {
            if (!all(endsWith(colnames(x)[cellIdx], rownames(value)))) {
                cli::cli_abort(
                    "Cell identifiers in {.var value} do not match to those in the object"
                )
            }
        }
        rownames(value) <- colnames(x)[cellIdx]
        colnames(value) <- paste0(name, "_", seq_len(ncol(value)))
        init[rownames(value), ] <- value
        x@dimReds[[name]] <- init
        if (is.null(asDefault)) {
            if (!is.null(x@uns$defaultDimRed)) asDefault <- FALSE
            else asDefault <- TRUE
        }
        if (isTRUE(asDefault)) defaultDimRed(x) <- name
        return(x)
    }
)

#' @export
#' @rdname liger-class
setMethod(
    "defaultDimRed",
    signature(x = "liger", useDatasets = "ANY"),
    function(x, useDatasets = NULL, cellIdx = cellIdx) {
        name <- x@uns$defaultDimRed
        if (is.null(name)) return(NULL)
        else dimRed(x, name = name, useDatasets = useDatasets, cellIdx = cellIdx)
    }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "defaultDimRed",
    signature(x = "liger", value = "character"),
    function(x, value) {
        if (length(value) != 1) {
            cli::cli_abort("Can only set one {.field dimRed} as default.")
        }
        value <- .findDimRedName(x, value, stopOnNull = TRUE)
        x@uns$defaultDimRed <- value
        return(x)
    }
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
                cli::cli_alert_warning(
                    "Not all variable features passed are found in datasets: {.val {names(x)[!checkResult]}}"
                )
            }
        }
        x
    }
)


#' @export
#' @rdname liger-class
setMethod("varUnsharedFeatures", signature(x = "liger"),
          function(x, dataset = NULL) {
              dataset <- .checkUseDatasets(x, dataset)
              vufList <- lapply(dataset, function(d) x@datasets[[d]]@varUnsharedFeatures)

              if (length(vufList) == 1) return(vufList[[1]])
              else {
                  names(vufList) <- dataset
                  return(vufList)
              }
          }
)

#' @export
#' @rdname liger-class
setReplaceMethod(
    "varUnsharedFeatures",
    signature(x = "liger", dataset = "ANY", check = "ANY", value = "character"),
    function(x, dataset, check = TRUE, value) {
        dataset <- .checkUseDatasets(x, dataset)
        x@datasets[[dataset]]@varUnsharedFeatures <- value
        if (isTRUE(check)) {
            if (!all(value %in% rownames(x@datasets[[dataset]]))) {
                cli::cli_alert_warning("Not all features passed are found in dataset {.val {dataset}}.")
            }
        }
        return(x)
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
    drs <- Reduce(cbind, model@dimReds)
    if (!is.null(drs)) df <- cbind(df, drs)
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
        cli::cli_abort("Can only combine {.cls liger} objects with {.fn c} method for now.")
    allNames <- unlist(lapply(objList, names))
    if (any(duplicated(allNames)))
        cli::cli_abort(
            c("x" = "Cannot combine {.cls liger} objects with duplicated dataset names.",
              "i" = "Dataset names of an individual {.cls liger} object can be modified with {.code names(x) <- value}.")
        )

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


#' @export
#' @rdname liger-class
setMethod(
    "commands",
    signature(x = "liger", funcName = "ANY", arg = "ANY"),
    function(x, funcName = NULL, arg = NULL) {
        if (is.null(funcName)) return(names(x@commands))
        cmdIdx <- c()
        for (n in funcName) {
            pat <- paste0("^", n)
            cmdIdx <- c(cmdIdx, grep(pat, names(x@commands)))
        }
        cmdIdx <- sort(unique(cmdIdx))
        result <- x@commands[cmdIdx]

        if (length(result) == 1) result <- result[[1]]

        if (!is.null(arg)) {
            if (is.list(result))
                result <- lapply(result, function(cmd) cmd@parameters[arg])
            else result <- unlist(result@parameters[arg])
        }
        return(result)
    }
)


#' @rdname peak
#' @export
setMethod(
    "rawPeak",
    signature(x = "liger", dataset = "character"),
    function(x, dataset) {
        atac <- dataset(x, dataset)
        rawPeak(atac)
    }
)

#' @rdname peak
#' @export
setReplaceMethod(
    "rawPeak",
    signature(x = "liger", dataset = "character"),
    function(x, dataset, check = TRUE, value) {
        ld <- dataset(x, dataset)
        rawPeak(ld, check = check) <- value
        x@datasets[[dataset]] <- ld
        x
    }
)

#' @rdname peak
#' @export
setMethod(
    "normPeak",
    signature(x = "liger", dataset = "character"),
    function(x, dataset) {
        atac <- dataset(x, dataset)
        normPeak(atac)
    }
)

#' @rdname peak
#' @export
setReplaceMethod(
    "normPeak",
    signature(x = "liger", dataset = "character"),
    function(x, dataset, check = TRUE, value) {
        ld <- dataset(x, dataset)
        normPeak(ld, check = check) <- value
        x@datasets[[dataset]] <- ld
        x
    }
)


#' @rdname coordinate
#' @export
setMethod(
    "coordinate",
    signature(x = "liger", dataset = "character"),
    function(x, dataset) {
        spatial <- dataset(x, dataset)
        coordinate(spatial)
    }
)

#' @rdname coordinate
#' @export
setReplaceMethod(
    "coordinate",
    signature(x = "liger", dataset = "character"),
    function(x, dataset, check = TRUE, value) {
        ld <- dataset(x, dataset)
        coordinate(ld, check = check) <- value
        x@datasets[[dataset]] <- ld
        x
    }
)
