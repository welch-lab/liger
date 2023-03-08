.log <- function(..., level = 1) {
    line <- paste(rep("...", level), collapse = "")
    pref <- paste0(date(), " ", line, " ")
    indent <- paste(rep(" ", nchar(pref)), collapse = "")
    content <- list(...)
    msg <- paste(content, collapse = "")
    msg <- gsub("\n", paste0("\n", indent), msg)
    message(pref, msg)
}

.checkObjVersion <- function(object) {
    if (inherits(object, "liger")) {
        if (!is.newLiger(object))
            stop("Old version of liger object detected. Please update the ",
                 "object with command:\nobject <- convertOldLiger(object)")
    }
}

.collapseLongNames <- function(x) {
    if (length(x) < 5) {
        return(paste(x, collapse = ", "))
    } else {
        head <- paste(x[seq(3)], collapse = ", ")
        tail <- x[length(x)]
        return(paste0(head, ", ..., ", tail))
    }
}

# Returns character names of datasets
.checkUseDatasets <- function(object, useDatasets = NULL,
                              modal = NULL) {
    if (!inherits(object, "liger"))
        stop("A liger object is required.")
    if (is.null(useDatasets)) {
        if (is.null(modal)) useDatasets <- names(object)
        else {
            selection <- sapply(names(object), function(d)
                inherits(dataset(object, d), .modalClassDict[[modal]])
            )
            useDatasets <- names(object)[selection]
        }
    } else {
        if (is.numeric(useDatasets)) {
            if (max(useDatasets) > length(object)) {
                stop("Numeric dataset index out of bound. Only ",
                     length(object), " datasets exist.")
            }
            useDatasets <- unique(useDatasets)
            useDatasets <- names(object)[useDatasets]
        } else if (is.logical(useDatasets)) {
            if (length(useDatasets) != length(object)) {
                stop("Logical dataset subscription does not match the number ",
                     "of datasets (", length(object), ").")
            }
            useDatasets <- names(object)[useDatasets]
        } else if (is.character(useDatasets)) {
            if (any(!useDatasets %in% names(object))) {
                notFound <- useDatasets[!useDatasets %in% names(object)]
                stop("Specified dataset name(s) not found: ",
                     paste(notFound, collapse = ", "))
            }
        } else {
            stop("Please use a proper numeric/logical/character vector to ",
                 "select dataset to use.")
        }
        if (!is.null(modal)) {
            passing <- sapply(useDatasets, function(d) {
                inherits(dataset(object, d), .modalClassDict[[modal]])
            })
            if (!all(passing))
                stop("Not all specified datasets are of `",
                     .modalClassDict[[modal]], "` class: ",
                     paste(useDatasets[!passing], collapse = ", "))
        }
    }
    useDatasets
}

# Check if the selection of cellMeta variable is valid, in terms of existence
# and class and return all result
.fetchCellMetaVar <- function(
        object,
        variables,
        cellIdx = NULL,
        checkCategorical = FALSE,
        drop = TRUE,
        returnList = FALSE
) {
    df <- cellMeta(object, variables, cellIdx, as.data.frame = TRUE,
                    drop = FALSE)
    if (isTRUE(checkCategorical)) {
        passing <- sapply(variables, function(v) {
            vec <- df[[v]]
            if (!is.null(dim(vec))) return(FALSE)
            if (is.factor(vec)) return(TRUE)
            if (is.character(vec)) {
                if (length(unique(vec)) > 50)
                    warning("Categorical variable selection `", v,
                            "' has more than 100 unique values.",
                            immediate. = TRUE)
                return(TRUE)
            }
            if (is.numeric(vec)) {
                if (length(unique(vec)) > 50)
                    warning("Categorical variable selection `", v,
                            "` has more than 100 unique values.",
                            immediate. = TRUE)
                return(FALSE)
            }
        })
        if (!all(passing)) {
            notPassed <- variables[!passing]
            stop("The following selected variables are not considered as ",
                 "categorical. Please use something else or try converting ",
                 "them to factor class to force passing checks.\n",
                 paste(notPassed, collapse = ", "))
        }
    }
    if (isTRUE(returnList)) {
        return(as.list(df))
    }
    df <- df[,,drop = drop]
    return(df)
}

# Used for checking the subscription of cells or features and returns numeric
# subscription
.idxCheck <- function(object, idx, orient) {
    if (orient == "feature") {
        getNames <- rownames
        getNumber <- nrow
    } else if (orient == "cell") {
        getNames <- colnames
        getNumber <- ncol
    } else {
        stop("`orient` should be either 'feature' or 'cell'.")
    }
    paramName <- paste0("`", orient, "Idx`")
    if (is.character(idx)) {
        if (any(!idx %in% getNames(object))) {
            notFound <- paste(idx[!idx %in% getNames(object)],
                              collapse = ", ")
            stop(paramName, " not found in object: ", notFound)
        }
        name <- seq(getNumber(object))
        names(name) <- getNames(object)
        idx <- name[idx]
        names(idx) <- NULL
    } else if (is.logical(idx)) {
        if (length(idx) != getNumber(object)) {
            stop("Length of logical ", paramName, " does not match to ",
                 "number of ", orient, "s in `object`.")
        }
        idx <- which(idx)
    } else if (is.numeric(idx)) {
        if (max(idx) > getNumber(object)) {
            stop("Numeric ", paramName, " subscription out of bound.")
        }
    } else if (is.null(idx)) {
        idx <- seq(getNumber(object))
    } else {
        stop("Please use character, logical or numeric subscription ",
             "for ", paramName, ".")
    }
    return(idx)
}

.checkLDSlot <- function(object, slot) {
    if (!inherits(object, "ligerDataset"))
        stop("Please use a ligerDataset object.")
    avail <- c("rawData", "normData", "scaleData")
    if (is.null(slot)) {
        slot <- avail
    } else {
        if (any(!slot %in% avail)) {
            notFound <- slot[!slot %in% avail]
            stop("Specified slot not availalble: ",
                 paste(notFound, collapse = ", "), ". Use one or more from ",
                 '"rawData", "normData" or "scaleData"')
        }
        if ("rawData" %in% slot &&
            is.null(rawData(object))) {
            stop("`rawData` is not available for use.")
        }
        if ("normData" %in% slot &&
            is.null(normData(object))) {
            stop("`normData` is not available for use.")
        }
        if ("scaleData" %in% slot &&
            is.null(scaleData(object))) {
            stop("`scaleData` is not available for use.")
        }
    }
    slot
}

.checkInit <- function(
    m, nCells, nGenes, k, type = c("W", "H", "V", "A", "V")
) {
    type <- match.arg(type)
    # Order of checklist:
    # islist = 0, rowGene = 0, colK = 0, rowK = 0, colCell = 0)
    if (type == "W") checklist <- c(0, 1, 1, 0, 0)
    if (type == "H") checklist <- c(1, 0, 0, 1, 1)
    if (type == "V") checklist <- c(1, 1, 1, 0, 0)
    if (checklist[1]) {
        if (!is.list(m))
            stop("`", type, ".init` should be a list of matrices")
        if (length(m) != length(nCells))
            stop("Number of matrices in `", type, ".init` should match number",
                 " of datasets in `object`")
        isMat <- sapply(m, is.matrix)
        if (!all(isMat)) {
            stop(sum(!isMat), " elements in `", type, ".init` is not a matrix.")
        }
        isValid <- sapply(seq_along(m), function(i)
            .checkInit.mat(m[[i]], nCells[i], nGenes, k, checklist))
        if (!all(isValid))
            stop("Not all matrices in `", type,
                 ".init` has valid dimensionality.")
    } else {
        if (!is.matrix(m))
            stop("`", type, ".init` should be a matrix")
        if (!.checkInit.mat(m, sum(nCells), nGenes, k, checklist))
            stop("`", type, ".init` does not have valid dimensionality.")
    }
    m
}

.checkInit.mat <- function(m, nCells, nGenes, k, checklist) {
    if (checklist[2]) if (nrow(m) != nGenes) return(FALSE)
    if (checklist[3]) if (ncol(m) != k) return(FALSE)
    if (checklist[4]) if (nrow(m) != k) return(FALSE)
    if (checklist[5]) if (ncol(m) != nCells) return(FALSE)
    return(TRUE)
}

.checkValidFactorResult <- function(object, useDatasets = NULL) {
    result <- TRUE
    useDatasets <- .checkUseDatasets(object, useDatasets)
    if (is.null(object@W)) {
        warning("W matrix does not exist.")
        result <- FALSE
    } else {
        # nGenes <- nrow(object@W)
        k <- ncol(object@W)

        for (d in useDatasets) {
            ld <- dataset(object, d)
            nCells <- ncol(ld)
            if (is.null(ld@V)) {
                warning("V matrix does not exist for dataset '", d, "'.")
                result <- FALSE
            } else {
                if (!identical(dim(ld@V), dim(object@W))) {
                    warning("Dimensionality of V matrix for dataset '", d,
                            "' does not match with W matrix.")
                    result <- FALSE
                }
            }
            if (is.null(ld@H)) {
                warning("H matrix does not exist for dataset '", d, "'.")
                result <- FALSE
            } else {
                if (!identical(dim(ld@H), c(k, nCells))) {
                    warning("Dimensionality of H matrix for dataset '", d,
                            "' is not valid")
                    result <- FALSE
                }
            }
        }
        if (k != object@uns$factorization$k)
            warning("Number of factors does not match with object `k` slot. ")
    }
    if (isFALSE(result))
        stop("Cannot detect valid existing factorization result. ",
             "Please run factorization first.")
}

# !!!MaintainerDeveloperNOTE:
#
# When renaming an argument, for example
#   `foo <- function(object, use.raw = FALSE) {}`
# We want to change it to
#   `foo <- function(object, useRaw = FALSE) {}`
#
# Please write in this way:
#   `foo <- function(object, useRaw = FALSE, use.raw = useRaw)`
# and at the beginning of the function definition body, add:
# ```
# .deprecateArgs(list(use.raw = "useRaw",
#                     another.old.arg.name = "correspondingNewArgName"))
# ```
.deprecateArgs <- function(replace = NULL, defunct = NULL) {
    # This retrieves the exact user call
    call <- match.call(definition = sys.function(-1),
                       call = sys.call(-1))
    callArgs <- rlang::call_args(call)
    parentFuncName <- as.list(call)[[1]]
    # This gives access to variable in the function operation environment
    p <- parent.frame()
    for (old in names(replace)) {
        new <- replace[[old]]
        if (old %in% names(callArgs)) {
            # User set old arg in their call
            what <- paste0(parentFuncName, "(", old, ")")
            with <- paste0(parentFuncName, "(", new, ")")
            lifecycle::deprecate_warn("1.2.0", what, with, always = TRUE)
            if (!new %in% names(callArgs)) {
                # If no setting with new argument is found in user call
                p[[new]] <- p[[old]]
            }
        }
    }
    for (old in defunct) {
        if (old %in% names(callArgs)) {
            lifecycle::deprecate_stop("1.2.0",
                                      paste0(parentFuncName, "(", old, ")"))
        }
    }
    rm(list = c(names(replace), defunct), envir = p)
    invisible(NULL)
}

# Used in gg plotting function to check for dependency
# n - Number of points to be shown
.checkRaster <- function(n, raster = NULL) {
    pkgAvail <- requireNamespace("scattermore", quietly = TRUE)
    if (!is.null(raster) && !is.logical(raster)) {
        stop("Please use `NULL` or logical value for `raster`.")
    }
    if (is.null(raster)) {
        # Automatically decide whether to rasterize depending on number of cells
        if (n > 1e5) {
            if (pkgAvail) {
                raster <- TRUE
                .log("NOTE: Points are rasterized as number of cells/nuclei ",
                     "plotted exceeds 100,000.\n",
                     "Use `raster = FALSE` or `raster = TRUE` to force plot ",
                     "in vector form or not.")
            } else {
                raster <- FALSE
                warning("Number of cells/nuclei plotted exceeds 100,000. ",
                        "Rasterizing the scatter plot is recommended but ",
                        "package \"scattermore\" is not available. ")
            }
        } else {
            raster <- FALSE
        }
    } else if (isTRUE(raster)) {
        if (!pkgAvail) {
            stop("Package \"scattermore\" needed for rasterizing the scatter ",
                 "plot. Please install it by command:\n",
                 "BiocManager::install('scattermore')",
                 call. = FALSE)
        }
    }
    return(raster)
}
