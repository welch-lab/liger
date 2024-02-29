.log <- function(..., level = 1) {
    line <- paste(rep("...", level), collapse = "")
    pref <- paste0(date(), " ", line, " ")
    indent <- paste(rep(" ", nchar(pref)), collapse = "")
    content <- list(...)
    content <- lapply(content, as.character)
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
        droplevels = TRUE,
        returnList = FALSE
) {
    df <- cellMeta(object, columns = variables, cellIdx = cellIdx, as.data.frame = TRUE,
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
    for (v in colnames(df)) {
        if (is.factor(df[[v]]) && isTRUE(droplevels)) {
            df[[v]] <- droplevels(df[[v]])
        }
        if (is.character(df[[v]]) && isTRUE(checkCategorical)) {
            df[[v]] <- factor(df[[v]])
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
        if (max(idx) > getNumber(object) || min(idx) < 1) {
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

.checkValidFactorResult <- function(object, useDatasets = NULL,
                                    checkV = TRUE) {
    result <- TRUE
    useDatasets <- .checkUseDatasets(object, useDatasets)
    if (is.null(object@W)) stop("W matrix does not exist.")
    k <- ncol(object@W)

    for (d in useDatasets) {
        ld <- dataset(object, d)
        nCells <- ncol(ld)
        if (isTRUE(checkV)) {
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
    if (isFALSE(result))
        stop("Cannot detect valid existing factorization result. ",
             "Please run factorization first. Check warnings.")
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
                       call = sys.call(-1), envir = parent.frame(2))
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
            lifecycle::deprecate_warn("1.99.0", what, with, always = TRUE)
            if (!new %in% names(callArgs)) {
                # If no setting with new argument is found in user call
                p[[new]] <- p[[old]]
            }
        }
    }
    for (old in defunct) {
        if (old %in% names(callArgs)) {
            lifecycle::deprecate_warn("1.99.0",
                                      paste0(parentFuncName, "(", old, ")"), details = "Ignored.")
            # lifecycle::deprecate_stop("1.99.0",
            #                           paste0(parentFuncName, "(", old, ")"))
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

.checkArgLen <- function(arg, n, repN = TRUE, class = NULL, .stop = TRUE) {
    argname <- deparse(substitute(arg))
    if (!is.null(arg)) {
        if (length(arg) == 1 && isTRUE(repN)) {
            arg <- rep(arg, n)
        }
        if (length(arg) != n) {
            if (isTRUE(.stop))
                stop("`", argname, "` has to be a vector of length ", n)
            else {
                warning("`", argname, "` has to be a vector of length ", n)
            }
        }
    }
    if (!is.null(class)) {
        allClassCheck <- sapply(class, function(x) methods::is(arg, x))
        if (!any(allClassCheck)) {
            if (isTRUE(.stop))
                stop("`", argname, "` has to be of class '",
                     paste(class, collapse = "', '"), "'")
            else {
                warning("`", argname, "` has to be of class '",
                        paste(class, collapse = "', '"), "'")
            }
        }
    }
    return(arg)
}

# Format "not found" string. When we need `need` elements from some source
# `from` format the string of ", " separeted list of not found elements.
.nfstr <- function(need, from) {
    nf <- need[!need %in% from]
    paste(nf, collapse = ", ")
}




.getSeuratData <- function(object, layer, slot, assay = NULL) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat package has to be installed in advance.")
    }
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        stop("SeuratObject package has to be installed in advance.")
    }
    assayObj <- Seurat::GetAssay(object, assay = assay)
    if (!"layers" %in% methods::slotNames(assayObj)) {
        # An old version Assay object,
        layer <- slot
    }
    if (utils::packageVersion("SeuratObject") >= package_version("4.9.9")) {
        layers <- SeuratObject::Layers(object, assay = assay, search = layer)
        if (length(layers) == 0) {
            stop("Layer '", layer, "' not found in object.")
        } else if (length(layers) == 1) {
            data <- SeuratObject::LayerData(object, assay = assay, layer = layers)
        } else {
            data <- lapply(layers, function(l) {
                SeuratObject::LayerData(object, assay = assay, layer = l)
            })
            names(data) <- layers
        }
    } else {
        data <- SeuratObject::GetAssayData(object, assay = assay, slot = layer)
    }
    return(data)
}

.setSeuratData <- function(object, layer, save, slot, value, assay = NULL,
                           denseIfNeeded = FALSE) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat package has to be installed in advance.")
    }
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        stop("SeuratObject package has to be installed in advance.")
    }
    assayObj <- Seurat::GetAssay(object, assay = assay)
    if (!"layers" %in% methods::slotNames(assayObj)) {
        # An old version Assay object,
        save <- slot
        if (isTRUE(denseIfNeeded)) value <- as.matrix(value)
    }
    if (utils::packageVersion("SeuratObject") >= package_version("4.9.9")) {
        if (!is.list(value)) {
            SeuratObject::LayerData(object, assay = assay, layer = save) <- value
        } else {
            # Assume this come from .getSeuratData and is guaranteed to be
            # existing layers
            olayer <- layer
            layer <- names(value)
            if (length(save) != length(layer)) {
                save <- make.unique(gsub(olayer, save, layer))
            }
            for (i in seq_along(layer)) {
                SeuratObject::LayerData(object, assay = assay, layer = save[i]) <- value[[layer[i]]]
            }
        }
    } else {
        object <- SeuratObject::SetAssayData(object, assay = assay,
                                             slot = save, new.data = value)
    }
    return(object)
}

# plyr::mapvalues
mapvalues <- function(x, from, to, warn_missing = TRUE) {
    if (length(from) != length(to)) {
        stop("`from` and `to` vectors are not the same length.")
    }
    if (!is.atomic(x) && !is.null(x)) {
        stop("`x` must be an atomic vector or NULL.")
    }

    if (is.factor(x)) {
        # If x is a factor, call self but operate on the levels
        levels(x) <- mapvalues(levels(x), from, to, warn_missing)
        return(x)
    }

    mapidx <- match(x, from)
    mapidxNA  <- is.na(mapidx)

    # index of items in `from` that were found in `x`
    from_found <- sort(unique(mapidx))
    if (warn_missing && length(from_found) != length(from)) {
        message("The following `from` values were not present in `x`: ",
                paste(from[!(1:length(from) %in% from_found) ], collapse = ", "))
    }

    x[!mapidxNA] <- to[mapidx[!mapidxNA]]
    x
}

.DataFrame.as.data.frame <- function(x)
{
    # Copied from Bioconductor package S4Vectors:::.as.data.frame.DataFrame
    # Removed some lines not necessary for liger
    row.names <- rownames(x)
    if (!is.null(row.names)) row.names <- make.unique(row.names)
    else if (ncol(x) == 0L) row.names <- seq_len(nrow(x))
    x_colnames <- colnames(x)
    df_list <- lapply(stats::setNames(seq_along(x), x_colnames), function(j) {
        col <- x[[j]]
        if (is.data.frame(col)) return(col)
        protect <- !methods::is(col, "AsIs") && is.list(col) && !is.object(col)
        if (protect) col <- I(col)  # set AsIs class to protect column
        # Doing all this copy-paste to avoid the R 4.3.2 new deprecation on the
        # direct call of as.data.frame.<class>. Yet I still don't understand how
        # the warning is invoked. I suspect that it is related to BiocGenerics
        # which has S4 generics for as.data.frame. R base as.data.frame.<class>
        # implements the check by seeing if the function that calls itself is
        # identical to the generic base::as.data.frame. It could be that
        # BiocGenerics::as.data.frame overwrote it in S4Vectors namespace.
        # ALL FOR THE NEXT LINE BELOW.
        df <- as.data.frame(col)
        if (protect) df[[1L]] <- unclass(df[[1L]])  # drop AsIs class
        if (is.null(colnames(col)) && ncol(df) == 1L)
            colnames(df) <- x_colnames[[j]]
        df
    })
    do.call(data.frame,
            c(df_list, list(row.names = row.names, stringsAsFactors = FALSE)))
}


splitRmMiss <- function(x, y) {
    y <- factor(y)
    y <- droplevels(y)
    matList <- lapply(levels(y), function(lvl) {
        idx <- y == lvl
        xsub <- x[, idx, drop = FALSE]
        xsub[rowSums(xsub) > 0, , drop = FALSE]
    })
    names(matList) <- levels(y)
    return(matList)
}
