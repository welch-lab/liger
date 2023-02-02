.log <- function(..., level = 1) {
    line <- paste(rep("...", level), collapse = "")
    pref <- paste0(date(), " ", line, " ")
    indent <- paste(rep(" ", nchar(pref)), collapse = "")
    content <- list(...)
    msg <- paste(content, collapse = "")
    msg <- gsub("\n", paste0("\n", indent), msg)
    message(pref, msg)
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
.checkUseDatasets <- function(object, useDatasets = NULL) {
    if (!inherits(object, "liger"))
        stop("A liger object is required.")
    if (is.null(useDatasets)) useDatasets <- names(object)
    else {
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
    }
    useDatasets
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
    avail <- c("raw.data", "norm.data", "scale.data")
    if (is.null(slot)) slot <- avail
    else {
        if (any(!slot %in% avail)) {
            notFound <- slot[!slot %in% avail]
            stop("Specified slot not availalble: ",
                 paste(notFound, collapse = ", "), ". Use one or more from ",
                 '"raw.data", "norm.data" or "scale.data"')
        }
        if ("raw.data" %in% slot &
            is.null(raw.data(object))) {
            stop("`raw.data` is not available for use.")
        }
        if ("norm.data" %in% slot &
            is.null(norm.data(object))) {
            stop("`norm.data` is not available for use.")
        }
        if ("scale.data" %in% slot &
            is.null(scale.data(object))) {
            stop("`scale.data` is not available for use.")
        }
    }
    slot
}

.checkInit <- function(m, nCells, nGenes, k, type = c("W", "H", "V", "A", "V")) {
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

.checkValidFactorResult <- function(object) {
    result <- TRUE
    if (is.null(object@W)) {
        warning("W matrix does not exist.")
        result <- FALSE
    } else {
        nGenes <- nrow(object@W)
        k <- ncol(object@W)
        for (d in names(object)) {
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
#                     another.old.arg.name = "correspondingNewArgName"),
#                call = rlang::call_args(match.call()))
# ```
# NEVER FORGET the `call` argument which provides the access to know whether an
# argument is specified by users or left missing.
.deprecateArgs <- function(replace = NULL, defunct = NULL, call = NULL) {
    parentFuncName <- as.list(sys.call(-1))[[1]]
    p <- parent.frame()
    for (old in names(replace)) {
        new <- replace[[old]]
        if (old %in% names(call)) {
            # User set old arg in their call
            what <- paste0(parentFuncName, "(", old, ")")
            with <- paste0(parentFuncName, "(", new, ")")
            lifecycle::deprecate_warn("1.2.0", what, with, always = TRUE)
            if (!new %in% names(call)) {
                # If no setting with new argument is found in user call
                p[[new]] <- p[[old]]
            }
        }
    }
    for (old in names(defunct))
        lifecycle::deprecate_stop("1.2.0", old)

    invisible(NULL)
}

# Retrieve a vector of feature values with length equal to nCell
.retrieveCellFeatureOld <- function(object, feature, slot) {
    if (slot != "cell.meta") {
        # Should be a list if we are not using H.norm or W
        dataList <- getMatrix(object, slot)
        if (slot != "H.norm")
            exp <- unlist(lapply(dataList, function(m) {
                if (inherits(m, "H5D")) {
                    if (length(m$dims) == 2) {
                        # Dense matrix, e.g. scale.data
                        # 1. Get the numeric index of feature
                        # 2. extract

                    } else if (length(m$dims) == 1) {
                        # Most likely the `x` of a sparse matrix
                        # 1. go through chunks to subscribe

                    }
                }
                else m[feature,]
            }), use.names = FALSE)
        else exp <- dataList[, feature]
        featureVal <- exp
    } else {
        featureVal <- object[[feature]]
    }
}

