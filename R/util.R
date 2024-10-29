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

cli_or <- function(x) cli::cli_vec(x, list("vec-last" = " or "))

.checkObjVersion <- function(object) {
    if (inherits(object, "liger")) {
        if (!is.newLiger(object))
            cli::cli_abort(
                "Old version of liger object is detected. Please run:
                {.code object <- convertOldLiger(object)}"
            )
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
        cli::cli_abort("A liger object is required.")
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
                cli::cli_abort(
                    "Numeric dataset index out of bound. Only {length(object)}
                    dataset{?s} exist.")
            }
            useDatasets <- unique(useDatasets)
            useDatasets <- names(object)[useDatasets]
        } else if (is.logical(useDatasets)) {
            if (length(useDatasets) != length(object)) {
                cli::cli_abort(
                "Logical dataset subscription does not match the number
                of datasets ({length(object)}).")
            }
            useDatasets <- names(object)[useDatasets]
        } else if (is.character(useDatasets)) {
            if (any(!useDatasets %in% names(object))) {
                cli::cli_abort(
                    "Specified dataset name(s) not found:
                    {.val {useDatasets[!useDatasets %in% names(object)]}}"
                )
            }
        } else {
            cli::cli_abort(
                "Please use a proper numeric/logical/character vector to
                 select dataset to use.")
        }
        if (!is.null(modal)) {
            passing <- sapply(useDatasets, function(d) {
                inherits(dataset(object, d), .modalClassDict[[modal]])
            })
            if (!all(passing))
                cli::cli_abort(
                    "Not all specified datasets are of
                    {(.modalClassDict[[modal]])} class:
                    {.val {useDatasets[!passing]}}"
                )
        }
    }
    useDatasets
}

.findDimRedName <- function(
        object,
        name,
        returnFirst = FALSE,
        stopOnNull = TRUE) {
    if (length(names(object@dimReds)) == 0) {
        cli::cli_abort("No {.field dimRed} available")
    }
    if (is.null(name)) {
        if (returnFirst) return(names(object@dimReds)[1])
        if (stopOnNull) {
            cli::cli_abort("No {.field dimRed} name specified.")
        } else {
            return(NULL)
        }
    }
    if (length(name) == 0) {
        if (stopOnNull) {
            cli::cli_abort("No {.field dimRed} name specified.")
        } else {
            return(NULL)
        }
    }
    if (is.character(name) || is.numeric(name)) {
        if (length(name) != 1) {
            cli::cli_abort("Only one {.field dimRed} can be retrieved at a time.")
        }
        if (is.character(name) && !name %in% names(object@dimReds)) {
            cli::cli_abort("Specified {.field dimRed} {.val {name}} does not exist in the object.")
        }
        if (is.numeric(name)) {
            if (name > length(object@dimReds)) {
                cli::cli_abort("Specified {.field dimRed} index {.val {name}} is out of range.")
            }
            name <- names(object@dimReds)[name]
        }
    } else {
        if (length(name) != length(object@dimReds)) {
            cli::cli_abort(
                c("x" = "{.cls logical} {.var name} specification has wrong length.",
                  "i" = "Should be of length {length(x@dimReds)}")
            )
        }
        name <- names(object@dimReds)[name]
        if (length(name) == 0) {
            if (stopOnNull) {
                cli::cli_abort("No {.field dimRed} name specified.")
            } else {
                return(NULL)
            }
        }
    }
    return(name)
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
    df <- cellMeta(object, columns = variables, cellIdx = cellIdx,
                   as.data.frame = TRUE, drop = FALSE)
    if (isTRUE(checkCategorical)) {
        passing <- sapply(variables, function(v) {
            vec <- df[[v]]
            if (!is.null(dim(vec))) return(FALSE)
            if (is.factor(vec)) return(TRUE)
            if (is.character(vec)) {
                if (length(unique(vec)) > 50)
                    cli::cli_alert_warning(
                        "Categorical variable selection `{v}` has more than 50 unique values."
                    )
                return(TRUE)
            }
            if (is.numeric(vec)) {
                if (length(unique(vec)) > 50)
                    cli::cli_alert_warning(
                        "Categorical variable selection `{v}` has more than 50 unique values."
                    )
                return(FALSE)
            }
        })
        passing <- unlist(passing)
        if (!all(passing)) {
            cli::cli_abort(
                "The following selected variables are not considered as
                categorical. Please use something else or try converting
                them to factor class to force passing checks.
                {.val {variables[!passing]}}"
            )
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
            cli::cli_abort(
                "{paramName} not found in object: {notFound}"
            )
        }
        name <- seq(getNumber(object))
        names(name) <- getNames(object)
        idx <- name[idx]
        names(idx) <- NULL
    } else if (is.logical(idx)) {
        if (length(idx) != getNumber(object)) {
            cli::cli_abort(
                "Length of logical {paramName} does not match to number of {orient}s in `object`."
            )
        }
        idx <- which(idx)
    } else if (is.numeric(idx)) {
        if (max(idx) > getNumber(object) || min(idx) < 1) {
            cli::cli_abort(
                "Numeric {paramName} subscription out of bound."
            )
        }
    } else if (is.null(idx)) {
        idx <- seq(getNumber(object))
    } else {
        cli::cli_abort(
            "Please use character, logical or numeric subscription for {paramName}."
        )
    }
    return(idx)
}

.checkLDSlot <- function(object, slot) {
    avail <- c("rawData", "normData", "scaleData")
    if (is.null(slot)) {
        slot <- avail
    } else {
        if (any(!slot %in% avail)) {
            notFound <- slot[!slot %in% avail]
            cli::cli_abort(
                "Specified slot not availalble: {.val {notFound}}.
                Use one or more from \"rawData\", \"normData\" or \"scaleData\""
            )
        }
        if ("rawData" %in% slot && is.null(rawData(object))) {
            cli::cli_abort("`rawData` is not available for use.")
        }
        if ("normData" %in% slot && is.null(normData(object))) {
            cli::cli_abort("`normData` is not available for use.")
        }
        if ("scaleData" %in% slot && is.null(scaleData(object))) {
            cli::cli_abort("`scaleData` is not available for use.")
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
            cli::cli_abort(
                "{.var {type}Init} should be a list of {.cls matrix}."
            )
        if (length(m) != length(nCells))
            cli::cli_abort(
                "Number of matrices in {.var {type}Init} should match number of datasets in {.var object}."
            )
        isMat <- sapply(m, is.matrix)
        if (!all(isMat)) {
            cli::cli_abort("{sum(!isMat)} elements in {.var {type}Init} is not {.cls matrix}.")
        }
        isValid <- sapply(seq_along(m), function(i)
            .checkInit.mat(m[[i]], nCells[i], nGenes, k, checklist))
        if (!all(isValid))
            cli::cli_abort("Not all matrices in {.var {type}Init} has valid dimensionality.")
    } else {
        if (!is.matrix(m))
            cli::cli_abort("{.var {type}Init} should be a {.cls matrix}.")
        if (!.checkInit.mat(m, sum(nCells), nGenes, k, checklist))
            cli::cli_abort("{.var {type}Init} does not have valid dimensionality.")
    }
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
    if (is.null(object@W)) cli::cli_abort("W matrix does not exist.")
    k <- ncol(object@W)

    for (d in useDatasets) {
        ld <- dataset(object, d)
        nCells <- ncol(ld)
        if (isTRUE(checkV)) {
            if (is.null(ld@V)) {
                cli::cli_alert_danger("V matrix does not exist for dataset {.val {d}}.")
                result <- FALSE
            } else {
                if (!identical(dim(ld@V), dim(object@W))) {
                    cli::cli_alert_danger(
                        "Dimensionality of V matrix for dataset {.val {d}} does not match with W matrix."
                    )
                    result <- FALSE
                }
            }
        }
        if (is.null(ld@H)) {
            cli::cli_alert_danger("H matrix does not exist for dataset {.val {d}}.")
            result <- FALSE
        } else {
            if (!identical(dim(ld@H), c(k, nCells))) {
                cli::cli_alert_danger(
                    "Dimensionality of H matrix for dataset {.val {d}} is not valid."
                )
                result <- FALSE
            }
        }
    }
    if (is.null(object@uns$factorization)) {
        rlang::warn("No recorded factorization parameter found in object.",
                    .frequency = "once", .frequency_id = "inmf_param",
                    use_cli_format = TRUE)
    } else {
        if (k != object@uns$factorization$k)
            cli::cli_alert_danger(
                "Number of factors does not match with recorded parameter."
            )
    }

    if (isFALSE(result))
        cli::cli_abort(
            c(x = "Cannot detect valid existing factorization result. ",
              i = "Please run factorization first. Check warnings.")
        )
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
    parentFuncName <- rlang::call_name(call)
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
        cli::cli_abort("Please use `NULL` or logical value for `raster`.")
    }
    if (is.null(raster)) {
        # Automatically decide whether to rasterize depending on number of cells
        if (n > 1e5) {
            if (pkgAvail) {
                raster <- TRUE
                cli::cli_alert_info(
                    "Points are rasterized since number of cells/nuclei exceeds
                    100,000.
                    Use {.code raster = FALSE} or {.code raster = TRUE} to force
                    plot in vectorized form or not."
                )
            } else {
                raster <- FALSE
                warning("Number of cells/nuclei plotted exceeds 100,000. ",
                        "Rasterizing the scatter plot is recommended but ",
                        "package {.pkg scattermore} is not available. ")
            }
        } else {
            raster <- FALSE
        }
    } else if (isTRUE(raster)) {
        if (!pkgAvail) {
            cli::cli_abort(
                "Package {.pkg scattermore} is needed for rasterizing the scatter
                plot. Please install it by command:
                {.code install.packages('scattermore')}"
            )
        }
    }
    return(raster)
}

.checkArgLen <- function(arg, n, repN = TRUE, class = NULL, .stop = TRUE) {
    if (is.null(arg)) return(arg)
    argname <- deparse(substitute(arg))
    if (length(arg) == 1 && isTRUE(repN)) {
        arg <- rep(arg, n)
    }
    if (length(arg) != n) {
        classTxt <- ifelse(is.null(class), "", " ")
        if (isTRUE(.stop))
            if (!is.null(class)) {
                cli::cli_abort(
                    c("{.var {argname}} has to be a length {ifelse(repN, paste0('1 or ', n), n)} object of class {.cls {class}}.",
                      "i" = "length: {length(arg)}; class: {.cls {class(arg)}}")
                )
            } else {
                cli::cli_abort(
                    c("{.var {argname}} has to be a length {ifelse(repN, paste0('1 or ', n), n)} object.",
                      "i" = "length: {length(arg)}; class: {.cls {class(arg)}}")
                )
            }
        else {
            if (!is.null(class)) {
                cli::cli_alert_warning(
                    c("{.var {argname}} has to be a length {ifelse(repN, paste0('1 or ', n), n)} object of class {.cls {class}}.",
                      i = "Using it anyway.")
                )
            } else {
                cli::cli_alert_warning(
                    c("{.var {argname}} has to be a length {ifelse(repN, paste0('1 or ', n), n)} object.",
                      i = "Using it anyway.")
                )
            }
        }
    }
    if (!is.null(class)) {
        allClassCheck <- sapply(class, function(x) methods::is(arg, x))
        if (!any(allClassCheck)) {
            class <- cli::cli_vec(class, list("vec-quote" ))
            if (isTRUE(.stop))
                cli::cli_abort(
                    c("{.var {argname}} has to be of class {.cls {class}}",
                      "i" = "Given class is {.cls {class(arg)}}")
                )
            else {
                cli::cli_alert_warning(
                    c("{.var {argname}} has to be of class {.cls {class}}. Using it anyway.")
                )
            }
        }
    }
    return(arg)
}

.getSeuratData <- function(object, layer, slot, assay = NULL) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        cli::cli_abort(
            "Package {.pkg Seurat} is needed for this function to work.
            Please install it by command:
            {.code install.packages('Seurat')}"
        )
    }
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        cli::cli_abort(
            "Package {.pkg SeuratObject} is needed for this function to work.
            Please install it by command:
            {.code install.packages('SeuratObject')}"
        )
    }
    assayObj <- Seurat::GetAssay(object, assay = assay)
    if (!"layers" %in% methods::slotNames(assayObj)) {
        # An old version Assay object,
        layer <- slot
    }
    if (utils::packageVersion("SeuratObject") >= package_version("4.9.9")) {
        layers <- SeuratObject::Layers(object, assay = assay, search = layer)
        if (length(layers) == 0) {
            cli::cli_abort("Layer {.val {layer}} not found in object.")
        } else if (length(layers) == 1) {
            data <- SeuratObject::LayerData(object, assay = assay, layer = layers)
        } else {
            data <- lapply(layers, function(l) {
                SeuratObject::LayerData(object, assay = assay, layer = l)
            })
            names(data) <- layers
        }
    } else { # nocov start
        cli::cli_alert_info("Using old Seurat package. Upgrade is recommended.")
        data <- SeuratObject::GetAssayData(object, assay = assay, slot = layer)
    } # nocov end
    return(data)
}

.setSeuratData <- function(object, layer, save, slot, value, assay = NULL,
                           denseIfNeeded = FALSE) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        cli::cli_abort(
            "Package {.pkg Seurat} is needed for this function to work.
            Please install it by command:
            {.emph install.packages('Seurat')}"
        )
    }
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        cli::cli_abort(
            "Package {.pkg SeuratObject} is needed for this function to work.
            Please install it by command:
            {.emph install.packages('SeuratObject')}"
        )
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
    } else { # nocov start
        cli::cli_alert_info("Using old {.pkg Seurat} package. Upgrade is recommended.")
        object <- SeuratObject::SetAssayData(object, assay = assay,
                                             slot = save, new.data = value)
    } # nocov end
    return(object)
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
        return(xsub)
    })
    names(matList) <- levels(y)
    return(matList)
}

searchH <- function(object, useRaw = NULL) {
    if (is.null(useRaw)) {
        # By default, look for quantile-normed H
        H <- getMatrix(object, "H.norm")
        if (is.null(H)) {
            # If not found, look for raw H
            Ht <- Reduce(cbind, getMatrix(object, "H"))
            if (is.null(Ht)) {
                cli::cli_abort(
                    "No cell factor loading available.
                    Please run {.fn runIntegration} and {.fn alignFactors} first."
                )
            } else {
                useRaw <- TRUE
                H <- t(Ht)
            }
        } else {
            useRaw <- FALSE
        }
    } else {
        if (isTRUE(useRaw)) {
            Ht <- Reduce(cbind, getMatrix(object, "H"))
            if (is.null(Ht)) {
                cli::cli_abort(
                    "Raw cell factor loading requested but not found.
                    Please run {.fn runIntegration}."
                )
            } else {
                H <- t(Ht)
            }
        } else {
            H <- getMatrix(object, "H.norm")
            if (is.null(H)) {
                cli::cli_abort(
                    "Aligned cell factor loading requested but
                    not found. Please run {.fn alignFactors} after
                    {.fn runIntegration}."
                )
            }
            useRaw <- FALSE
        }
    }
    return(list(H = H, useRaw = useRaw))
}


.pivot_longer <- function(
        data,
        cols,
        names_to = "name",
        values_to = "value"
) {
    if (is.numeric(cols) || is.logical(cols)) cols <- colnames(data)[cols]
    if (!is.character(cols))
        cli::cli_abort("`cols` should be a character vector.")
    keeps <- setdiff(colnames(data), cols)
    blocks <- lapply(cols, function(col) {
        len <- nrow(data)
        blockData <- list()
        blockData[[names_to]] <- rep(col, len)
        blockData[[values_to]] <- data[[col]]
        for (keep in keeps) {
            blockData[[keep]] <- data[[keep]]
        }
        as.data.frame(blockData)
    })
    do.call(rbind, blocks)
}
