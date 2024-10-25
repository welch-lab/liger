help.ligerToH5AD <- function() {
    cat(
        "Write liger object to H5AD files\n",
        "\n",
        "Usage:\n",
        "  ligerToH5AD(\n",
        "    object,\n",
        "    filename,\n",
        "    useSlot = c(\"scaleData\", \"normData\", \"rawData\"),\n",
        "    saveRaw = useSlot != \"rawData\",\n",
        "    overwrite = FALSE,\n",
        "    verbose = getOption(\"ligerVerbose\", TRUE)\n",
        "  )\n",
        "\n",
        "Arguments:\n",
        "  object:    A liger object\n",
        "  filename:  A character string, the path to the H5AD file to be written\n",
        "  useSlot:   Character scalar, which type of data is going to be stored\n",
        "             to `adata.X`. Default \"scaleData\", choose from\n",
        "             \"scaleData\", \"normData\", or \"rawData\".\n",
        "  saveRaw:   Logical, whether to save rawData to `adata.raw.X`.\n",
        "             Default TRUE when `useSlot` is not \"rawData\", otherwise\n",
        "             FALSE.\n",
        "  overwrite: Logical, whether to overwrite the file if it exists.\n",
        "  verbose:   Logical. Whether to show information of the progress. Default\n",
        "             `getOption(\"ligerVerbose\")` which is TRUE if users have not set.\n",
        "\n",
        "Value:\n",
        "  No return value, an H5AD file is written to disk\n",
        "\n",
        "Examples:\n",
        "  ligerToH5AD(pbmc, filename = tempfile(fileext = \".h5ad\"))\n"
    )
}


#' Write liger object to H5AD files
#' @param object A \linkS4class{liger} object
#' @param filename A character string, the path to the H5AD file to be written
#' @param useSlot Character scalar, which type of data is going to be stored
#' to \code{adata.X}. Default \code{"scaleData"}, choose from
#' \code{"scaleData"}, \code{"normData"}, or \code{"rawData"}.
#' @param saveRaw Logical, whether to save rawData to \code{adata.raw.X}.
#' Default \code{TRUE} when \code{useSlot} is not \code{"rawData"}, otherwise
#' \code{FALSE}.
#' @param overwrite Logical, whether to overwrite the file if it exists.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @return No return value, an H5AD file is written to disk
#' @export
#' @examples
#' ligerToH5AD(pbmc, filename = tempfile(fileext = ".h5ad"))
ligerToH5AD <- function(
        object,
        filename,
        useSlot = c("scaleData", "normData", "rawData"),
        saveRaw = useSlot != "rawData",
        overwrite = FALSE,
        verbose = getOption("ligerVerbose", TRUE)
) {
    if (file.exists(filename)) {
        if (isTRUE(overwrite)) {
            file.remove(filename)
        } else {
            cli::cli_abort("H5AD file exists at {.file {normalizePath(filename)}}")
        }
    }
    useSlotCheckOrder <- c("scaleData", "normData", "rawData")
    useSlot <- match.arg(useSlot)
    rownames <- "_index"
    dfile <- hdf5r::H5File$new(
        filename = filename,
        mode = ifelse(test = overwrite, yes = "w", no = "w-")
    )

    # Work on expression matrices
    Xslot <- NULL
    if (all(!sapply(rliger::getMatrix(object, useSlot), is.null))) {
        # If the specified slot is all available, use it
        Xslot <- useSlot
        Xmat <- rliger::mergeSparseAll(rliger::getMatrix(object, useSlot))
    } else {
        cli::cli_alert_warning(
            "Requested {.field {useSlot}} is not available, trying other slots."
        )
        # Otherwise, check by priority order
        for (i in useSlotCheckOrder) {
            if (all(!sapply(getMatrix(object, i), is.null))) {
                # If the specified slot is all available, use it
                Xslot <- i
                Xmat <- rliger::mergeSparseAll(rliger::getMatrix(object, i))
                break
            }
        }
    }
    if (is.null(Xslot)) {
        cli::cli_abort("No data available to be added to {.field X}")
    }

    if (isTRUE(saveRaw) && Xslot != "rawData") {
        rawSlot <- "rawData"
    } else {
        rawSlot <- NULL
    }

    if (isTRUE(verbose)) {
        cliID <- cli::cli_process_start(sprintf("Adding %s to X", Xslot))
    }
    .writeMatrixToH5AD(x = Xmat, dfile = dfile, dname = "X")
    .H5AD.addEncoding(dfile = dfile, dname = 'X')
    if (isTRUE(verbose)) cli::cli_process_done(cliID)

    # Work on var, as per liger object design, we don't provide feature metadata
    # for merged dataset for now.
    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Adding var")
    varDF <- data.frame(name = rownames(Xmat), row.names = rownames(Xmat))
    .writeDataFrameToH5AD(df = varDF, dfile = dfile, dname = "var")
    if (isTRUE(verbose)) cli::cli_process_done(cliID)

    # Add raw
    if (!is.null(rawSlot) && isTRUE(saveRaw)) {
        if (isTRUE(verbose)) {
            cliID <- cli::cli_process_start(sprintf("Adding %s to raw", rawSlot))
        }
        rawMat <- rliger::mergeSparseAll(rliger::rawData(object))
        dfile$create_group(name = 'raw')
        .writeMatrixToH5AD(x = rawMat, dfile = dfile, dname = "raw/X")
        .H5AD.addEncoding(dfile = dfile, dname = 'raw/X')

        # Similarly, Add meta.features
        rawVarDF <- data.frame(name = rownames(rawMat),
                               row.names = rownames(rawMat))
        .writeDataFrameToH5AD(df = rawVarDF, dfile = dfile, dname = "raw/var")
    }

    # Add cell metadata
    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Adding obs")
    .writeDataFrameToH5AD(
        df = cellMeta(object, as.data.frame = TRUE),
        dfile = dfile,
        dname = "obs"
    )
    if (isTRUE(verbose)) cli::cli_process_done(cliID)


    # Add dimensional reduction information
    obsm <- dfile$create_group(name = 'obsm')
    varm <- dfile$create_group(name = 'varm')
    reductions <- object@dimReds
    h <- rliger::getMatrix(object, "H")
    if (all(!sapply(h, is.null))) {
        H <- Reduce(cbind, h)
        H <- t(H)
        reductions[["H"]] <- H
    }
    if (!is.null(object@H.norm)) reductions[["H_norm"]] <- object@H.norm
    for (reduc in names(reductions)) {
        newname <- paste0("X_", tolower(reduc))
        if (isTRUE(verbose)) {
            cliID <- cli::cli_process_start(
                sprintf("Adding low-dim representation %s to obsm as %s",
                        reduc, newname)
            )
        }
        obsm$create_dataset(
            name = newname,
            robj = t(reductions[[reduc]]),
            dtype = .H5AD.guessDType(reductions[[reduc]]),
            chunk_dims = NULL
        )
        if (isTRUE(verbose)) cli::cli_process_done(cliID)
    }

    if (!is.null(object@W) &&
        Xslot == "scaleData") {
        W <- object@W
        if (!identical(rownames(W), rownames(Xmat))) {
            cli::cli_alert_info(
                "Feature names in {.field W} do not match those in 'X' ({.field {Xslot}}), skipping."
            )
        }
        if (isTRUE(verbose)) {
            cliID <- cli::cli_process_start("Adding factor feature loadings to varm as 'W'")
        }
        varm$create_dataset(
            name = "W",
            robj = t(W),
            dtype = .H5AD.guessDType(t(W)),
            chunk_dims = NULL
        )
        if (isTRUE(verbose)) cli::cli_process_done(cliID)
    }
    dfile$close_all()
}

.writeMatrixToH5AD <- function(
        x,
        dfile,
        dname
) {
    dnamePaths <- unlist(strsplit(dname, '/'))
    # Recursively create groups
    for (i in seq_along(dnamePaths)) {
        search <- paste0(dnamePaths[1:i], collapse = '/')
        if (!dfile$exists(name = search)) {
            dfile$create_group(name = search)
        }
    }
    dfile[[dname]]$create_dataset(
        name = 'data',
        robj = x@x,
        dtype = .H5AD.guessDType(x@x)
    )
    dfile[[dname]]$create_dataset(
        name = 'indices',
        robj = x@i,
        dtype = .H5AD.guessDType(x@i)
    )
    dfile[[dname]]$create_dataset(
        name = 'indptr',
        robj = x@p,
        dtype = .H5AD.guessDType(x@p)
    )
    dfile[[dname]]$create_attr(
        attr_name = 'shape',
        robj = rev(dim(x)),
        dtype = .H5AD.guessDType(dim(x))
    )
    return(invisible(NULL))
}

.writeDataFrameToH5AD <- function(
        df,
        dfile,
        dname
) {
    dnamePaths <- unlist(strsplit(dname, '/'))
    # Recursively create groups
    for (i in seq_along(dnamePaths)) {
        search <- paste0(dnamePaths[1:i], collapse = '/')
        if (!dfile$exists(name = search)) {
            dfile$create_group(name = search)
        }
    }
    # Add index
    dfile[[dname]]$create_dataset(
        name = '_index',
        robj = rownames(df),
        dtype = .H5AD.guessDType(rownames(df))
    )
    dfile[[dname]]$create_attr(
        attr_name = '_index',
        robj = '_index',
        dtype = .H5AD.guessDType('index'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    # Add columns
    for (i in colnames(df)) {
        if (is.factor(df[[i]])) {
            dfile[[dname]]$create_dataset(
                name = i,
                robj = as.integer(df[[i]]) - 1L,
                dtype = .H5AD.guessDType(as.integer(df[[i]]) - 1L)
            )
            if (!dfile[[dname]]$exists(name = '__categories')) {
                dfile[[dname]]$create_group(name = '__categories')
            }
            dfile[[dname]][['__categories']]$create_dataset(
                name = i,
                robj = levels(df[[i]]),
                dtype = .H5AD.guessDType(levels(df[[i]])),
                chunk_dims = NULL
            )
            dfile[[dname]][['__categories']][[i]]$create_attr(
                attr_name = "ordered",
                robj = FALSE,
                dtype = .H5AD.guessDType(FALSE)
            )
            ref <- .H5.create_reference(dfile[[dname]][['__categories']][[i]])
            dfile[[dname]][[i]]$create_attr(
                attr_name = 'categories',
                robj = ref,
                dtype = .H5AD.guessDType(ref),
                space = hdf5r::H5S$new(type = "scalar")
            )
        } else {
            dfile[[dname]]$create_dataset(
                name = i,
                robj = df[[i]],
                dtype = .H5AD.guessDType(df[[i]])
            )
        }
    }
    dfile[[dname]]$create_attr(
        attr_name = 'column-order',
        robj = colnames(df),
        dtype = .H5AD.guessDType(colnames(df))
    )
    encoding.info <- c('type' = 'dataframe', 'version' = '0.1.0')
    names(encoding.info) <- paste0('encoding-', names(encoding.info))
    for (i in seq_along(encoding.info)) {
        attr.name <- names(encoding.info)[i]
        attr.value <- encoding.info[i]
        if (dfile[[dname]]$attr_exists(attr_name = attr.name)) {
            dfile[[dname]]$attr_delete(attr_name = attr.name)
        }
        dfile[[dname]]$create_attr(
            attr_name = attr.name,
            robj = attr.value,
            dtype = .H5AD.guessDType(attr.value),
            space = hdf5r::H5S$new(type = "scalar")
        )
    }
    return(invisible(NULL))
}

# self - Only support an H5D object for now
.H5.create_reference <- function(self, ...) {
    space <- self$get_space()
    do.call("[", c(list(space), list(...)))
    ref_type <- hdf5r::h5const$H5R_OBJECT
    ref_obj <- hdf5r::H5R_OBJECT$new(1, self)
    res <- .Call("R_H5Rcreate", ref_obj$ref, self$id, ".", ref_type,
                 space$id, FALSE, PACKAGE = "hdf5r")
    if (res$return_val < 0) {
        stop("Error creating object reference")
    }
    ref_obj$ref <- res$ref
    return(ref_obj)
}

.H5AD.addEncoding <- function(dfile, dname) {
    encoding.info <- c('type' = 'csr_matrix', 'version' = '0.1.0')
    names(encoding.info) <- paste0('encoding-', names(encoding.info))
    if (inherits(dfile[[dname]], 'H5Group')) {
        for (i in seq_along(encoding.info)) {
            attr.name <- names(encoding.info)[i]
            attr.value <- encoding.info[i]
            if (dfile[[dname]]$attr_exists(attr_name = attr.name)) {
                dfile[[dname]]$attr_delete(attr_name = attr.name)
            }
            dfile[[dname]]$create_attr(
                attr_name = attr.name,
                robj = attr.value,
                dtype = .H5AD.guessDType(x = attr.value),
                space = hdf5r::H5S$new(type = "scalar")
            )
        }
    }
    return(invisible(NULL))
}

.H5AD.guessDType <- function(x, stype = "utf8", ...) {
    dtype <- hdf5r::guess_dtype(x = x, ...)
    if (inherits(dtype, "H5T_STRING")) {
        dtype <- .H5AD.stringType(stype = stype)
    }
    else if (inherits(dtype, "H5T_COMPOUND")) {
        cpd.dtypes <- dtype$get_cpd_types()
        for (i in seq_along(cpd.dtypes)) {
            if (inherits(cpd.dtypes[[i]], "H5T_STRING")) {
                cpd.dtypes[[i]] <- .H5AD.stringType(stype = stype)
            }
        }
        dtype <- hdf5r::H5T_COMPOUND$new(
            labels = dtype$get_cpd_labels(),
            dtypes = cpd.dtypes,
            size = dtype$get_size()
        )
    }
    else if (inherits(dtype, "H5T_LOGICAL")) {
        dtype <- hdf5r::guess_dtype(x = .boolToInt(x = x), ...)
    }
    return(dtype)
}

.H5AD.stringType <- function(stype = c("utf8", "ascii7"))
{
    stype <- match.arg(arg = stype)
    switch(
        EXPR = stype,
        utf8 = hdf5r::H5T_STRING$new(size = Inf)$set_cset(
            cset = hdf5r::h5const$H5T_CSET_UTF8
        ),
        ascii7 = hdf5r::H5T_STRING$new(size = 7L)
    )
}

.boolToInt <- function(x) {
    x <- as.integer(x)
    x[is.na(x)] <- 2L
    return(x)
}

