#' Write liger object to H5AD files
#' @description
#' Create an H5AD file from a \linkS4class{liger} object. This function writes
#' only raw counts to \code{adata.X}, while normalized and scaled expression
#' data will not be written, because LIGER use different normalization and
#' scaling strategy than most of the other tools utilizing H5AD format.
#'
#' Supports for single sparse matrices or internal \linkS4class{ligerDataset}
#' objects are also provided if there is a need to convert single datasets.
#' @param object One of \linkS4class{liger}, \linkS4class{ligerDataset} or
#' \link[Matrix]{dgCMatrix-class} object.
#' @param filename A character string, the path to the H5AD file to be written
#' @param obs External data.frame that contains metadata of the cells but does
#' not embed inside the object. Rownames must be identical to the colnames of
#' object.
#' @param var External data.frame that contains metadata of the features but
#' does not embed inside the object. Rownames must be identical to the rownames
#' of object.
#' @param overwrite Logical, whether to overwrite the file if it exists.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param ... Arguments passed down to S3 methods
#' @return No return value, an H5AD file is written to disk with the following
#' specification, assuming the file is loaded to \code{adata} in Python:
#' \itemize{
#'   \item{\code{adata.X} - Raw count CSR matrix, outer joined with all
#'   datasets}
#'   \item{\code{adata.obs} - Cell metadata, with exactly same content of
#'   \code{cellMeta(object)}}
#'   \item{\code{adata.var} - Feature metadata containing only the feature names
#'   as the index of \code{pd.DataFrame}.}
#'   \item{\code{adata.obsm['X_inmf_aligned']} - The integrated embedding,
#'   aligned cell factor loading matrix, the primary output of LIGER, if
#'   available.}
#'   \item{\code{adata.obsm['X_inmf']} - The raw cell factor loading matrix, if
#'   available.}
#'   \item{\code{adata.obsm['<dimRedName>']} - The dimensional reduction matrix,
#'   such as UMAP or TSNE, if available.}
#'   \item{\code{adata.uns['inmf']['W']} - The shared factor feature loading
#'   matrix, if available.}
#'   \item{\code{adata.uns['inmf']['V']['<datasetName>']} - The dataset-specific
#'   factor feature loading matrix, if available.}
#'   \item{\code{adata.uns['inmf']['features']} - The variable features being
#'   used for factorization, supposed to match to the second shape of W and V,
#'   if available.}
#'   \item{\code{adata.uns['inmf']['lambda']} - The hyperparameter lambda used,
#'   the regularization parameter for the factorization, if available.}
#'   \item{\code{adata.uns['inmf']['k']} - The number of factors used for the
#'   factorization, if available.}
#' }
#' @export
#' @examples
#' writeH5AD(pbmc, filename = tempfile(fileext = ".h5ad"))
writeH5AD <- function(
        object,
        ...
) {
    UseMethod('writeH5AD', object)
}

#' @export
#' @rdname writeH5AD
#' @method writeH5AD dgCMatrix
writeH5AD.dgCMatrix <- function(
        object,
        filename,
        obs = NULL,
        var = NULL,
        overwrite = FALSE,
        verbose = getOption('ligerVerbose', TRUE),
        ...
) {
    if (file.exists(filename)) {
        if (isTRUE(overwrite)) {
            file.remove(filename)
        } else {
            cli::cli_abort("H5AD file exists at {.file {normalizePath(filename)}}")
        }
    }
    if (!is.null(obs) &&
        !identical(rownames(obs), colnames(object))) {
        cli::cli_abort("The rownames of {.field obs} should be identical to the colnames of {.field object}.")
    }
    if (!is.null(var) &&
        !identical(rownames(var), rownames(object))) {
        cli::cli_abort("The rownames of {.field var} should be identical to the rownames of {.field object}.")
    }
    dfile <- hdf5r::H5File$new(filename = filename, mode = 'w')
    # Work on expression matrices
    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Adding `adata.X`")
    .writeMatrixToH5AD(x = object, dfile = dfile, dname = "X")
    if (isTRUE(verbose)) cli::cli_process_done(cliID)

    # Add cell metadata
    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Adding `adata.obs`")
    obsDF <- obs %||% data.frame(row.names = colnames(object))
    .writeDataFrameToH5AD(df = obsDF, dfile = dfile, dname = "obs")
    if (isTRUE(verbose)) cli::cli_process_done(cliID)

    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Adding `adata.var`")
    varDF <- var %||% data.frame(row.names = rownames(object))
    .writeDataFrameToH5AD(df = varDF, dfile = dfile, dname = "var")
    if (isTRUE(verbose)) cli::cli_process_done(cliID)

    dfile$close_all()
    return(invisible(NULL))
}

#' @export
#' @rdname writeH5AD
#' @method writeH5AD ligerDataset
writeH5AD.ligerDataset <- function(
        object,
        filename,
        obs = NULL,
        overwrite = FALSE,
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    if (file.exists(filename)) {
        if (isTRUE(overwrite)) {
            file.remove(filename)
        } else {
            cli::cli_abort("H5AD file exists at {.file {normalizePath(filename)}}")
        }
    }
    if (!is.null(obs) &&
        !identical(rownames(obs), colnames(object))) {
        cli::cli_abort("The rownames of {.field obs} should be identical to the colnames of {.field object}.")
    }
    dfile <- hdf5r::H5File$new(filename = filename, mode = 'w')
    # Work on expression matrices
    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Adding `adata.X`")
    .writeMatrixToH5AD(x = rawData(object), dfile = dfile, dname = "X")
    if (isTRUE(verbose)) cli::cli_process_done(cliID)

    # Add cell metadata
    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Adding `adata.obs`")
    obsDF <- obs %||% data.frame(row.names = colnames(object))
    .writeDataFrameToH5AD(df = obsDF, dfile = dfile, dname = "obs")
    if (isTRUE(verbose)) cli::cli_process_done(cliID)

    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Adding `adata.var`")
    varDF <- as.data.frame(featureMeta(object))
    .writeDataFrameToH5AD(df = varDF, dfile = dfile, dname = "var")
    if (isTRUE(verbose)) cli::cli_process_done(cliID)

    dfile$close_all()
    return(invisible(NULL))
}

#' @export
#' @rdname writeH5AD
#' @method writeH5AD liger
writeH5AD.liger <- function(
        object,
        filename,
        overwrite = FALSE,
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    if (file.exists(filename)) {
        if (isTRUE(overwrite)) {
            file.remove(filename)
        } else {
            cli::cli_abort("H5AD file exists at {.file {normalizePath(filename)}}")
        }
    }
    dfile <- hdf5r::H5File$new(filename = filename, mode = 'w')

    # Work on expression matrices
    if (isTRUE(verbose)) {
        cliID <- cli::cli_process_start("Adding merged raw counts to `adata.X`")
    }
    X <- mergeSparseAll(rawData(object))
    .writeMatrixToH5AD(x = X, dfile = dfile, dname = "X")
    if (isTRUE(verbose)) cli::cli_process_done(cliID)

    # Add cell metadata
    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Adding `adata.obs`")
    .writeDataFrameToH5AD(
        df = cellMeta(object, as.data.frame = TRUE),
        dfile = dfile,
        dname = "obs"
    )
    if (isTRUE(verbose)) cli::cli_process_done(cliID)

    # Work on var, as per liger object design, we don't provide feature metadata
    # for merged dataset for now.
    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Adding empty `adata.var`")
    varDF <- data.frame(
        row.names = rownames(X)
    )
    .writeDataFrameToH5AD(df = varDF, dfile = dfile, dname = "var")
    if (isTRUE(verbose)) cli::cli_process_done(cliID)

    # Add dimensional reduction information
    # Check if any to be written first
    obsmList <- list()
    if (!is.null(getMatrix(object, "H.norm"))) {
        if (isTRUE(verbose)) {
            cli::cli_alert_info("Planning aligned cell factor loading to `adata.obsm['X_inmf_aligned']`")
        }
        # Transposed due to H5AD automatically aligns column-majored storage in
        # HDF5 to row-majored arrays in NumPy array. i.e. an inner array for
        # a cell (obs entry) stays the same.
        obsmList$X_inmf_aligned <- t(getMatrix(object, "H.norm"))
    }

    if (!any(sapply(getMatrix(object, "H"), is.null))) {
        if (isTRUE(verbose)) {
            cli::cli_alert_info("Planning raw cell factor loading to `adata.obsm['X_inmf']`")
        }
        HList <- getMatrix(object, 'H')
        obsmList$X_inmf <- Reduce(cbind, HList)
    }
    if (length(dimReds(object)) > 0) {
        for (i in seq_along(dimReds(object))) {
            dimRedName <- names(dimReds(object))[i]
            if (isTRUE(verbose)) {
                cli::cli_alert_info(
                    "Planning `dimRed(object, {dimRedName})` to `adata.obsm['{dimRedName}']`"
                )
            }
            # Transposed due to the same reason as H.norm above.
            obsmList[[dimRedName]] <- t(dimReds(object)[[i]])
        }
    }
    if (length(obsmList) == 0 && isTRUE(verbose))
        cli::cli_alert_info("No low-dimensional representation to be written.")
    if (length(obsmList) > 0) {
        if (isTRUE(verbose))
            cliID <- cli::cli_process_start("Adding obsm")
        .createMapping(dfile, 'obsm')
        for (i in seq_along(obsmList)) {
            obsmName <- names(obsmList)[i]
            dfile[['obsm']]$create_dataset(
                name = obsmName,
                robj = obsmList[[i]],
                dtype = .H5AD.guessDType(obsmList[[i]]),
                chunk_dims = NULL
            )
        }
        if (isTRUE(verbose)) cli::cli_process_done(cliID)
    }


    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Adding uns")
    .createMapping(dfile, 'uns')
    .createMapping(dfile[['uns']], 'inmf')
    unsList <- object@uns$factorization
    if (!is.null(getMatrix(object, 'W'))) {
        unsList$W <- getMatrix(object, 'W')
        unsList$features <- rownames(unsList$W)
    }
    if (length(unsList) > 0) {
        for (i in seq_along(unsList)) {
            unsName <- names(unsList)[i]
            dfile[['uns/inmf']]$create_dataset(
                name = unsName,
                robj = unsList[[i]],
                dtype = .H5AD.guessDType(unsList[[i]]),
                chunk_dims = NULL
            )
        }
    }
    Vlist <- getMatrix(object, 'V')
    Vlist <- Vlist[!sapply(Vlist, is.null)]
    if (length(Vlist) > 0) {
        .createMapping(dfile[['uns/inmf']], 'V')
        for (i in seq_along(Vlist)) {
            VName <- names(Vlist)[i]
            dfile[['uns/inmf/V']]$create_dataset(
                name = VName,
                robj = Vlist[[i]],
                dtype = .H5AD.guessDType(Vlist[[i]]),
                chunk_dims = NULL
            )
        }
    }
    if (isTRUE(verbose)) cli::cli_process_done(cliID)


    dfile$close_all()
    return(invisible(NULL))
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
        attr_name = 'encoding-type',
        robj = 'csr_matrix',
        dtype = .H5AD.guessDType('csr_matrix'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-version',
        robj = '0.1.0',
        dtype = .H5AD.guessDType('0.1.0'),
        space = hdf5r::H5S$new(type = "scalar")
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
        path <- paste0(dname, '/', i)
        if (is.factor(df[[i]])) {
            # Writing categorical
            .writeCategroricalToH5AD(df[[i]], dfile, path)
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
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-type',
        robj = 'dataframe',
        dtype = .H5AD.guessDType('dataframe'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-version',
        robj = '0.2.0',
        dtype = .H5AD.guessDType('0.2.0'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    return(invisible(NULL))
}

.writeCategroricalToH5AD <- function(x, dfile, dname) {
    dfile$create_group(name = dname)
    dfile[[dname]]$create_dataset(
        name = 'categories',
        robj = levels(x),
        dtype = .H5AD.guessDType(levels(x)),
        chunk_dims = NULL
    )
    dfile[[dname]]$create_dataset(
        name = "codes",
        robj = as.integer(x) - 1L,
        dtype = .H5AD.guessDType(as.integer(x) - 1L)
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-type',
        robj = 'categorical',
        dtype = .H5AD.guessDType('categorical'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-version',
        robj = '0.2.0',
        dtype = .H5AD.guessDType('0.2.0'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    dfile[[dname]]$create_attr(
        attr_name = 'ordered',
        robj = FALSE,
        dtype = .H5AD.guessDType(FALSE),
        space = hdf5r::H5S$new(type = "scalar")
    )
    return(invisible(NULL))
}

.createMapping <- function(dfile, dname) {
    dfile$create_group(name = dname)
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-type',
        robj = 'dict',
        dtype = .H5AD.guessDType('dict'),
        space = hdf5r::H5S$new(type = "scalar")
    )
    dfile[[dname]]$create_attr(
        attr_name = 'encoding-version',
        robj = '0.1.0',
        dtype = .H5AD.guessDType('0.1.0'),
        space = hdf5r::H5S$new(type = "scalar")
    )
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

.boolToInt <- function(x) {
    x <- as.integer(x)
    x[is.na(x)] <- 2L
    return(x)
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
