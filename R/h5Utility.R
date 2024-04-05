#' Apply function to chunks of H5 data in ligerDataset object
#' @description h5 calculation wrapper, that runs specified calculation with
#' on-disk matrix in chunks
#' @param object A \linkS4class{ligerDataset} object.
#' @param FUN A function that is applied to each chunk. See detail for
#' restrictions.
#' @param init Initialized result if it need to be updated iteratively. Default
#' \code{NULL}.
#' @param useData The slot name of the data to be processed. Choose from
#' \code{"rawData"}, \code{"normData"}, \code{"scaleData"}. Default
#' \code{"rawData"}.
#' @param chunkSize Number if columns to be included in each chunk.
#' Default \code{1000}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param ... Other arguments to be passed to \code{FUN}.
#' @details The \code{FUN} function has to have the first four arguments ordered
#' by:
#' \enumerate{
#' \item \bold{chunk data:} A sparse matrix
#' (\code{\link[Matrix]{dgCMatrix-class}}) containing maximum \code{chunkSize}
#' columns.
#' \item \bold{x-vector index:} The index that subscribes the vector of \code{x}
#' slot of a dgCMatrix, which points to the values in each chunk. Mostly used
#' when need to write a new sparse matrix to H5 file.
#' \item \bold{cell index:} The column index of each chunk out of the whole
#' original matrix
#' \item \bold{Initialized result:} A customized object, the value passed to
#' \code{H5Apply(init)} argument will be passed here in the first iteration. And
#' the returned value of \code{FUN} will be iteratively passed here in next
#' chunk iterations. So it is important to keep the object structure of the
#' returned value consistent with \code{init}.
#' }
#' No default value to these four arguments should be pre-defined because
#' \code{H5Apply} will automatically generate the input.
H5Apply <- function(
        object,
        FUN,
        init = NULL,
        useData = c("rawData", "normData"),
        chunkSize = 1000,
        verbose = getOption("ligerVerbose"),
        ...
) {
    fun.args <- list(...)
    useData <- match.arg(useData)
    h5meta <- h5fileInfo(object)
    numCells <- ncol(object)
    numFeatures <- nrow(object)

    prev_end_col <- 1
    prev_end_data <- 1
    numChunks <- ceiling(numCells / chunkSize)
    ind <- 0
    h5file <- h5meta$H5File
    colptr <- h5file[[h5meta$indptrName]]
    rowind <- h5file[[h5meta$indicesName]]
    data <- h5file[[h5meta[[useData]]]]
    if (isTRUE(verbose))
        cliID <- cli::cli_progress_bar(name = "HDF5 chunk processing", type = "iter",
                                       total = numChunks, clear = FALSE)
        # pb <- utils::txtProgressBar(0, numChunks, style = 3)
    for (i in seq(numChunks)) {
        start <- (i - 1)*chunkSize + 1
        end <- if (i*chunkSize > ncol(object)) ncol(object) else i*chunkSize
        colptrStart <- start
        colptrEnd <- end + 1
        chunkColptr <- colptr[colptrStart:colptrEnd]
        nnzStart <- chunkColptr[1] + 1
        nnzEnd <- chunkColptr[length(chunkColptr)]
        chunkData <- data[nnzStart:nnzEnd] # This step is freaking slow
        chunkRowind <- rowind[nnzStart:nnzEnd] # This step is freaking slow

        chunkColptr <- chunkColptr - chunkColptr[1]

        chunk <- Matrix::sparseMatrix(i = chunkRowind + 1, p = chunkColptr,
                                      x = chunkData,
                                      dims = c(numFeatures, end - start + 1),
                                      dimnames = list(rownames(object),
                                                      colnames(object)[start:end]))
        init <- do.call(FUN, c(list(chunk, nnzStart:nnzEnd,
                                    start:end, init),
                               fun.args))
        # if (isTRUE(verbose)) utils::setTxtProgressBar(pb, i)
        if (isTRUE(verbose)) cli::cli_progress_update(id = cliID, set = i)
    }
    # Break a new line otherwise next message comes right after the "%" sign.
    if (isTRUE(verbose)) cat("\n")
    init
}

# Safely add new H5 Data to the HDF5 file in a ligerDataset object
safeH5Create <- function(object,
                         dataPath,
                         dims,
                         dtype = "double",
                         chunkSize = dims) {
    if (inherits(object, "ligerDataset") && isH5Liger(object)) {
        h5file <- getH5File(object)
    } else if (inherits(object, "H5File")) {
        h5file <- object
    }

    # Check/Create H5Group ####
    # Inspect given `dataPath` b/c `hdf5r` does not allow creating dataset w/
    # "group" path(s)
    dataPath <- .checkArgLen(dataPath, n = 1, class = "character")
    dataPath <- trimws(dataPath, whitespace = "/")
    dataPath <- strsplit(dataPath, "/")[[1]]
    if (length(dataPath) > 1) {
        dataPath <- list(groups = dataPath[seq(length(dataPath) - 1)],
                         data = dataPath[length(dataPath)])
    } else {
        dataPath <- list(groups = NULL, data = dataPath)
    }
    if (!is.null(dataPath$groups)) {
        for (depth in seq_along(dataPath$groups)) {
            pathNow <- paste(dataPath$groups[seq(depth)], collapse = "/")
            if (!h5file$exists(pathNow)) {
                h5file$create_group(pathNow)
            }
        }
    }
    dataPath <- paste0(paste(dataPath$groups, collapse = "/"),
                       "/", dataPath$data)

    # Now we can work on the H5D data itself ####
    if (!h5file$exists(dataPath)) {
        # If specified data does not exist yet, just simply create the link
        h5file$create_dataset(
            name = dataPath,
            dims = dims,
            dtype = hdf5r::h5types[[dtype]],
            chunk_dims = chunkSize
        )
    } else {
        # If it already exist, need to check the dimension and extend as needed.
        originalDims <- h5file[[dataPath]]$dims
        if (length(dims) == 1) {
            if (originalDims < dims) {
                hdf5r::extendDataSet(h5file[[dataPath]], dims)
            } else if (originalDims > dims) {
                h5file$link_delete(dataPath)
                h5file$create_dataset(
                    name = dataPath,
                    dims = dims,
                    dtype = hdf5r::h5types[[dtype]],
                    chunk_dims = chunkSize
                )
            }
        } else if (length(dims) == 2) {
            # Only check for the 1st dim for now (number of feature)
            if (originalDims[1] < dims[1]) {
                hdf5r::extendDataSet(h5file[[dataPath]], dims)
            } else if (originalDims[1] > dims[1]) {
                h5file$link_delete(dataPath)
                h5file$create_dataset(
                    name = dataPath,
                    dims = dims,
                    dtype = hdf5r::h5types[[dtype]],
                    chunk_dims = chunkSize
                )
            }
        }
    }
}

#' Restore links (to HDF5 files) for reloaded liger/ligerDataset object
#' @description When loading the saved liger object with HDF5 data in a new R
#' session, the links to HDF5 files would be closed. This function enables
#' the restoration of those links so that new analyses can be carried out.
#' @param object \linkS4class{liger} or \linkS4class{ligerDataset} object.
#' @param filePath Paths to HDF5 files. A single character path for
#' \linkS4class{ligerDataset} input or a list of paths named by the datasets for
#' \linkS4class{liger} object input. Default \code{NULL} looks for the path(s)
#' of the last valid loading.
#' @return \code{object} with restored links.
#' @rdname restoreH5Liger
#' @export
#' @examples
#' h5Path <- system.file("extdata/ctrl.h5", package = "rliger")
#' tempPath <- tempfile(fileext = ".h5")
#' file.copy(from = h5Path, to = tempPath)
#' lig <- createLiger(list(ctrl = tempPath))
#' # Now it is actually an invalid object! which is equivalent to what users
#' # will get with `saveRDS(lig, "object.rds"); lig <- readRDS("object.rds")``
#' closeAllH5(lig)
#' lig <- restoreH5Liger(lig)
restoreH5Liger <- function(object, filePath = NULL) {
    if (!inherits(object, "liger") && !inherits(object, "ligerDataset")) {
        cli::cli_abort("Please specify a {.cls liger} or {.cls ligerDataset} object to restore.")
    }
    if (inherits(object, "ligerDataset")) {
        if (isTRUE(methods::validObject(object, test = TRUE))) {
            return(object)
        }
        h5.meta <- h5fileInfo(object)
        if (is.null(filePath)) filePath <- h5.meta$filename
        if (is.null(filePath)) {
            cli::cli_abort("No filename identified.")
        }
        if (!file.exists(filePath)) {
            cli::cli_abort("HDF5 file path does not exist: {.file {filePath}}")
        }
        cliID <- cli::cli_process_start("Restoring HDF5 link from: {.file {filePath}}")
        h5file <- hdf5r::H5File$new(filePath, mode = "r+")
        h5.meta$filename <- h5file$filename
        pathChecks <- unlist(lapply(h5.meta[4:10], function(x) {
            if (!is.null(x)) h5file$link_exists(x)
            else TRUE
        }))
        if (any(!pathChecks)) {
            info.name <- names(pathChecks)[!pathChecks]
            paths <- unlist(h5.meta[info.name])
            errMsg_cli <- paste0("HDF5 info {.val ", info.name, "} not found at path: {.val ", paths, "}")
            lapply(errMsg_cli, cli::cli_alert_danger)
            cli::cli_abort(
                "Cannot restore this dataset."
            )
            # errorMsg <- paste(paste0('HDF5 info "', info.name,
            #                          '" not found at path: "', paths, '"'),
            #                   collapse = "\n  ")
            # stop(errorMsg)
        }
        barcodes <- h5file[[h5.meta$barcodesName]]
        if (identical(barcodes, colnames(object))) {
            cli::cli_abort("Barcodes in the HDF5 file do not match to object.")
        }
        features <- h5file[[h5.meta$genesName]]
        if (identical(features, rownames(object))) {
            cli::cli_abort("Features in the HDF5 file do not match to object.")
        }
        # All checks passed!
        h5.meta$H5File <- h5file
        h5fileInfo(object, check = FALSE) <- h5.meta
        rawData(object, check = FALSE) <- h5file[[h5.meta$rawData]]
        if (!is.null(h5.meta$normData))
            normData(object, check = FALSE) <- h5file[[h5.meta$normData]]
        if (!is.null(h5.meta$scaleData)) {
            scaleData(object, check = FALSE) <- h5file[[h5.meta$scaleData]]
        }
        methods::validObject(object)
        cli::cli_process_done(id = cliID)
    } else {
        # Working for liger object
        if (!is.null(filePath)) {
            if (!is.list(filePath) || is.null(names(filePath)))
                cli::cli_abort(
                    "{.var filePath} has to be named list of {.cls liger} objects."
                )
        }
        for (d in names(object)) {
            if (isH5Liger(object, d)) {
                path <- NULL
                if (d %in% names(filePath)) {
                    if (!hdf5r::is.h5file(filePath[[d]])) {
                        cli::cli_alert_danger("Path for dataset {.val {d}} is not an HDF5 file: {.file {filePath[[d]]}}")
                    } else path <- filePath[[d]]
                }
                cliID <- cli::cli_process_start("Restoring dataset {.val {d}}")
                datasets(object, check = FALSE)[[d]] <-
                    restoreH5Liger(dataset(object, d), filePath[[d]])
                cli::cli_process_done(id = cliID)
            }
        }
    }
    return(object)
}

#' @note
#' \code{restoreOnlineLiger} will be deprecated for clarifying the terms used
#' for data structure.
#' @rdname restoreH5Liger
#' @export
#' @param file.path Will be deprecated with \code{restoreOnlineLiger}. The same
#' as \code{filePath}.
restoreOnlineLiger <- function(object, file.path = NULL) {
    lifecycle::deprecate_warn("1.99.0", "restoreOnlineLiger()",
                              "restoreH5Liger(object, filePath)")
    restoreH5Liger(object, file.path)
}

.inspectH5Path <- function(path) {
    if (length(path) != 1 || !is.character(path)) {
        cli::cli_abort("{.var path} has to be a single {.cls character}.")
    }
    path <- trimws(path, whitespace = "/")
    path <- strsplit(path, "/")[[1]]
    if (length(path) > 1) {
        list(folder = path[seq(length(path) - 1)],
             data = path[length(path)])
    } else {
        list(folder = NULL, data = path)
    }
}

#' Close all links (to HDF5 files) of a liger object
#' @description When need to interact with the data embedded in HDF5 files out
#' of the currect R session, the HDF5 files has to be closed in order to be
#' available to other processes.
#' @param object liger object.
#' @return Nothing is returned.
#' @export
#' @rdname closeAllH5
closeAllH5 <- function(object) {
    UseMethod("closeAllH5", object)
}

#' @rdname closeAllH5
#' @export
#' @method closeAllH5 liger
closeAllH5.liger <- function(object) {
    for (dn in names(object)) {
        ld <- dataset(object, dn)
        closeAllH5(ld)
    }
}

#' @rdname closeAllH5
#' @export
#' @method closeAllH5 ligerDataset
closeAllH5.ligerDataset <- function(object) {
    h5file <- getH5File(object)
    if (!is.null(h5file)) {
        path <- h5file$filename
        cli::cli_alert_info("Closing H5 file: {.file {path}}")
        h5file$close_all()
    }
    return(invisible(NULL))
}

.H5GroupToH5SpMat <- function(obj, dims) {
    groupPath <- obj$get_obj_name()
    RcppPlanc::H5SpMat(filename = obj$get_filename(),
                       valuePath = paste0(groupPath, "/data"),
                       rowindPath = paste0(groupPath, "/indices"),
                       colptrPath = paste0(groupPath, "/indptr"),
                       nrow = dims[1], ncol = dims[2])
}

# .H5DToH5Mat <- function(obj) {
#     RcppPlanc::H5Mat(filename = obj$get_filename(),
#                      dataPath = obj$get_obj_name())
# }





#' Write in-memory data into H5 file
#' @rdname writeH5
#' @description
#' This function writes in-memory data into H5 file by default in 10x cellranger
#' HDF5 output format. The main goal of this function is to allow users to
#' integrate large H5-based dataset, that cannot be fully loaded into memory,
#' with other data already loaded in memory using \code{\link{runOnlineINMF}}.
#' In this case, users can write the smaller in-memory data to H5 file instead
#' of loading subset of the large H5-based dataset into memory, where
#' information might be lost.
#'
#' Basing on the goal of the whole workflow, the data will always be written
#' in a CSC matrix format and colnames/rownames are always required.
#'
#' The default method coerces the input to a \linkS4class{dgCMatrix}. Methods
#' for other container classes tries to extract proper data and calls the
#' default method.
#' @param x An object with in-memory data to be written into H5 file.
#' @param file A character string of the file path to be written.
#' @param overwrite Logical, whether to overwrite the file if it already exists.
#' Default \code{FALSE}.
#' @param indicesPath,indptrPath,dataPath The paths inside the H5 file where
#' the \linkS4class{dgCMatrix} constructor \code{i}, \code{p}, and \code{x} will
#' be written to, respectively. Default using cellranger convention
#' \code{"matrix/indices"}, \code{"matrix/indptr"}, and \code{"matrix/data"}.
#' @param shapePath The path inside the H5 file where the shape of the matrix
#' will be written to. Default \code{"matrix/shape"}.
#' @param barcodesPath The path inside the H5 file where the barcodes/colnames
#' will be written to. Default \code{"matrix/barcodes"}. Skipped if the object
#' does not have colnames.
#' @param featuresPath The path inside the H5 file where the features/rownames
#' will be written to. Default \code{"matrix/features/name"}. Skipped if the
#' object does not have rownames.
#' @param useDatasets For liger method. Names or indices of datasets to be
#' written to H5 files. Required.
#' @param ... Arguments passed to other S3 methods.
#' @export
#' @seealso
#' \href{https://www.10xgenomics.com/cn/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices}{10X cellranger H5 matrix detail}
#' @return Nothing is returned. H5 file will be created on disk.
#' @examples
#' raw <- rawData(pbmc, "ctrl")
#' writeH5(raw, tempfile(pattern = "ctrl_", fileext = ".h5"))
writeH5 <- function(x, file, ...) {
    UseMethod("writeH5", x)
}

#' @rdname writeH5
#' @export
#' @method writeH5 default
writeH5.default <- function(x, file, ...) {
    x <- methods::as(x, "CsparseMatrix")
    writeH5(x, file, ...)
}

#' @rdname writeH5
#' @export
#' @method writeH5 dgCMatrix
writeH5.dgCMatrix <- function(
        x,
        file,
        overwrite = FALSE,
        indicesPath = "matrix/indices",
        indptrPath = "matrix/indptr",
        dataPath = "matrix/data",
        shapePath = "matrix/shape",
        barcodesPath = "matrix/barcodes",
        featuresPath = "matrix/features/name",
        ...
) {
    if (file.exists(file)) {
        if (isFALSE(overwrite)) {
            cli::cli_abort(
                c("x" = "File already exists at: {.file {normalizePath(file)}}.",
                  "i" = "Use {.code overwrite = TRUE} to overwrite.")
            )
        } else {
            file.remove(file)
        }
    }

    if (any(sapply(dimnames(x), is.null))) {
        cli::cli_abort("Both rownames and colnames are required.")
    }

    dataType <- typeof(x@x)
    dataDtype <- switch(
        dataType,
        double = "double",
        integer = "uint32_t",
        cli::cli_abort(
            "Unsupported data type in the sparse matrix: {.val {dataType}}"
        )
    )
    h5file <- hdf5r::H5File$new(file, mode = "w")
    safeH5Create(
        object = h5file,
        dataPath = indicesPath,
        dims = length(x@i),
        dtype = "uint32_t",
        chunkSize = 4096
    )
    safeH5Create(
        object = h5file,
        dataPath = indptrPath,
        dims = length(x@p),
        dtype = "uint32_t",
        chunkSize = 2048
    )
    safeH5Create(
        object = h5file,
        dataPath = dataPath,
        dims = length(x@x),
        dtype = dataDtype,
        chunkSize = 4096
    )
    safeH5Create(
        object = h5file,
        dataPath = shapePath,
        dims = 2,
        dtype = "uint64_t"
    )
    h5file[[indicesPath]][] <- x@i
    h5file[[indptrPath]][] <- x@p
    h5file[[dataPath]][] <- x@x
    h5file[[shapePath]][] <- dim(x)

    safeH5Create(
        object = h5file,
        dataPath = barcodesPath,
        dims = ncol(x),
        dtype = "char"
    )
    h5file[[barcodesPath]][] <- colnames(x)

    safeH5Create(
        object = h5file,
        dataPath = featuresPath,
        dims = nrow(x),
        dtype = "char"
    )
    h5file[[featuresPath]][] <- rownames(x)

    h5file$close()
    invisible(NULL)
}

#' @rdname writeH5
#' @export
#' @method writeH5 ligerDataset
writeH5.ligerDataset <- function(x, file, ...) {
    raw <- rawData(x)
    if (is.null(raw)) {
        cli::cli_abort("No {.code rawData(x)} available.")
    }
    writeH5(raw, file, ...)
}

#' @rdname writeH5
#' @export
#' @method writeH5 liger
writeH5.liger <- function(x, file, useDatasets, ...) {
    useDatasets <- .checkUseDatasets(x, useDatasets)
    file <- .checkArgLen(file, n = length(useDatasets), repN = FALSE,
                         class = "character", .stop = TRUE)
    for (i in seq_along(useDatasets)) {
        if (isH5Liger(x, i)) {
            cli::cli_alert_warning(
                "Dataset {.val {useDatasets[i]}} is H5 based, file located at: {.file {getH5File(x, useDatasets[i])$filename}}. Skipped."
            )
            next
        }
        raw <- rawData(x, useDatasets[i])
        if (is.null(raw)) {
            cli::cli_abort("No {.code rawData(x, '{useDatasets[i]}')} available.")
        }
        tryCatch(
            expr = {
                writeH5(raw, file[i], ...)
            },
            error = function(e) {
                cli::cli_alert_danger("Failed to write dataset {.val {useDatasets[i]}} to H5 file at {.file {file[i]}}. Continuing with the others.")
                msg <- format(e)
                cli::cli_alert_danger("Call {.code {msg[['call']]}}: {msg[['message']]}")
            }
        )
    }
    invisible(NULL)
}
