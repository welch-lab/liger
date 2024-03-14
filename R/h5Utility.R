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
        Sys.sleep(0.1)
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
#' lig <- createLiger(list(ctrl = h5Path))
#' # Now it is actually an invalid object! which is equivalent to what users
#' # will get with `saveRDS(lig, "object.rds"); lig <- readRDS("object.rds")``
#' lig <- closeAllH5(lig)
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
#' @return \code{object} with links closed.
#' @export
closeAllH5 <- function(object) {
    if (!isH5Liger(object)) return(object)
    for (dn in names(object)) {
        if (!isH5Liger(object, dn)) next
        h5f <- getH5File(object, dn)
        h5f$close_all()
    }
    return(object)
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
