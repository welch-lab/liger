#ctrl.h5 <- hdf5r::H5File$new("pbmcs_ctrl.h5")

#' Apply function to chunks of H5 data in ligerDataset object
#' @description h5 calculation wrapper, that run user specified calculation with
#' on-disk matrix in chunks
#' @param h5file H5File object
#' @param FUN A function where the first argument must be the chunk of data, the
#' return value must be a vector.
H5Apply <- function(
        object,
        FUN,
        init = NULL,
        useData = c("raw.data", "norm.data", "scale.data"),
        chunkSize = 1000,
        verbose = TRUE,
        ...) {
    fun.args <- list(...)
    useData <- match.arg(useData)
    h5meta <- h5file.info(object)
    numCells <- ncol(object)
    numFeatures <- nrow(object)

    prev_end_col <- 1
    prev_end_data <- 1
    numChunks <- ceiling(numCells / chunkSize)
    ind <- 0
    h5file <- h5meta$H5File
    indptr <- h5file[[h5meta$indptr.name]]
    indices <- h5file[[h5meta$indices.name]]
    data <- h5file[[h5meta[[useData]]]]
    if (isTRUE(verbose)) pb <- utils::txtProgressBar(0, numChunks, style = 3)

    while (prev_end_col < numCells) {
        ind <- ind + 1
        # Calculate the proper position/index for extract the data of the chunk
        if (numCells - prev_end_col < chunkSize) {
            chunkSize <- numCells - prev_end_col + 1
        }
        start_inds <- indptr[prev_end_col:(prev_end_col + chunkSize)]
        sparseXIdx <- (prev_end_data):(tail(start_inds, 1))
        cellIdx <- prev_end_col:(prev_end_col + chunkSize - 1)
        row_inds <- indices[sparseXIdx]
        counts <- data[sparseXIdx]
        # Construct sparse matrix of the chunk
        chunkData <- Matrix::sparseMatrix(
            i = row_inds[seq_along(counts)] + 1,
            p = start_inds[seq(chunkSize + 1)] - prev_end_data + 1,
            x = counts,
            dims = c(numFeatures, chunkSize)
        )
        # Process chunk data with given function and additional arguments if
        # applicable. Then insert value to initialized output data.
        init <- do.call(FUN, c(list(chunkData,
                                    sparseXIdx,
                                    cellIdx,
                                    init),
                               fun.args))
        ############################
        # Post-processing updates on indices for the next iteration
        prev_end_col <- prev_end_col + chunkSize
        prev_end_data <- prev_end_data + length(chunkData@x)

        if (isTRUE(verbose)) utils::setTxtProgressBar(pb, ind)
    }
    # Break a new line otherwise next message comes right after the "%" sign.
    if (isTRUE(verbose)) cat("\n")
    init
}

#' Safely add new H5 Data to the HDF5 file in a ligerDataset object
#' @noRd
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
    if (length(dataPath) != 1 | !is.character(dataPath)) {
        stop("`path` has to be a single character.")
    }
    dataPath <- trimws(dataPath, whitespace = "/")
    dataPath <- strsplit(dataPath, "/")[[1]]
    if (length(dataPath) > 1) {
        dataPath <- list(groups = dataPath[seq(length(dataPath) - 1)],
                         data = dataPath[length(dataPath)])
    } else {
        dataPath <- list(groups = NULL, data = dataPath)
    }
    if (!is.null(dataPath$groups)) {
        for (depth in seq(length(dataPath$groups))) {
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
#' session, the links to HDF5 files would be corrupted. This functions enables
#' the restoration of those links so that new analyses can be carried out.
#' @param object liger or ligerDataset object.
#' @param file.path Paths to HDF5 files. A single character path for
#' ligerDataset input or a list of paths named by the datasets for liger object
#' input. Default \code{NULL} looks for the path(s) of the last valid loading.
#' @return \code{object} with restored links.
#' @export
restoreH5Liger <- function(object, filePath = NULL) {
    if (!inherits(object, "liger") & !inherits(object, "ligerDataset")) {
        stop("Please specify a liger or ligerDataset object to restore.")
    }
    if (inherits(object, "ligerDataset")) {
        if (isTRUE(validObject(object, test = TRUE))) {
            return(object)
        }
        h5.meta <- h5file.info(object)
        if (is.null(filePath)) filePath <- h5.meta$filename
        if (!file.exists(filePath)) {
            stop("HDF5 file path does not exist:\n",
                 filePath)
        }
        message(date(), " ... Restoring ligerDataset from: ", filePath)
        h5file <- hdf5r::H5File$new(filePath, mode = "r+")
        h5.meta$filename <- h5file$filename
        pathChecks <- unlist(lapply(h5.meta[4:10], function(x) {
            if (!is.null(x)) h5file$link_exists(x)
            else TRUE
        }))
        if (any(!pathChecks)) {
            info.name <- names(pathChecks)[!pathChecks]
            paths <- unlist(h5.meta[info.name])
            errorMsg <- paste(paste0('HDF5 info "', info.name,
                                     '" not found at path: "', paths, '"'),
                              collapse = "\n  ")
            stop(errorMsg)
        }
        barcodes <- h5file[[h5.meta$barcodes.name]]
        if (identical(barcodes, colnames(object))) {
            stop("Barcodes in the HDF5 file do not match to object.")
        }
        features <- h5file[[h5.meta$genes.name]]
        if (identical(features, rownames(object))) {
            stop("Features in the HDF5 file do not match to object.")
        }
        # All checks passed!
        h5.meta$H5File <- h5file
        h5file.info(object, check = FALSE) <- h5.meta
        raw.data(object, check = FALSE) <- h5file[[h5.meta$raw.data]]
        if (!is.null(h5.meta$norm.data))
            norm.data(object, check = FALSE) <- h5file[[h5.meta$norm.data]]
        if (!is.null(h5.meta$scale.data))
            scale.data(object, check = FALSE) <- h5file[[h5.meta$scale.data]]
        validObject(object)
    } else {
        # Working for liger object
        if (!is.null(filePath)) {
            if (!is.list(filePath) || is.null(names(filePath)))
                stop("`filePath` has to be a named list for liger object.")
            if (!any(names(filePath) %in% names(object)))
                stop("names of `filePath` must be found in `names(object)`.")
        }
        for (d in names(object)) {
            if (isH5Liger(object, d)) {
                path <- NULL
                if (d %in% names(filePath)) {
                    if (!hdf5r::is.h5file(filePath[[d]]))
                        warning("Path for dataset \"", d, "\" is not an HDF5 file: ",
                                filePath[[d]])
                    else path <- filePath[[d]]
                }
                dataset(object, d, qc = FALSE) <-
                    restoreH5Liger(dataset(object, d), filePath[[d]])
            }
        }
    }
    return(object)
}

.inspectH5Path <- function(path) {
    if (length(path) != 1 | !is.character(path)) {
        stop("`path` has to be a single character.")
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

