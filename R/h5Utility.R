#ctrl.h5 <- hdf5r::H5File$new("pbmcs_ctrl.h5")

#' Apply function to chunks of H5 data
#' @description h5 calculation wrapper, that run user specified calculation with on-disk
#' matrix in chunks
#' @param h5file H5File object
#' @param FUN A function where the first argument must be the chunk of data, the
#' return value must be a vector.
H5Apply <- function(
        h5file,
        init,
        FUN,
        dataPath = "matrix/data",
        indicesPath = "matrix/indices",
        indptrPath = "matrix/indptr",
        cellIDPath = "matrix/barcodes",
        featureIDPath = "matrix/features/name",
        chunkSize = 1000,
        #vectorType = c("integer", "numeric", "character"),
        verbose = TRUE,
        ...) {
    fun.args <- list(...)
    #vectorType <- match.arg(vectorType)

    numCells <- h5file[[cellIDPath]]$dims
    numFeatures <- h5file[[featureIDPath]]$dims
    #initVec <- do.call(vectorType, list(numCells))
    #initVec <- c()

    prev_end_col <- 1
    prev_end_data <- 1
    numChunks <- ceiling(numCells / chunkSize)
    ind = 0

    if (isTRUE(verbose)) pb <- utils::txtProgressBar(0, numChunks, style = 3)

    while (prev_end_col < numCells) {
        ind <- ind + 1
        # Calculate the proper position/index for extract the data of the chunk
        if (numCells - prev_end_col < chunkSize) {
            chunkSize <- numCells - prev_end_col + 1
        }
        start_inds <-
            h5file[[indptrPath]][prev_end_col:(prev_end_col + chunkSize)]
        sparseXIdx <- (prev_end_data):(tail(start_inds, 1))
        row_inds <- h5file[[indicesPath]][sparseXIdx]
        counts <- h5file[[dataPath]][sparseXIdx]
        # Construct sparse matrix of the chunk
        chunkData <- Matrix::sparseMatrix(
            i = row_inds[seq_along(counts)] + 1,
            p = start_inds[seq(chunkSize + 1)] - prev_end_data + 1,
            x = counts,
            dims = c(numFeatures, chunkSize)
        )
        # Process chunk data with given function and additional arguments if
        # applicable. Then insert value to initialized output data.
        #chunkResult <- FUN(chunkData, ...)
        init <- FUN(chunkData, sparseXIdx, init, ...)
        #if (!is.vector(chunkResult) | length(chunkResult) != chunkSize) {
        #    stop("Chunk ", ind, " did not return proper vector result.")
        #}
        #cellIdx <- seq(prev_end_col, prev_end_col + chunkSize - 1)
        #initVec[cellIdx] <- chunkResult
        #initVec <- c(initVec, chunkResult)
        ############################
        # Post-processing updates on indices for the next iteration
        prev_end_col <- prev_end_col + chunkSize
        prev_end_data <- prev_end_data + length(chunkData@x)

        if (isTRUE(verbose)) utils::setTxtProgressBar(pb, ind)
    }
    init
}

#' Safely add new H5 Data to the HDF5 file in a ligerDataset object
#' @noRd
safeH5Create <- function(object,
                        dataPath,
                        dims,
                        dtype = "double",
                        chunkSize = dims) {
    h5file <- getH5File(object)
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
            }
        } else if (length(dims) == 2) {
            # Only check for the 1st dim for now (number of feature)
            if (originalDims[1] < dims[1]) {
                hdf5r::extendDataSet(h5file[[dataPath]], dims)
            }
        }
    }
}

#' Restore links (to HDF5 files) for reloaded online liger object
#' @description When loading the saved online liger object in a new R session,
#' the links to HDF5 files may be corrupted. This functions enables the
#' restoration of those links so that new analyses can be carried out.
#' @param object liger or ligerDataset object.
#' @param file.path List of paths to hdf5 files.
#' @return \code{object} with restored links.
#' @export
restoreOnlineLiger <- function(object, filePath = NULL) {
    if (!inherits(object, "liger") & !inherits(object, "ligerDataset")) {
        stop("Please specify a liger or ligerDataset object to restore.")
    }
    if (inherits(object, "ligerDataset")) {
        if (isTRUE(validObject(object, test = TRUE))) {
            return(object)
        }
        h5.meta <- object@h5file.info
        if (is.null(filePath)) filePath <- h5.meta$filename
        if (!file.exists(filePath)) {
            stop("HDF5 file path does not exist:\n",
                 filePath)
        }
        message(date(), " ... Restoring ligerDataset from: ", filePath)
        h5file <- hdf5r::H5File$new(filePath, mode = "r+")
        #h5.meta$H5File <- h5file
        h5.meta$filename <- h5file$filename
        pathChecks <- unlist(lapply(h5.meta[4:10], h5file$link_exists))
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
        object@h5file.info <- h5.meta
        raw.data(object, check = FALSE) <- h5file[[h5.meta$raw.data]]
        norm.data(object, check = FALSE) <- h5file[[h5.meta$norm.data]]
        scale.data(object, check = FALSE) <- h5file[[h5.meta$scale.data]]
        validObject(object)
    }
    return(object)
}


normH5 <- function(h5file) {
    num_genes <- h5file[["matrix/features/name"]]$dims
    num_cells <- h5file[["matrix/barcodes"]]$dims

    # Initialize result
    results <- list(
        nUMI = c(),
        geneSumSq = rep(0, num_genes),
        geneMeans = rep(0, num_genes)
    )
    # Use safe create here in practice
    #h5file$create_dataset(
    #    name = "norm_data2",
    #    dims = h5file[["matrix/data"]]$dims,
    #    dtype = hdf5r::h5types$double,
    #    chunk_dims = 1000
    #)
    # Chunk run
    results <- H5Apply(h5file, init = results,
                       function(chunk, sparseXIdx, values) {
        values$nUMI <- c(values$nUMI, colSums(chunk))
        norm.data <- Matrix.column_norm(chunk)
        h5file[["norm_data2"]][sparseXIdx] <- norm.data@x
        row_sums <- rowSums(norm.data)
        values$geneSumSq <- values$geneSumSq + rowSums(norm.data * norm.data)
        values$geneMeans <<- values$geneMeans + row_sums
        values
    })
    results$geneMeans <- results$geneMeans / num_cells
    return(results)
}

#res1 <- normH5(ctrl.h5)

#' @useDynLib rliger, .registration = TRUE
calcGeneVar.h5 <- function(h5file) {
    num_genes <- h5file[["matrix/features/name"]]$dims
    num_cells <- h5file[["matrix/barcodes"]]$dims

    gene_vars <- rep(0, num_genes)
    gene_means <- res1$geneMeans
    print(gene_means[1:10])
    gene_vars <- H5Apply(h5file, gene_vars,
                         function(chunk, sparseXIdx, values) {
        values + sumSquaredDeviations(chunk, gene_means)
    }, dataPath = "norm.data")
    gene_vars / (num_cells - 1)
}

#gene_vars <- calcGeneVar.h5(ctrl.h5)
# EXAMPLE ##########
#library(Matrix)
# Read data
#ctrl.h5 <- hdf5r::H5File$new("pbmcs_ctrl.h5")
#num_genes <- ctrl.h5[["matrix/features/name"]]$dims
#num_cells <- ctrl.h5[["matrix/barcodes"]]$dims

# Initialize result
#nUMI.test <- c()
#gene_sum_sq = rep(0, num_genes)
#gene_means = rep(0, num_genes)
#ctrl.h5$create_dataset(
#    name = "norm_data2",
#    dims = ctrl.h5[["matrix/data"]]$dims,
#    dtype = hdf5r::h5types$double,
#    chunk_dims = 1000
#)

# Chunk run
#H5Apply(ctrl.h5, function(chunk, sparseXIdx) {
#    nUMI.test <<- c(nUMI.test, colSums(chunk))
#    norm.data <- Matrix.column_norm(chunk)
#    ctrl.h5[["norm_data2"]][sparseXIdx] <<- norm.data@x
#    row_sums <- rowSums(norm.data)
#    gene_sum_sq <<- gene_sum_sq + rowSums(norm.data * norm.data)
#    gene_means <<- gene_means + row_sums
#})
#gene_means <- gene_means / num_cells

# Check
#all(ctrl.h5[["norm.data"]][] == ctrl.h5[["norm_data2"]][])

#ctrl.dgC <- readRDS("pbmcs_ctrl.RDS")
#nUMI.true <- colSums(ctrl.dgC)
#all(nUMI.test == nUMI.true)

#ctrl.dgC.norm <- Matrix.column_norm(ctrl.dgC)
#gene_means.true <- rowMeans(ctrl.dgC.norm)
#all(gene_means.true == gene_means)
#all(round(gene_means.true, 10) == round(gene_means, 10))

#gene_sum_sq.true <- rowSums(ctrl.dgC.norm * ctrl.dgC.norm)
#all(round(gene_sum_sq, 10) == round(gene_sum_sq.true, 10))
