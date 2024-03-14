.cbind.ligerDataset.mem <- function(args) {
    # TODO: Now by default merging the three slots
    # NEED TO add checks on existence and think about adding NA if any slot is
    # missing for one dataset
    libraryNames <- NULL
    if (!is.null(names(args))) libraryNames <- names(args)
    if (any(is.na(libraryNames)))
        libraryNames[is.na(libraryNames)] <- paste0("Dataset_",
                                                    which(is.na(libraryNames)))
    newColnames <- unlist(lapply(args, colnames), use.names = FALSE)
    if (!is.null(libraryNames))
        newColnames <- paste0(rep(libraryNames, times = lapply(args, ncol)),
                              "_", newColnames)
    if (all(sapply(args, function(ld) !is.null(rawData(ld)))))
        rawData <- mergeSparseAll(lapply(args, rawData),
                                  libraryNames = libraryNames)
    else rawData <- NULL
    if (all(sapply(args, function(ld) !is.null(normData(ld)))))
        normData <- mergeSparseAll(lapply(args, normData),
                                   libraryNames = libraryNames)
    else normData <- NULL
    # if (all(sapply(args, function(ld) !is.null(scaleData(ld)))))
    #     scaleData <- mergeSparseAll(lapply(args, scaleData),
    #                                 libraryNames = libraryNames)
    # else scaleData <- NULL
    createLigerDataset(rawData = rawData, normData = normData)#,
                       # scaleData = scaleData)
}

.cbind.ligerDataset.h5 <- function(args) stop("Not implemented yet")

#' Merge matrices while keeping the union of rows
#' @rdname mergeSparseAll
#' @description \code{mergeSparseAll} takes in a list of DGEs, with genes as
#' rows and cells as columns, and merges them into a single DGE. Also adds
#' \code{libraryNames} to colnames from each DGE if expected to be overlap
#' (common with 10X barcodes). Values in \code{rawData} or \code{normData}
#' slot of a \linkS4class{ligerDataset} object can be processed with this.
#'
#' For a list of dense matrices, usually the values in \code{scaleData} slot of
#' a \linkS4class{ligerDataset} object, please use \code{mergeDenseAll} which
#' works in the same way.
#' @param datalist List of dgCMatrix for \code{mergeSparseAll} or a list of
#' matrix for \code{mergeDenseAll}.
#' @param libraryNames Character vector to be added as the prefix for the
#' barcodes in each matrix in \code{datalist}. Length should match with the
#' number of matrices. Default \code{NULL} do not modify the barcodes.
#' @param mode Whether to take the \code{"union"} or \code{"intersection"} of
#' features when merging. Default \code{"union"}.
#' @return dgCMatrix or matrix with all barcodes in \code{datalist} as columns
#' and the union of genes in \code{datalist} as rows.
#' @export
#' @examples
#' rawDataList <- getMatrix(pbmc, "rawData")
#' merged <- mergeSparseAll(rawDataList, libraryNames = names(pbmc))
mergeSparseAll <- function(datalist, libraryNames = NULL,
                           mode = c("union", "intersection")) {
    # Use summary to convert the sparse matrices into three-column indexes where
    # i are the row numbers, j are the column numbers, and x are the nonzero
    # entries
    mode <- match.arg(mode)
    col_offset <- 0
    if (mode == "union") {
        allGenes <- Reduce(union, lapply(datalist, rownames))
    } else if (mode == "intersection") {
        allGenes <- Reduce(intersect, lapply(datalist, rownames))
    }
    allCells <- c()
    for (i in 1:length(datalist)) {
        curr <- datalist[[i]]
        curr_s <- summary(curr)
        # Now, alter the indexes so that the two 3-column matrices can be
        # properly merged.
        # First, make the current and full column numbers non-overlapping.
        curr_s[, 2] <- curr_s[, 2] + col_offset

        # Update full cell names
        if (!is.null(libraryNames)) {
            cellnames <- paste0(libraryNames[i], "_", colnames(curr))
        } else {
            cellnames <- colnames(curr)
        }
        allCells <- c(allCells, cellnames)

        # Next, change the row (gene) indexes so that they index on the union of the
        # gene sets, so that proper merging can occur.
        idx <- match(rownames(curr), allGenes)
        newgenescurr <- idx[curr_s[, 1]]
        curr_s[, 1] <- newgenescurr
        #curr_s <- curr_s[!is.na(curr_s[,1]),]
        # Now bind the altered 3-column matrices together, and convert into a single sparse matrix.
        if (!exists("full_mat")) {
            full_mat <- curr_s
        } else {
            full_mat <- rbind(full_mat, curr_s)
        }
        col_offset <- length(allCells)
    }
    full_mat <- full_mat[!is.na(full_mat[,1]),]
    M <- Matrix::sparseMatrix(
        i = full_mat[, 1],
        j = full_mat[, 2],
        x = full_mat[, 3],
        dims = c(length(allGenes),
                 length(allCells)),
        dimnames = list(allGenes,
                        allCells)
    )
    return(M)
}

#' @rdname mergeSparseAll
#' @export
mergeDenseAll <- function(datalist, libraryNames = NULL) {
    if (length(datalist) == 0) return(NULL)
    if (length(datalist) == 1) return(datalist[[1]])
    allRows <- unique(unlist(lapply(datalist, rownames), use.names = FALSE))
    dfList <- lapply(seq_along(datalist), function(i) {
        dat <- datalist[[i]]
        if (!is.null(libraryNames))
            colnames(dat) <- paste0(libraryNames[i], "_", colnames(dat))
        if (nrow(dat) == 0)
            dat <- matrix(NA, nrow = length(allRows), ncol = ncol(dat),
                          dimnames = list(allRows, colnames(dat)))
        dat <- as.data.frame(dat)
        rnCol <- data.frame(rn = rownames(dat))
        cbind(rnCol, dat)
    })

    # Could have use reduce. but `Reduce` doesn't allow additional arguments,
    # and I don't want to introduce extra dependency like "purrr"
    result <- dfList[[1]]
    for (i in seq(2, length(dfList))) {
        result <- dplyr::full_join(result, dfList[[i]], by = "rn")
    }
    rownames(result) <- result$rn
    result$rn <- NULL
    as.matrix(result)
}

#' Merge hdf5 files
#' @description This function merges hdf5 files generated from different
#' libraries (cell ranger by default) before they are preprocessed through Liger
#' pipeline.
#' @param file.list List of path to hdf5 files.
#' @param library.names Vector of library names (corresponding to file.list)
#' @param new.filename String of new hdf5 file name after merging (default
#' new.h5).
#' @param format.type string of HDF5 format (10X CellRanger by default).
#' @param data.name Path to the data values stored in HDF5 file.
#' @param indices.name Path to the indices of data points stored in HDF5 file.
#' @param indptr.name Path to the pointers stored in HDF5 file.
#' @param genes.name Path to the gene names stored in HDF5 file.
#' @param barcodes.name Path to the barcodes stored in HDF5 file.
#' @return Directly generates newly merged hdf5 file.
#' @export
#' @examples
#' \dontrun{
#' # For instance, we want to merge two datasets saved in HDF5 files (10X
#' # CellRanger) paths to datasets: "library1.h5","library2.h5"
#' # dataset names: "lib1", "lib2"
#' # name for output HDF5 file: "merged.h5"
#' mergeH5(list("library1.h5","library2.h5"), c("lib1","lib2"), "merged.h5")
#' }
mergeH5 <- function(file.list,
                    library.names,
                    new.filename,
                    format.type = "10X",
                    data.name = NULL,
                    indices.name = NULL,
                    indptr.name = NULL,
                    genes.name = NULL,
                    barcodes.name = NULL) {
    h5_merged = hdf5r::H5File$new(paste0(new.filename, ".h5"), mode = "w")
    h5_merged$create_group("matrix")
    h5_merged$create_group("matrix/features")
    num_data_prev = 0
    num_indptr_prev = 0
    num_cells_prev = 0
    last_inptr = 0
    for (i in 1:length(file.list)) {
        h5file = hdf5r::H5File$new(file.list[[i]], mode = "r")
        if (format.type == "10X") {
            data = h5file[["matrix/data"]][]
            indices = h5file[["matrix/indices"]][]
            indptr = h5file[["matrix/indptr"]][]
            barcodes = paste0(library.names[i], "_", h5file[["matrix/barcodes"]][])
            genes = h5file[["matrix/features/name"]][]
        } else if (format.type == "AnnData") {
            data = h5file[["raw.X/data"]][]
            indices = h5file[["raw.X/indices"]][]
            indptr = h5file[["raw.X/indptr"]][]
            barcodes = paste0(library.names[i], "_", h5file[["obs"]][]$cell)
            genes = h5file[["raw.var"]][]$index

        } else {
            data = h5file[[data.name]][]
            indices = h5file[[indices.name]][]
            indptr = h5file[[indptr.name]][]
            barcodes = paste0(library.names[i], "_", h5file[[barcodes.name]][])
            genes = h5file[[genes.name]][]
        }

        if (i != 1)
            indptr = indptr[2:length(indptr)]
        num_data = length(data)
        num_indptr = length(indptr)
        num_cells = length(barcodes)
        indptr = indptr + last_inptr
        last_inptr = indptr[num_indptr]
        if (i == 1) {
            h5_merged[["matrix/data"]] = data
            h5_merged[["matrix/indices"]] = indices
            h5_merged[["matrix/indptr"]] = indptr
            h5_merged[["matrix/barcodes"]] = barcodes
            h5_merged[["matrix/features/name"]] = genes
        } else {
            h5_merged[["matrix/data"]][(num_data_prev + 1):(num_data_prev + num_data)] = data
            h5_merged[["matrix/indices"]][(num_data_prev + 1):(num_data_prev + num_data)] = indices
            h5_merged[["matrix/indptr"]][(num_indptr_prev + 1):(num_indptr_prev + num_indptr)] = indptr
            h5_merged[["matrix/barcodes"]][(num_cells_prev + 1):(num_cells_prev + num_cells)] = barcodes
        }
        num_data_prev = num_data_prev + num_data
        num_indptr_prev = num_indptr_prev + num_indptr
        num_cells_prev = num_cells_prev + num_cells
        h5file$close_all()
    }
    h5_merged$close_all()
}
