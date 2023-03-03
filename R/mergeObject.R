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
    if (all(sapply(args, function(ld) !is.null(raw.data(ld)))))
        rawData <- mergeSparseAll(lapply(args, raw.data),
                                  libraryNames = libraryNames)
    else rawData <- NULL
    if (all(sapply(args, function(ld) !is.null(norm.data(ld)))))
        normData <- mergeSparseAll(lapply(args, norm.data),
                                   libraryNames = libraryNames)
    else normData <- NULL
    if (all(sapply(args, function(ld) !is.null(scale.data(ld)))))
        scaleData <- mergeDenseAll(lapply(args, scale.data),
                                   libraryNames = libraryNames)
    else scaleData <- NULL
    createLigerDataset(raw.data = rawData, norm.data = normData,
                       scale.data = scaleData)
}

.cbind.ligerDataset.h5 <- function(args) stop("Not implemented yet")

#' Merge matrices while keeping the union of rows
#' @rdname mergeSparseAll
#' @description \code{mergeSparseAll} takes in a list of DGEs, with genes as
#' rows and cells as columns, and merges them into a single DGE. Also adds
#' \code{libraryNames} to colnames from each DGE if expected to be overlap
#' (common with 10X barcodes). Values in \code{raw.data} or \code{norm.data}
#' slot of a \linkS4class{ligerDataset} object can be processed with this.
#'
#' For a list of dense matrices, usually the values in \code{scale.data} slot of
#' a \linkS4class{ligerDataset} object, please use \code{mergeDenseAll} which
#' works in the same way.
#' @param datalist List of dgCMatrix for \code{mergeSparseAll} or a list of
#' matrix for \code{mergeDenseAll}.
#' @param libraryNames Character vector to be added as the prefix for the
#' barcodes in each matrix in \code{datalist}. Length should match with the
#' number of matrices. Default \code{NULL} do not modify the barcodes.
#' @return dgCMatrix or matrix with all barcodes in \code{datalist} as columns
#' and the union of genes in \code{datalist} as rows.
#' @export
mergeSparseAll <- function(datalist, libraryNames = NULL) {
    # Use summary to convert the sparse matrices into three-column indexes where
    # `i` are the row numbers, `j` are the column numbers, and x are the nonzero
    # entries
    col_offset <- 0
    allGenes <- unique(unlist(lapply(datalist, rownames)))
    allCells <- c()
    for (i in seq_along(datalist)) {
        curr <- datalist[[i]]
        curr_s <- summary(curr)

        # Now, alter the indexes so that the two 3-column matrices can be properly
        # merged.
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

        # Now bind the altered 3-column matrices together, and convert into a single sparse matrix.
        if (!exists("full_mat")) {
            full_mat <- curr_s
        } else {
            full_mat <- rbind(full_mat, curr_s)
        }
        col_offset <- length(allCells)
    }
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
        data.frame(rn = rownames(dat), dat)
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
