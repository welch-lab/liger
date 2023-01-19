#' Downsample datasets
#' @description This function mainly aims at downsampling datasets to a size
#' suitable for plotting.
#' @param object \linkS4class{liger} object
#' @param balance \code{"all"} for sampling \code{maxCells} cells from all
#' datasets specified by \code{useDatasets}. \code{"cluster"} for sampling
#' \code{maxCells} cells per cluster per dataset. \code{"dataset"} for
#' \code{maxCells} cells per dataset.
#' @param maxCells Max number of cells to sample from the grouping based on
#' \code{balance}.
#' @param useDatasets Index selection of datasets to consider. Default
#' \code{NULL} for using all datasets.
#' @param seed Random seed for reproducibility. Default \code{1}.
#' @param ... Arguments passed to \code{\link{subsetLiger}}, where
#' \code{cellIdx} is occupied by internal implementation.
#' @return Subset of \linkS4class{liger} \code{object}.
#' @export
downsample <- function(
    object,
    balance = c("all", "cluster", "dataset"),
    maxCells = 1000,
    useDatasets = NULL,
    seed = 1,
    ...
) {
    # TODO: multi-variable category balancing
    balance <- match.arg(balance)
    set.seed(seed)
    useDatasets <- .checkUseDatasets(object, useDatasets)
    selected <- c()
    if (balance == "all") {
        useCells <- which(object$dataset %in% useDatasets)
        maxCells <- min(maxCells, length(useCells))
        selected <- sort(sample(useCells, maxCells))
    } else if (balance == "cluster") {
        if (!"cluster" %in% names(cell.meta(object))) {
            stop('"cluster" not found in `cell.meta(object)`')
        }
        if (!is.factor(object$cluster)) {
            object$cluster <- factor(object$cluster)
        }
        for (d in useDatasets) {
            for (c in levels(object$cluster)) {
                useCells <- which(object$dataset == d & object$cluster == c)
                maxCells <- min(maxCells, useCells)
                selected <- c(selected , sample(useCells, maxCells))
            }
        }
        selected <- sort(selected)
    } else if (balance == "dataset") {
        for (d in useDatasets) {
            ld <- dataset(object, d)
            maxCells <- min(maxCells, ncol(ld))
            selected.name <- sample(colnames(ld), maxCells)
            selected <- c(selected, selected.name)
        }
        # Requires that all barcodes in liger object is unique
        selected <- which(colnames(object) %in% selected)
    }
    subsetLiger(object = object, cellIdx = selected, ...)
}

#' [Deprecated] See \code{\link{downsample}}
#' @description This function mainly aims at downsampling datasets to a size
#' suitable for plotting.
#' @param object \linkS4class{liger} object
#' @param slot.use Only create subset from one or more of \code{"raw.data"},
#' \code{"norm.data"} and \code{"scale.data"}. Default \code{NULL} subsets the
#' whole \code{object} including downstream results.
#' @param balance \code{"all"} for sampling \code{maxCells} cells from all
#' datasets specified by \code{useDatasets}. \code{"cluster"} for sampling
#' \code{maxCells} cells per cluster per dataset. \code{"dataset"} for
#' \code{maxCells} cells per dataset.
#' @param max.cells Max number of cells to sample from the grouping based on
#' \code{balance}.
#' @param chunk Integer. Number of maximum number of cells in each chunk,
#' Default \code{1000}.
#' @param datasets.use Index selection of datasets to consider. Default
#' \code{NULL} for using all datasets.
#' @param genes.use Character vector. Subset features to this specified range.
#' Default \code{NULL} does not subset features.
#' @param rand.seed Random seed for reproducibility. Default \code{1}.
#' @param verbose Logical. Whether to show the progress. Default \code{TRUE}.
#' @return Subset of \linkS4class{liger} \code{object}.
#' @seealso \code{\link{downsample}}, \code{\link{subsetLiger}},
#' \code{\link{subsetLigerDataset}}
#' @export
readSubset <- function(
        object,
        slot.use = "norm.data",
        balance = NULL,
        max.cells = 1000,
        chunk = 1000,
        datasets.use = NULL,
        genes.use = NULL,
        rand.seed = 1,
        verbose = TRUE
) {
    .Deprecated("downsample")
    if (!is.null(balance)) balance <- match.arg(balance)
    else balance <- "all"
    downsample(object = object, balance = balance,
               maxCells = max.cells, useDatasets = datasets.use,
               useGenes = genes.use, useSlot = slot.use, seed = rand.seed,
               chunkSize = chunk, verbose = verbose)
}
