#' Downsample datasets
#' @description This function mainly aims at downsampling datasets to a size
#' suitable for plotting or expensive in-memmory calculation.
#'
#' Users can balance the sample size of categories of interests with
#' \code{balance}. Multi-variable specification to \code{balance} is supported,
#' so that at most \code{maxCells} cells will be sampled from each combination
#' of categories from the variables. For example, when two datasets are
#' presented and three clusters labeled across them, there would then be at most
#' \eqn{2 \times 3 \times maxCells} cells being selected. Note that
#' \code{"dataset"} will automatically be added as one variable when balancing
#' the downsampling. However, if users want to balance the downsampling solely
#' basing on dataset origin, users have to explicitly set \code{balance =
#' "dataset"}.
#' @param object \linkS4class{liger} object
#' @param balance Character vector of categorical variable names in
#' \code{cellMeta} slot, to subsample \code{maxCells} cells from each
#' combination of all specified variables. Default \code{NULL} samples
#' \code{maxCells} cells from the whole object.
#' @param maxCells Max number of cells to sample from the grouping based on
#' \code{balance}.
#' @param useDatasets Index selection of datasets to include Default
#' \code{NULL} for using all datasets.
#' @param seed Random seed for reproducibility. Default \code{1}.
#' @param returnIndex Logical, whether to only return the numeric index that can
#' subset the original object instead of a subset object. Default \code{FALSE}.
#' @param ... Arguments passed to \code{\link{subsetLiger}}, where
#' \code{cellIdx} is occupied by internal implementation.
#' @return By default, a subset of \linkS4class{liger} \code{object}.
#' Alternatively when \code{returnIndex = TRUE}, a numeric vector to be used
#' with the original object.
#' @export
#' @examples
#' # Subsetting an object
#' pbmc <- downsample(pbmc)
#' # Creating a subsetting index
#' sampleIdx <- downsample(pbmcPlot, balance = "leiden_cluster",
#'                         maxCells = 10, returnIndex = TRUE)
#' plotClusterDimRed(pbmcPlot, cellIdx = sampleIdx)
downsample <- function(
    object,
    balance = NULL,
    maxCells = 1000,
    useDatasets = NULL,
    seed = 1,
    returnIndex = FALSE,
    ...
) {
    set.seed(seed)
    useDatasets <- .checkUseDatasets(object, useDatasets)
    selected <- c()
    if (is.null(balance)) {
        useCells <- which(object$dataset %in% useDatasets)
        maxCells <- min(maxCells, length(useCells))
        selected <- sort(sample(useCells, maxCells))
    } else {
        balance <- unique(c("dataset", balance))
        vars <- .fetchCellMetaVar(object, balance, checkCategorical = TRUE,
                                  droplevels = TRUE, drop = FALSE,
                                  cellIdx = object$dataset %in% useDatasets)

        vars <- vars %>%
            dplyr::group_by_at(.vars = balance) %>%
            dplyr::count()
        for (i in seq(nrow(vars))) {
            comb <- vars[i,]
            name <- names(comb)[seq_along(balance)]
            value <- as.vector(t(as.data.frame(comb)))[seq_along(balance)]
            subscrTxt <- paste0(
                "which(",
                paste0("object$", name, ' == "', value, '"', collapse = " & "),
                ")"
            )
            useCells <- eval(parse(text = subscrTxt))
            if (maxCells < comb[["n"]])
                selected <- c(selected , sample(useCells, maxCells))
            else
                selected <- c(selected, useCells)
        }
        selected <- sort(selected)
    }
    if (isTRUE(returnIndex)) return(selected)
    else return(subsetLiger(object = object, cellIdx = selected, ...))
}

#' `r lifecycle::badge("superseded")` See \code{\link{downsample}}
#' @description This function mainly aims at downsampling datasets to a size
#' suitable for plotting.
#' @param object \linkS4class{liger} object
#' @param slot.use Only create subset from one or more of \code{"rawData"},
#' \code{"normData"} and \code{"scaleData"}. Default \code{NULL} subsets the
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
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @return Subset of \linkS4class{liger} \code{object}.
#' @seealso \code{\link{downsample}}, \code{\link{subsetLiger}},
#' \code{\link{subsetLigerDataset}}
#' @export
readSubset <- function(
        object,
        slot.use = "normData",
        balance = NULL,
        max.cells = 1000,
        chunk = 1000,
        datasets.use = NULL,
        genes.use = NULL,
        rand.seed = 1,
        verbose = getOption("ligerVerbose", TRUE)
) {
    .Deprecated("downsample") # nocov start
    if (!is.null(balance)) balance <- match.arg(balance)
    else balance <- "all"
    downsample(object = object, balance = balance,
               maxCells = max.cells, useDatasets = datasets.use,
               useGenes = genes.use, useSlot = slot.use, seed = rand.seed,
               chunkSize = chunk, verbose = verbose) # nocov end
}

