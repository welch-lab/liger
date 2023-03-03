#' General QC for liger object
#' @description Calculate number of UMIs, number of detected features and
#' percentage of feature subset (e.g. mito) expression per cell.
#' @param object \linkS4class{liger} object with \code{rawData} available in
#' each \linkS4class{ligerDataset} embedded
#' @param mito,ribo,hemo Whether to calculate the expression percentage of
#' mitochondrial, ribosomal or hemoglobin genes, respectively. Default
#' \code{TRUE}.
#' @param features Feature names matching the feature subsets that users want to
#' calculate the expression percentage with. A vector for a single subset, or a
#' named list for multiple subset. Default \code{NULL}.
#' @param pattern Regex patterns for matching the feature subsets that users
#' want to calculate the expression percentage with. A vector for a single
#' subset, or a named list for multiple subset. Default \code{NULL}.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to be included for QC. Default
#' \code{names(object)} (i.e. all datasets).
#' @param chunkSize Integer number of cells to include in a chunk when working
#' on HDF5 based dataset. Default \code{1000}
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{TRUE}.
#' @return Updated \code{object} with \code{nUMI}, \code{nGene} updated
#' in \code{cellMeta(object)}, as well as expression percentage value for each
#' feature subset.
#' @export
runGeneralQC <- function(
        object,
        mito = TRUE,
        ribo = TRUE,
        hemo = TRUE,
        features = NULL,
        pattern = NULL,
        useDatasets = names(object),
        chunkSize = 1000,
        verbose = TRUE
) {
    useDatasets <- .checkUseDatasets(object, useDatasets)
    # Process the the two arguments all into one named list of feature names
    # before exactly calculate the percentage
    featureSubsets <- list()
    allFeatures <- unique(unlist(lapply(datasets(object), rownames),
                                 use.names = FALSE))

    # Work on the presets
    if (isTRUE(mito))
        featureSubsets$mito <- grep("^MT-", allFeatures, value = TRUE)
    if (isTRUE(ribo))
        featureSubsets$ribo <- grep("^RP[SL]", allFeatures, value = TRUE)
    if (isTRUE(hemo))
        featureSubsets$hemo <- grep("^HB[^(P)]", allFeatures, value = TRUE)

    # Then process the user specified gene sets
    if (!is.null(features)) {
        if (is.list(features)) {
            featureSubsets <- c(featureSubsets, features)
        } else if (is.vector(features)) {
            featureSubsets[["featureSubset_name"]] <- features
        }
    }
    if (!is.null(pattern)) {
        if (is.list(pattern)) {
            pattern <- lapply(pattern, function(x) {
                grep(x, allFeatures, value = TRUE)
            })
            featureSubsets <- c(featureSubsets, pattern)
        } else if (is.vector(pattern)) {
            pattern <- grep(pattern, allFeatures, value = TRUE)
            featureSubsets[["featureSubsets_pattern"]] <- pattern
        }
    }

    # Start calculation on each dataset
    newResultNames <- c("nUMI", "nGene", names(featureSubsets))
    # Not using S4 cellMeta() method below because no need to do so
    for (nrn in newResultNames) object[[nrn]] <- 0
    for (d in useDatasets) {
        ld <- dataset(object, d)
        if (isTRUE(verbose)) .log('calculating QC for dataset "', d, '"')
        if (isH5Liger(ld))
            results <- runGeneralQC.h5(
                ld,
                featureSubsets = featureSubsets,
                chunkSize = chunkSize,
                verbose = verbose
            )
        else
            results <- runGeneralQC.Matrix(
                ld,
                featureSubsets = featureSubsets,
                verbose = verbose
            )
        object@cellMeta[object$dataset == d, newResultNames] <- results$cell
        featureMeta(ld, check = FALSE)$nCell <- results$feature
        datasets(object, check = FALSE)[[d]] <- ld
    }

    return(object)
}

#' Calculate general QC on H5 based ligerDataset object
#' @param object ligerDataset object
#' @param featureSubsets Named list passed from \code{runGeneralQC}
#' @param chunkSize Integer
#' @return data.frame
#' @noRd
runGeneralQC.h5 <- function(
        object,
        featureSubsets = NULL,
        chunkSize = 1000,
        verbose = TRUE) {
    allFeatures <- rownames(object)
    # Initialize results
    cell <- data.frame(row.names = colnames(object))
    cell$nUMI <- 0
    cell$nGene <- 0
    for (i in names(featureSubsets)) {
        cell[[i]] <- 0
    }
    nCell <- rep(0, nrow(object))
    rowIndices <- lapply(featureSubsets, function(x) allFeatures %in% x)

    # Calculate in only one iteration
    H5Apply(
        object,
        init = list(cell = cell, feature = nCell),
        useData = "rawData",
        chunkSize = chunkSize,
        verbose = verbose,
        FUN = function(chunk, sparseXIdx, cellIdx, values) {
            nUMI <- colSums(chunk)
            values$cell$nUMI[cellIdx] <- nUMI
            nonzero <- chunk > 0
            values$cell$nGene[cellIdx] <- colSums(nonzero)
            for (fs in names(rowIndices)) {
                values$cell[[fs]][cellIdx] <-
                    colSums(chunk[rowIndices[[fs]],]) / nUMI * 100
            }
            values$feature <- values$feature + Matrix::rowSums(nonzero)
            return(values)
        }
    )
}

#' Calculate general QC on in memory matrix based ligerDataset object
#' @param object ligerDataset object
#' @param featureSubsets  Named list passed from \code{runGeneralQC}
#' @return data.frame
#' @noRd
runGeneralQC.Matrix <- function(
        object,
        featureSubsets = NULL,
        verbose = TRUE) {
    nUMI <- Matrix::colSums(rawData(object))
    nonzero <- rawData(object) > 0
    nGene <- Matrix::colSums(nonzero)
    nCell <- Matrix::rowSums(nonzero)
    results <- data.frame(nUMI = nUMI, nGene = nGene,
                          row.names = colnames(object))
    if (length(featureSubsets) > 0) {
        percentages <- lapply(featureSubsets, function(x) {
            row.idx <- rownames(object) %in% x
            colSums(rawData(object)[row.idx,]) / colSums(rawData(object)) * 100
        })
        results <- cbind(results, as.data.frame(percentages))
    }
    list(cell = results, feature = nCell)
}

#' Remove missing cells or features from liger object
#' @param object \linkS4class{liger} object
#' @param orient Choose to remove non-expressing features (\code{"feature"}),
#' empty barcodes (\code{"cell"}), or both of them (\code{"both"}). Default
#' \code{"both"}.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to be processed. Default
#' \code{NULL} removes empty entries from all datasets.
#' @param filenameSuffix When subsetting H5 based datasets to new H5 files, this
#' suffix will be added to all the filenames. Default \code{"removeMissing"}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{TRUE}.
#' @param ... Arguments passed to \code{\link{subsetLigerDataset}}
#' @return Updated (subset) \code{object}.
#' @export
removeMissing <- function(
        object,
        orient = c("both", "feature", "cell"),
        useDatasets = NULL,
        filenameSuffix = "removeMissing",
        verbose = TRUE,
        ...
) {
    if (!inherits(object, "liger")) {
        stop("Please use a `liger` object.")
    }
    orient <- match.arg(orient)
    useDatasets <- .checkUseDatasets(object, useDatasets)
    rmFeature <- ifelse(orient %in% c("both", "feature"), TRUE, FALSE)
    rmCell <- ifelse(orient %in% c("both", "cell"), TRUE, FALSE)
    datasets.new <- list()
    subsetted <- c()
    for (d in useDatasets) {
        ld <- dataset(object, d)
        if (rmFeature) {
            featureIdx <- which(featureMeta(ld)$nCell > 0)
        } else {
            featureIdx <- seq(nrow(ld))
        }
        if (length(featureIdx) == nrow(ld)) rmFeature <- FALSE
        if (rmCell) {
            cellIdx <- colnames(object)[object$dataset == d & object$nGene > 0]
            cellIdx <- which(colnames(ld) %in% cellIdx)
        } else {
            cellIdx <- seq(ncol(ld))
        }
        if (length(cellIdx) == ncol(ld)) rmCell <- FALSE
        subsetted <- c(subsetted, any(c(rmFeature, rmCell)))
        if (any(c(rmFeature, rmCell))) {
            if (isTRUE(verbose)) .log("Removing missing in dataset: ", d)
            datasets.new[[d]] <- subsetLigerDataset(
                ld,
                featureIdx = featureIdx,
                cellIdx = cellIdx,
                filenameSuffix = filenameSuffix,
                verbose = verbose, ...
            )
        } else {
            datasets.new[[d]] <- ld
        }
    }
    if (any(subsetted)) {
        methods::new(
            "liger",
            datasets = datasets.new,
            cellMeta = cellMeta(object, cellIdx = object$nGene > 0,
                                  drop = FALSE),
            varFeatures = character(),
            H.norm = object@H.norm[cellIdx, , drop = FALSE]
        )
    } else {
        object
    }
}

#' Generate violin plot of QC metrics of each dataset
#'
#' By default, a violin plot grouped by dataset variable, overlaid by a box plot
#' would be generated. Options are available for adding jittered dot plot
#' behind.
#' @param object \linkS4class{liger} object, or a data.frame like object of cell
#' metadata, where \code{"dataset"} variable must exist.
#' @param metric Metric to use. Should be found in \code{cellMeta(object)}.
#' @param groupBy Group data by this categorical variable.
#' @param colorBy Color the dot, violin and box plot with this variable.
#' @param dotPlot,violinPlot,boxPlot Whether to add corresponding plot(s).
#' @param dotColor Specify a color to set uniform color to all dots. Default
#' \code{"black"}. Set \code{NULL} then colored by \code{colorBy}.
#' @param dotSize Size of dots.
#' @param violinAlpha,boxAlpha The transparency of violins or boxes,
#' respectively. \code{1} is opaque and \code{0} is transparent.
#' @param ... Theme setting parameters. Check \code{\link{.ggplotLigerTheme}}
#' for more options and details.
#' @export
#' @rdname plotQCMetric
plotQCMetric <- function(
        object,
        metric,
        groupBy = "dataset",
        colorBy = NULL,
        dotPlot = FALSE,
        dotColor = "black",
        dotSize = 1,
        violinPlot = TRUE,
        violinAlpha = .8,
        boxPlot = FALSE,
        boxAlpha = 0,
        xlab = groupBy,
        ylab = metric,
        ...
) {
    if (is.null(colorBy))
        p <- ggplot2::ggplot(object, ggplot2::aes(x = .data[[groupBy]],
                                                  y = .data[[metric]]))
    else
        p <- ggplot2::ggplot(object, ggplot2::aes(x = .data[[groupBy]],
                                                  y = .data[[metric]],
                                                  color = .data[[colorBy]]))
    if (isTRUE(dotPlot)) {
        if (is.null(dotColor))
            p <- p + ggplot2::geom_jitter(size = dotSize,
                                          stroke = 0.1,
                                          height = 0)
        else
            p <- p + ggplot2::geom_jitter(color = dotColor,
                                          size = dotSize,
                                          stroke = 0.1,
                                          height = 0)
    }

    if (isTRUE(violinPlot)) p <- p + ggplot2::geom_violin()

    if (isTRUE(boxPlot)) p <- p + ggplot2::geom_boxplot(alpha = boxAlpha)

    p <- .ggplotLigerTheme(p, xlab = xlab, ylab = ylab, ...)
    p
}

#' @export
#' @rdname plotQCMetric
plotTotalCountViolin <- function(
    object,
    groupBy = "dataset",
    ...
) {
    plotCellViolin(object, y = "nUMI", groupBy = groupBy,
                   ylab = "Total counts", ...)
}

#' @export
#' @rdname plotQCMetric
plotGeneDetectedViolin <- function(
    object,
    groupBy = "dataset",
    ...
) {
    callArgs <- names(rlang::call_args(match.call()))

    plotCellViolin(object, y = "nGene", groupBy = groupBy,
                   ylab = "Number of Genes Detected", ...)
}
