#' General QC for liger object
#' @description Calculate number of UMIs, number of detected features and
#' percentage of feature subset (e.g. mito) expression per cell.
#' @param object \linkS4class{liger} object with \code{raw.data} available in
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
#' in \code{cell.meta(object)}, as well as expression percentage value for each
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
    # Not using S4 cell.meta() method below because no need to do so
    for (nrn in newResultNames) object[[nrn]] <- 0
    for (d in useDatasets) {
        ld <- dataset(object, d)
        if (isTRUE(verbose))
            message(date(), ' ... calculating QC for dataset "', d, '"')
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
        object@cell.meta[object$dataset == d, newResultNames] <- results$cell
        feature.meta(ld, check = FALSE)$nCell <- results$feature
        datasets(object)[[d]] <- ld
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
        useData = "raw.data",
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
    nUMI <- Matrix::colSums(raw.data(object))
    nonzero <- raw.data(object) > 0
    nGene <- Matrix::colSums(nonzero)
    nCell <- Matrix::rowSums(nonzero)
    results <- data.frame(nUMI = nUMI, nGene = nGene,
                          row.names = colnames(object))
    if (length(featureSubsets) > 0) {
        percentages <- lapply(featureSubsets, function(x) {
            row.idx <- rownames(object) %in% x
            colSums(raw.data(object)[row.idx,]) / colSums(raw.data(object)) * 100
        })
        results <- cbind(results, as.data.frame(percentages))
    }
    list(cell = results, feature = nCell)
}

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
    for (d in useDatasets) {
        if (isTRUE(verbose)) message(date(),
                                     " ... Removing missing in dataset: ", d)
        ld <- dataset(object, d)
        if (rmFeature) {
            featureIdx <- feature.meta(ld)$nCell > 0
        } else {
            featureIdx <- seq(nrow(ld))
        }
        if (rmCell) {
            cellIdx <- colnames(object)[object$dataset == d & object$nGene > 0]
        } else {
            cellIdx <- seq(ncol(ld))
        }
        datasets.new[[d]] <- subsetLigerDataset(ld, featureIdx = featureIdx,
                                                cellIdx = cellIdx,
                                                filenameSuffix = filenameSuffix,
                                                verbose = verbose, ...)
    }
    methods::new(
        "liger",
        datasets = datasets.new,
        cell.meta = cell.meta(object)[object$nGene > 0, , drop = FALSE],
        var.features = character(),
        H.norm = object@H.norm[cellIdx, , drop = FALSE],
        version = packageVersion("rliger")
    )
}

#' Generate violin plot of QC metrics of each dataset
#'
#' By default, a violin plot grouped by dataset variable, overlaid by a box plot
#' would be generated. Options are available for adding jittered dot plot
#' behind.
#' @param object \linkS4class{liger} object, or a data.frame like object of cell
#' metadata, where \code{"dataset"} variable must exist.
#' @param metric Metric to use. Should be found in \code{cell.meta(object)}.
#' @param groupBy Group data by this categorical variable.
#' @param colorBy Color the dot, violin and box plot with this variable.
#' @param dotPlot,violinPlot,boxPlot Whether to add corresponding plot(s).
#' @param dotColor Specify a color to set uniform color to all dots. Default
#' \code{"black"}. Set \code{NULL} then colored by \code{colorBy}.
#' @param dotSize Size of dots.
#' @param violinAlpha,boxAlpha The transparency of violins or boxes,
#' respectively. \code{1} is opaque and \code{0} is transparent.
#' @param legend Whether to show the legend.
#' @param xlab,ylab X and Y axis title text.
#' @param baseSize One-parameter control of all text sizes. Individual text
#' element size can be controled by other size arguments. "Title" sizes are
#' 2 points larger than "text" sizes.
#' @param xTextSize,yTextSize,legendTextSize Size of axis texts and legend text.
#' Default \code{NULL} controls by \code{baseSize}.
#' @param xTitleSize,yTitleSize,legendTitleSize Size of axis titles and legend
#' title. Default \code{NULL} controls by \code{baseSize}.
#' @param ... Arguments passed to \code{plotQCMetric}
#' @export
#' @rdname plotQCMetric
#' @importFrom rlang .data
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
        boxPlot = TRUE,
        boxAlpha = 0,
        legend = TRUE,
        xlab = "Dataset",
        ylab = metric,
        baseSize = 10,
        xTextSize = NULL,
        xTitleSize = NULL,
        yTextSize = NULL,
        yTitleSize = NULL,
        legendTextSize = NULL,
        legendTitleSize = NULL
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

    p <- p + ggplot2::theme_classic()

    xText <- yText <- legendText <- baseSize
    xTitle <- yTitle <- legendTitle <- baseSize + 2
    if (!is.null(xTextSize)) xText <- xTextSize
    if (!is.null(xTitleSize)) xTitle <- xTitleSize
    if (!is.null(yTextSize)) yText <- yTextSize
    if (!is.null(yTitleSize)) yTitle <- yTitleSize
    if (!is.null(legendTextSize)) legendText <- legendTextSize
    if (!is.null(legendTitleSize)) legendTitle <- legendTitleSize

    p <- p + ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = xText),
        axis.title.x = ggplot2::element_text(size = xTitle),
        axis.text.y = ggplot2::element_text(size = yText),
        axis.title.y = ggplot2::element_text(size = yTitle),
        legend.text = ggplot2::element_text(size = legendText),
        legend.title = ggplot2::element_text(size = legendTitle)
    )

    if (!isTRUE(legend))
        p <- p + ggplot2::theme(legend.position = "none")
    if (!is.null(xlab))
        p <- p + xlab(xlab)
    else
        p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank())
    if (!is.null(ylab))
        p <- p + ylab(ylab)
    else
        p <- p + ggplot2::theme(axis.title.y = ggplot2::element_blank())
    p
}

#' @export
#' @rdname plotQCMetric
plotTotalCountViolin <- function(
    object,
    metric = "nUMI",
    ylab = "Total Counts",
    ...
) {
    plotQCMetric(object, metric = "nUMI", ylab = ylab, ...)
}

#' @export
#' @rdname plotQCMetric
plotGeneDetectedViolin <- function(
    object,
    metric = "nGene",
    ylab = "Number of Genes Detected",
    ...
) {
    plotQCMetric(object, metric = "nGene", ylab = ylab, ...)
}
