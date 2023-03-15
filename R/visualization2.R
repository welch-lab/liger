#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Scatter Plots of DimRed ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Generate Dimensionality Reduction Plot with Coloring
#' @description some text
#' @param object A \linkS4class{liger} object.
#' @param useCluster Name of variable in \code{cellMeta} slot. Default
#' \code{"leiden_cluster"}.
#' @param useDimRed Name of the variable storing dimensionality reduction result
#' in the \code{cellMeta} slot. Default \code{"UMAP"}.
#' @param combinePlots Logical, whether to utilize
#' \code{\link[cowplot]{plot_grid}} to combine multiple plots into one. Default
#' \code{TRUE} returns combined ggplot. \code{FALSE} returns a list of ggplot.
#' @param features,factors Name of genes or index of factors that need to be
#' visualized.
#' @param log Logical. Whether to log transform the normalized expression of
#' genes. Default \code{TRUE}.
#' @param scaleFactor Number to be multiplied with the normalized expression of
#' genes before log transformation. Default \code{1e4}. \code{NULL} for not
#' scaling.
#' @param zeroAsNA Logical, whether to swap all zero values to \code{NA} so
#' \code{naColor} will be used to represent non-expressing features. Default
#' \code{TRUE}.
#' @param trimHigh Number for highest cut-off to limit the outliers. Factor
#' loading above this value will all be trimmed to this value. Default
#' \code{0.03}.
#' @param colorPalette Name of viridis palette. See
#' \code{\link[viridisLite]{viridis}} for options. Default \code{"C"} ("plasma")
#' for gene expression and \code{"D"} ("viridis") for factor loading.
#' @param ... Additional graphic setting arguments passed to
#' \code{\link{plotCellScatter}}.
#' @return ggplot object when only one feature (e.g. cluster variable, gene,
#' factor) is set. List object when multiple of those are specified.
#' @seealso Please refer to \code{\link{plotCellScatter}},
#' \code{\link{.ggScatter}}, \code{\link{.ggplotLigerTheme}} for additional
#' graphic setting
#' @rdname plotDimRed
#' @export
#' @examples
#' data("pbmcPlot", package = "rliger")
#' plotClusterDimRed(pbmcPlot)
#' plotDatasetDimRed(pbmcPlot)
#' plotByDatasetAndCluster(pbmcPlot)
#' plotGeneDimRed(pbmcPlot, varFeatures(pbmcPlot)[1])
#' plotFactorDimRed(pbmcPlot, 2)
plotClusterDimRed <- function(
        object,
        useCluster = "leiden_cluster",
        useDimRed = "UMAP",
        ...) {
    xVar <- paste0(useDimRed, ".1")
    yVar <- paste0(useDimRed, ".2")
    plotCellScatter(object, x = xVar, y = yVar, colorBy = useCluster,
                    slot = "cellMeta", dotOrder = "shuffle", ...)
}

#' @rdname plotDimRed
#' @export
plotDatasetDimRed <- function(
        object,
        useDimRed = "UMAP",
        ...) {
    xVar <- paste0(useDimRed, ".1")
    yVar <- paste0(useDimRed, ".2")
    plotCellScatter(object, x = xVar, y = yVar, colorBy = "dataset",
                    slot = "cellMeta", labelText = FALSE,
                    dotOrder = "shuffle", ...)
}

#' @rdname plotDimRed
#' @export
plotByDatasetAndCluster <- function(
        object,
        useDimRed = "UMAP",
        useCluster = "leiden_cluster",
        combinePlots = TRUE,
        ...
) {
    plot <- list(
        dataset = plotDatasetDimRed(object, useDimRed = useDimRed, ...),
        cluster = plotClusterDimRed(object, useCluster = useCluster,
                                    useDimRed = useDimRed, ...)
    )
    if (isTRUE(combinePlots)) {
        plot <- cowplot::plot_grid(plotlist = plot, nrow = 1,
                                   align = "h", axis = "tblr")
    }
    return(plot)
}

#' @rdname plotDimRed
#' @export
plotGeneDimRed <- function(
        object,
        features,
        useDimRed = "UMAP",
        log = TRUE,
        scaleFactor = 1e4,
        zeroAsNA = TRUE,
        colorPalette = "C",
        ...
) {
    xVar <- paste0(useDimRed, ".1")
    yVar <- paste0(useDimRed, ".2")
    scaleFunc <- function(x) {
        if (!is.null(scaleFactor)) x <- scaleFactor*x
        if (isTRUE(log)) x <- log2(x + 1)
        x
    }
    plotCellScatter(object, x = xVar, y = yVar, colorBy = features,
                    slot = "normData", colorByFunc = scaleFunc,
                    dotOrder = "ascending", zeroAsNA = zeroAsNA,
                    colorPalette = colorPalette, ...)
}

#' @rdname plotDimRed
#' @export
plotFactorDimRed <- function(
        object,
        factors,
        useDimRed = "UMAP",
        trimHigh = 0.03,
        zeroAsNA = TRUE,
        colorPalette = "D",
        ...
) {
    xVar <- paste0(useDimRed, ".1")
    yVar <- paste0(useDimRed, ".2")
    plotCellScatter(object, x = xVar, y = yVar, colorBy = factors,
                    slot = "H.norm", dotOrder = "ascending",
                    trimHigh = trimHigh, zeroAsNA = zeroAsNA,
                    colorPalette = colorPalette, ...)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Violin Plots ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Visualize gene expression or cell metadata with violin plot
#' @param object A \linkS4class{liger} object.
#' @param gene Character gene names.
#' @param byDataset Logical, whether the violin plot should be splitted by
#' dataset. Default \code{TRUE}.
#' @param groupBy Names of available categorical variable in \code{cellMeta}
#' slot. Use \code{FALSE} for no grouping. Default \code{NULL} looks clustering
#' result but will not group if no clustering found.
#' @param ... Additional arguments passed to \code{\link{plotCellViolin}}.
#' @return ggplot if using a single gene and not splitting by dataset.
#' Otherwise, list of ggplot.
#' @export
#' @rdname plotViolin
#' @examples
#' data("pbmcPlot", package = "rliger")
#' plotGeneViolin(pbmcPlot, varFeatures(pbmcPlot)[1],
#'                groupBy = "leiden_cluster")
#' data("pbmc", package = "rliger")
#' plotTotalCountViolin(pbmc)
#' plotGeneDetectedViolin(pbmc, dot = TRUE, box = TRUE, colorBy = "dataset")
plotGeneViolin <- function(
        object,
        gene,
        byDataset = TRUE,
        groupBy = NULL,
        ...
) {
    splitBy <- NULL
    if (isTRUE(byDataset)) splitBy <- "dataset"

    if (is.null(groupBy)) {
        if ("leiden_cluster" %in% names(cellMeta(object)))
            groupBy <- "leiden_cluster"
        else if ("louvain_cluster" %in% names(cellMeta(object)))
            groupBy <- "louvain_cluster"
        else if ("H.norm_cluster" %in% names(cellMeta(object)))
            groupBy <- "H.norm_cluster"
    } else if (isFALSE(groupBy)) {
        groupBy <- NULL
    }

    plotList <- plotCellViolin(
        object,
        y = gene,
        slot = "normData",
        yFunc = function(x) log2(10000*x + 1),
        groupBy = groupBy,
        splitBy = splitBy,
        ...
    )

    if (!is.null(splitBy)) {
        datasetNames <- names(object)
        plotTitles <- rep(datasetNames, length(gene))
        for (i in seq_along(plotList)) {
            plotList[[i]] <- plotList[[i]] +
                ggplot2::ggtitle(paste0("Dataset: ", plotTitles[i]))
        }
    }

    return(plotList)
}

#' @export
#' @rdname plotViolin
plotTotalCountViolin <- function(
        object,
        groupBy = "dataset",
        ...
) {
    plotCellViolin(object, y = "nUMI", groupBy = groupBy,
                   ylab = "Total counts", ...)
}

#' @export
#' @rdname plotViolin
plotGeneDetectedViolin <- function(
        object,
        groupBy = "dataset",
        ...
) {
    plotCellViolin(object, y = "nGene", groupBy = groupBy,
                   ylab = "Number of Genes Detected", ...)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Proportion #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Visualize proportion across two categorical variables
#' @description \code{plotProportionBar} creates bar plots comparing the
#' cross-category proportion. \code{plotProportionDot} creates dot plots.
#' \code{plotClusterProportions} has variable pre-specified and calls the dot
#' plot. \code{plotProportion} produces a combination of both bar plots and dot
#' plot.
#'
#' Having package "ggrepel" installed can help adding tidier percentage
#' annotation on the pie chart.
#' @param object A \linkS4class{liger} object.
#' @param class1,class2 Each should be a single name of a categorical variable
#' available in \code{cellMeta} slot. Number of cells in each categories in
#' \code{class2} will be served as the denominator when calculating proportions.
#' @param method For bar plot, choose whether to draw \code{"stack"} or
#' \code{"group"} bar plot. Default \code{"stack"}.
#' @param showLegend,panelBorder,... ggplot theme setting arguments passed to
#' \code{\link{.ggplotLigerTheme}}.
#' @param inclRev Logical, for barplot, whether to reverse the specification for
#' \code{class1} and \code{class2} and produce two plots. Default \code{FALSE}.
#' @param combinePlot Logical, whether to combine the two plots with
#' \code{\link[cowplot]{plot_grid}} when two plots are created. Default
#' \code{TRUE}.
#' @param useCluster For \code{plotClusterProportions}. Same as \code{class1}
#' while \code{class2} is hardcoded with \code{"dataset"}.
#' @param labelSize,labelColor Settings on pie chart percentage label. Default
#' \code{4} and \code{"white"}.
#' @param return.plot \bold{defuncted}.
#' @return ggplot or list of ggplot
#' @rdname plotProportion
#' @export
#' @examples
#' data("pbmcPlot", package = "rliger")
#' plotProportion(pbmcPlot)
#' plotProportionBar(pbmcPlot, method = "group")
#' plotProportionPie(pbmcPlot)
plotProportion <- function(
        object,
        class1 = "leiden_cluster",
        class2 = "dataset",
        method = c("stack", "group", "pie"),
        ...
) {
    method <- match.arg(method)
    p1 <- plotProportionDot(object, class1 = class1, class2 = class2, ...)
    if (method %in% c("stack", "group")) {
        p2 <- plotProportionBar(object, class1 = class1, class2 = class2,
                                inclRev = TRUE, combinePlot = TRUE,
                                method = method, ...)
        cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(1, 2))
    } else {
        p2 <- plotProportionPie(object, class1 = class1, class2 = class2, ...)
        p3 <- plotProportionPie(object, class1 = class2, class2 = class1, ...)
        bottom <- cowplot::plot_grid(p2, p3, nrow = 1)
        cowplot::plot_grid(p1, bottom, nrow = 2, rel_heights = c(1, 2))
    }

}

#' @rdname plotProportion
#' @export
plotProportionDot <- function(
        object,
        class1 = "leiden_cluster",
        class2 = "dataset",
        showLegend = FALSE,
        panelBorder = TRUE,
        ...
) {
    if (length(class1) != 1 ||
        length(class2) != 1)
        stop("`class1` and `class2` must be name of one categorical variable ",
             "in `cellMeta` slot.")
    vars <- .fetchCellMetaVar(object, c(class1, class2),
                              checkCategorical = TRUE)
    freq <- table(vars)
    for (i in seq(ncol(freq))) freq[,i] <- freq[,i] / sum(freq[,i])
    freqDF <- as.data.frame(freq)
    colnames(freqDF)[3] <- "Proportion"
    p <- ggplot2::ggplot(freqDF,
                    ggplot2::aes(x = .data[[class1]], y = .data[[class2]],
                                 size = .data[["Proportion"]],
                                 color = .data[[class1]])) +
        ggplot2::geom_point() +
        ggplot2::theme(
            axis.line = ggplot2::element_blank(),
            plot.margin = grid::unit(c(0, 0, 0, 0), "cm")
        ) +
        ggplot2::scale_y_discrete(position = "right") +
        ggplot2::coord_fixed(ratio = 0.5)
    .ggplotLigerTheme(p, showLegend = showLegend, panelBorder = panelBorder,
                      ...)
}

#' @rdname plotProportion
#' @export
plotProportionBar <- function(
        object,
        class1 = "leiden_cluster",
        class2 = "dataset",
        method = c("stack", "group"),
        inclRev = FALSE,
        panelBorder = TRUE,
        combinePlot = TRUE,
        ...
) {
    if (length(class1) != 1 ||
        length(class2) != 1)
        stop("`class1` and `class2` must be name of one categorical variable ",
             "in `cellMeta` slot.")
    method <- match.arg(method)
    vars <- .fetchCellMetaVar(object, c(class1, class2),
                              checkCategorical = TRUE)
    freq <- table(vars)
    for (i in seq(ncol(freq))) freq[,i] <- freq[,i] / sum(freq[,i])
    freqDF <- as.data.frame(freq)
    colnames(freqDF)[3] <- "Proportion"
    p <- ggplot2::ggplot(freqDF, ggplot2::aes(x = .data[[class2]],
                                         y = .data[["Proportion"]],
                                         fill = .data[[class1]])) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::scale_x_discrete(expand = c(0, 0))
    if (method == "stack")
        p <- p + ggplot2::geom_col(position = "fill", width = 0.95)
    else if (method == "group")
        p <- p + ggplot2::geom_bar(position = "dodge", stat = "identity")
    p <- .ggplotLigerTheme(p, panelBorder = panelBorder, ...)

    if (isTRUE(inclRev)) {
        p2 <- plotProportionBar(object,
                                class1 = class2,
                                class2 = class1,
                                method = method,
                                inclRev = FALSE,
                                panelBorder = panelBorder, ...)
        if (isTRUE(combinePlot))
            return(cowplot::plot_grid(p, p2, align = "h", axis = "tblr"))
        else return(list(p, p2))
    } else {
        return(p)
    }
}

#' @rdname plotProportion
#' @export
plotClusterProportions <- function(
        object,
        useCluster = "leiden_cluster",
        return.plot = FALSE,
        ...
) {
    .deprecateArgs(defunct = "return.plot")
    plotProportionDot(object, class1 = useCluster, class2 = "dataset", ...)
}

#' @rdname plotProportion
#' @export
plotProportionPie <- function(
        object,
        class1 = "leiden_cluster",
        class2 = "dataset",
        labelSize = 4,
        labelColor = "white",
        ...
) {
    df <- .fetchCellMetaVar(object, c(class1, class2), drop = FALSE,
                            checkCategorical = TRUE) %>%
        dplyr::group_by(.data[[class2]]) %>%
        dplyr::count(.data[[class1]]) %>%
        dplyr::mutate(proportion = .data[["n"]] / sum(.data[["n"]])) %>%
        dplyr::mutate(cumsumProp = cumsum(.data[["proportion"]]) -
                          0.5*.data[["proportion"]]) %>%
        dplyr::ungroup()

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["cumsumProp"]],
                                     y = .data[[class2]],
                                     fill = .data[[class1]],
                                     label = sprintf("%1.1f%%",
                                                     100*.data[["proportion"]]),
                                     width = .data[["proportion"]])) +
        ggplot2::geom_tile(colour = "white", size = 0.3) +
        ggplot2::coord_polar()
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggplot2::geom_text(
            size = labelSize, color = labelColor,
            position = ggplot2::position_nudge(y = 0.25)
        )
    } else {
        p <- p + ggrepel::geom_text_repel(
            size = labelSize, color = labelColor, force = 0.001, max.overlaps = 4,
            position = ggplot2::position_nudge(y = 0.25)
        )
    }
    .ggplotLigerTheme(p, ...) +
        ggplot2::theme(
            axis.line = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        )
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Volcano plot ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create volcano plot for Wilcoxon test result
#' @description \code{plotVolcano} is a simple implementation and shares
#' most of arguments with other rliger plotting functions.
#' \code{plotEnhancedVolcano} is a wrapper function of
#' \code{\link[EnhancedVolcano]{EnhancedVolcano}}, which has provides
#' substantial amount of arguments for graphical control. However, that requires
#' the installation of package "EnhancedVolcano".
#' @rdname plotVolcano
#' @param result Data frame table returned by \code{\link{runWilcoxon}}
#' @param group Selection of one group available from \code{result$group}
#' @param logFCThresh Number for the threshold on the absolute value of the log2
#' fold change statistics. Default \code{1}.
#' @param padjThresh Number for the threshold on the adjusted p-value
#' statistics. Default \code{0.01}.
#' @param labelTopN Number of top differential expressed features to be labeled
#' on the top of the dots. Default \code{20}.
#' @param dotSize,dotAlpha Numbers for universal aesthetics control of dots.
#' Default \code{2} and \code{0.8}.
#' @param legendPosition Text indicating where to place the legend. Choose from
#' \code{"top"}, \code{"bottom"}, \code{"left"} or \code{"right"}. Default
#' \code{"top"}.
#' @param labelSize Size of labeled top features and line annotations. Default
#' \code{4}.
#' @param ... For \code{plotVolcano}, more theme setting arguments passed to
#' \code{\link{.ggplotLigerTheme}}. For \code{plotEnhancedVolcano}, arguments
#' passed to \code{\link[EnhancedVolcano]{EnhancedVolcano}}.
#' @return ggplot
plotVolcano <- function(
        result,
        group,
        logFCThresh = 1,
        padjThresh = 0.01,
        labelTopN = 20,
        dotSize = 2,
        dotAlpha = 0.8,
        legendPosition = "top",
        labelSize = 4,
        ...
) {
    if (!group %in% result$group) {
        stop("Selected group does not exist in `result`.")
    }
    result <- result[result$group == group, ]
    result <- result[order(abs(result$logFC), decreasing = TRUE), ]
    rownames(result) <- result$Gene
    # Prepare for coloring that shows the filtering
    result$Significance <- "Not significant"
    result$Significance[abs(result$logFC) > logFCThresh] <- "logFC"
    result$Significance[result$padj < padjThresh] <- "padj"
    result$Significance[abs(result$logFC) > logFCThresh &
                            result$padj < padjThresh] <- "padj & logFC"
    result$Significance <- factor(result$Significance,
                                  levels = c("Not significant",
                                             "logFC", "padj", "padj & logFC"))
    result$padj[result$padj == 0] <- min(result$padj[result$padj > 0]) / 10
    result$padj <- -log10(result$padj)
    # Prepare for Top result text labeling
    passIdx <- result$Significance == "padj & logFC"
    result$label <- NA
    if (!is.null(labelTopN) && !isFALSE(labelTopN)) {
        labelTopN <- min(labelTopN, length(which(passIdx)))
        if (labelTopN > 0) {
            labelIdx <- which(passIdx)[seq(labelTopN)]
            result$label[labelIdx] <- result$feature[labelIdx]
        }
    }
    # Prepare for lines that mark the cutoffs
    vlineLab <- data.frame(
        X = c(-logFCThresh, logFCThresh)
    )
    hlineLab <- data.frame(
        Y = c(-log10(padjThresh))
    )
    p <- .ggScatter(result, x = "logFC", y = "padj",
                    colorBy = "Significance", zeroAsNA = FALSE,
                    labelBy = "label",
                    xlab = "Log2 Fold Change",
                    ylab = "-log10 Adjusted P-value",
                    colorValues = c("black", "#ef2301", "#416ae1", "#238b22"),
                    legendPosition = legendPosition, ...) +
        ggplot2::xlim(-max(abs(result$logFC)), max(abs(result$logFC))) +
        ggplot2::geom_vline(data = vlineLab,
                            mapping = ggplot2::aes(xintercept = .data[["X"]]),
                            linetype = "longdash") +
        ggplot2::geom_hline(data = hlineLab,
                            ggplot2::aes(yintercept = .data[["Y"]]),
                            linetype = "longdash") +
        ggplot2::annotate("text",
                          x = c(logFCThresh + 3, -logFCThresh - 3),
                          y = c(-10, -10),
                          label = paste0(c("higher ", "lower "),
                                         "log2FC cutoff: ",
                                         c(logFCThresh, -logFCThresh)),
                          size = labelSize) +
        ggplot2::annotate("text",
                          x = -max(abs(result$logFC)) + 2,
                          y = 10,
                          label = paste("p-adj cutoff:", padjThresh),
                          size = labelSize)
    return(p)
}

#' @rdname plotVolcano
#' @export
plotEnhancedVolcano <- function(
        result,
        group,
        ...
) {
    if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
        stop("Package \"EnhancedVolcano\" needed for this function to work. ",
             "Please install it by command:\n",
             "BiocManager::install('EnhancedVolcano')",
             call. = FALSE)
    }
    result <- result[result$group == group, ]
    EnhancedVolcano::EnhancedVolcano(
        toptable = result,
        lab = result$feature, x = "logFC", y = "padj",
        ...
    )
}
