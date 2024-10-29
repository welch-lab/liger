#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Scatter Plots of DimRed ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname plotDimRed
#' @param object A \linkS4class{liger} object.
#' @param useCluster Name of variable in \code{cellMeta(object)}. Default
#' \code{NULL} uses default cluster.
#' @param useDimRed Name of the variable storing dimensionality reduction result
#' in the \code{cellMeta(object)}. Default \code{NULL} use default dimRed.
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
#' @return ggplot object when only one feature (e.g. cluster variable, gene,
#' factor) is set. List object when multiple of those are specified.
#' @export
#' @examples
#' plotClusterDimRed(pbmcPlot)
#' plotDatasetDimRed(pbmcPlot)
#' plotByDatasetAndCluster(pbmcPlot)
#' plotGeneDimRed(pbmcPlot, varFeatures(pbmcPlot)[1])
#' plotFactorDimRed(pbmcPlot, 2)
plotClusterDimRed <- function(
        object,
        useCluster = NULL,
        useDimRed = NULL,
        ...) {
    useCluster <- useCluster %||% object@uns$defaultCluster
    plotDimRed(object, colorBy = useCluster, useDimRed = useDimRed,
               slot = "cellMeta", dotOrder = "shuffle", ...)
}

#' @rdname plotDimRed
#' @export
plotDatasetDimRed <- function(
        object,
        useDimRed = NULL,
        ...) {
    plotDimRed(object, colorBy = "dataset", slot = "cellMeta",
               useDimRed = useDimRed, labelText = FALSE,
               dotOrder = "shuffle", ...)
}

#' @rdname plotDimRed
#' @export
plotByDatasetAndCluster <- function(
        object,
        useDimRed = NULL,
        useCluster = NULL,
        combinePlots = TRUE,
        ...
) {
    plot <- list(
        dataset = plotDatasetDimRed(
            object, useDimRed = useDimRed, ...
        ),
        cluster = plotClusterDimRed(
            object, useCluster = useCluster, useDimRed = useDimRed, ...
        )
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
        useDimRed = NULL,
        log = TRUE,
        scaleFactor = 1e4,
        zeroAsNA = TRUE,
        colorPalette = "C",
        ...
) {
    scaleFunc <- function(x) {
        if (!is.null(scaleFactor)) x <- scaleFactor*x
        if (isTRUE(log)) x <- log2(x + 1)
        x
    }
    plotDimRed(object, colorBy = features, useDimRed = useDimRed,
               slot = "normData", colorByFunc = scaleFunc,
               dotOrder = "ascending", zeroAsNA = zeroAsNA,
               colorPalette = colorPalette, ...)
}

#' @rdname plotDimRed
#' @export
plotPeakDimRed <- function(
        object,
        features,
        useDimRed = NULL,
        log = TRUE,
        scaleFactor = 1e4,
        zeroAsNA = TRUE,
        colorPalette = "C",
        ...
) {
    scaleFunc <- function(x) {
        if (!is.null(scaleFactor)) x <- scaleFactor*x
        if (isTRUE(log)) x <- log2(x + 1)
        x
    }
    plotDimRed(
        object, useDimRed = useDimRed, colorBy = features, slot = "normPeak",
        colorByFunc = scaleFunc, dotOrder = "ascending", zeroAsNA = zeroAsNA,
        colorPalette = colorPalette, ...
    )
}

#' @rdname plotDimRed
#' @export
plotFactorDimRed <- function(
        object,
        factors,
        useDimRed = NULL,
        trimHigh = 0.03,
        zeroAsNA = TRUE,
        colorPalette = "D",
        ...
) {
    plotDimRed(object, colorBy = factors, useDimRed = useDimRed,
               slot = "H.norm", dotOrder = "ascending",
               trimHigh = trimHigh, zeroAsNA = zeroAsNA,
               colorPalette = colorPalette, ...)
}

#' Comprehensive group splited cluster plot on dimension reduction with
#' proportion
#' @description
#' This function produces combined plot on group level (e.g. dataset, other
#' metadata variable like biological conditions). Scatter plot of dimension
#' reduction with cluster labeled is generated per group. Furthermore, a stacked
#' barplot of cluster proportion within each group is also combined with the
#' subplot of each group.
#' @param object A \linkS4class{liger} object with dimension reduction, grouping
#' variable and cluster assignment in \code{cellMeta(object)}.
#' @param useGroup Variable name of the group division in metadata. Default
#' \code{"dataset"}.
#' @param useCluster Name of variable in \code{cellMeta(object)}. Default
#' \code{NULL} uses default cluster.
#' @param useDimRed Name of the variable storing dimensionality reduction result
#' in \code{cellMeta(object)}. Default \code{NULL} use default dimRed.
#' @param combinePlot Whether to return combined plot. Default \code{TRUE}. If
#' \code{FALSE}, will return a list containing only the scatter plots.
#' @param droplevels Logical, whether to perform \code{\link{droplevels}()} on
#' the selected grouping variable. Default \code{TRUE} will not show groups that
#' are listed as categories but do not indeed have any cells.
#' @param relHeightMainLegend Relative heights of the main combination panel and
#' the legend at the bottom. Must be a numeric vector of 2 numbers. Default
#' \code{c(5, 1)}.
#' @param relHeightDRBar Relative heights of the scatter plot and the barplot
#' within each subpanel. Must be a numeric vector of 2 numbers. Default
#' \code{c(10, 1)}.
#' @param mainNRow,mainNCol Arrangement of the main plotting region, for number
#' of rows and columns. Default \code{NULL} will be automatically handled by
#' \code{\link[cowplot]{plot_grid}}.
#' @param legendNRow Arrangement of the legend, number of rows. Default
#' \code{1}.
#' @inheritDotParams .ggScatter shapeBy dotOrder dotSize dotAlpha raster labelText labelTextSize seed
#' @inheritDotParams .ggplotLigerTheme subtitle baseSize titleSize subtitleSize xTextSize xTitleSize yTextSize yTitleSize panelBorder colorValues naColor plotly
#' @return ggplot object when only one feature (e.g. cluster variable, gene,
#' factor) is set. List object when multiple of those are specified.
#' @export
#' @examples
#' plotGroupClusterDimRed(pbmcPlot)
plotGroupClusterDimRed <- function(
        object,
        useGroup = "dataset",
        useCluster = NULL,
        useDimRed = NULL,
        combinePlot = TRUE,
        droplevels = TRUE,
        relHeightMainLegend = c(5, 1),
        relHeightDRBar = c(10, 1),
        mainNRow = NULL,
        mainNCol = NULL,
        legendNRow = 1,
        ...
) {
    useCluster <- useCluster %||% object@uns$defaultCluster
    clusterVar <- cellMeta(object, useCluster)
    groupVar <- .fetchCellMetaVar(object, useGroup, checkCategorical = TRUE,
                                  droplevels = droplevels)
    plotList <- list()
    propPList <- list()
    for (lvl in levels(groupVar)) {
        mask <- groupVar == lvl
        clusterVarMasked <- clusterVar
        clusterVarMasked[!mask] <- NA
        clusterVarMasked <- droplevels(clusterVarMasked)
        tempVarName <- paste0(useCluster, "_", lvl)
        cellMeta(object, tempVarName) <- clusterVarMasked
        plotList[[lvl]] <- plotDimRed(
            object, colorBy = tempVarName, useDimRed = useDimRed,
            slot = "cellMeta", dotOrder = "shuffle", titles = lvl,
            legendColorTitle = "", legendPosition = "bottom",
            legendNRow = legendNRow, ...
        ) +
            ggplot2::theme(
                line = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank()
            )
        proportions <- table(clusterVarMasked) / sum(table(clusterVarMasked))
        propDF <- data.frame(group = lvl, proportions)
        # propDF$clusterVarMasked <- as.character()
        # Reverse row order so the bar plot shows in the same order as legend
        propDF <- propDF[rev(seq_len(nrow(propDF))),]
        propPList[[lvl]] <- ggplot2::ggplot(propDF, ggplot2::aes(y = .data[["group"]],
                                                                 x = .data[["Freq"]],
                                                                 fill = .data[["clusterVarMasked"]])) +
            ggplot2::geom_bar(stat = "identity")
        propPList[[lvl]] <- .ggplotLigerTheme(propPList[[lvl]], ...) +
            ggplot2::theme_void() +
            ggplot2::theme(legend.position = "none")
    }
    if (!isTRUE(combinePlot)) return(plotList)
    suppressWarnings({
        legend <- cowplot::get_legend(plotList[[1]])
        plotList <- lapply(plotList, function(gg) {
            return(gg + ggplot2::theme(legend.position = "none"))
        })
        plotList <- lapply(names(plotList), function(lvl) {
            cowplot::plot_grid(plotList[[lvl]], propPList[[lvl]], nrow = 2,
                               rel_heights = relHeightDRBar)
        })
        main <- cowplot::plot_grid(plotlist = plotList, ncol = mainNCol,
                                   nrow = mainNRow)
        final <- cowplot::plot_grid(main, legend, nrow = 2, rel_heights = relHeightMainLegend)
    })
    return(final)
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
#' @inheritDotParams plotCellViolin slot yFunc cellIdx titles
#' @inheritDotParams .ggCellViolin violin violinAlpha violinWidth box boxAlpha boxWidth dot dotColor dotSize xlabAngle raster seed
#' @inheritDotParams .ggplotLigerTheme subtitle xlab ylab legendFillTitle showLegend legendPosition baseSize titleSize subtitleSize xTextSize xTitleSize yTextSize yTitleSize legendTextSize legendTitleSize legendNRow legendNCol colorLabels colorValues panelBorder plotly
#' @param ... Additional arguments passed to \code{\link{plotCellViolin}}.
#' @return ggplot if using a single gene and not splitting by dataset.
#' Otherwise, list of ggplot.
#' @export
#' @rdname plotViolin
#' @examples
#' plotGeneViolin(pbmcPlot, varFeatures(pbmcPlot)[1],
#'                groupBy = "leiden_cluster")
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

    groupBy <- groupBy %||% object@uns$defaultCluster
    if (isFALSE(groupBy)) groupBy <- NULL

    plotList <- plotCellViolin(
        object,
        y = gene,
        slot = "normData",
        yFunc = function(x) log2(10000*x + 1),
        groupBy = groupBy,
        splitBy = splitBy,
        ...
    )

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


#' Create violin plot for multiple genes grouped by clusters
#' @description
#' Make violin plots for each given gene grouped by cluster variable and stack
#' along y axis.
#'
#' @details
#' If \code{xlab} need to be set, set \code{xlabAngle} at the same time. This is
#' due to that the argument parsing mechanism will partially match it to main
#' function arguments before matching the \code{...} arguments.
#' @param object A \linkS4class{liger} object.
#' @param gene Character vector of gene names.
#' @param groupBy The name of an available categorical variable in
#' \code{cellMeta} slot. This forms the main x-axis columns. Use \code{FALSE}
#' for no grouping. Default \code{NULL} looks clustering result but will not
#' group if no clustering is found.
#' @param colorBy The name of another categorical variable in \code{cellMeta}
#' slot. This split the main grouping columns and color the violins. Default
#' \code{NULL} will not split and color the violins.
#' @param box Logical, whether to add boxplot. Default \code{FALSE}.
#' @param boxAlpha Numeric, transparency of boxplot. Default \code{0.1}.
#' @param yFunc Function to transform the y-axis. Default is
#' \code{log1p(x*1e4)}. Set to \code{NULL} for no transformation.
#' @param showLegend Whether to show the legend. Default \code{FALSE}.
#' @param xlabAngle Numeric, counter-clockwise rotation angle in degrees of X
#' axis label text. Default \code{40}.
#' @inheritDotParams .ggplotLigerTheme title subtitle xlab ylab legendFillTitle legendPosition baseSize titleSize xTitleSize yTitleSize legendTitleSize subtitleSize xTextSize yTextSize legendTextSize yFacetSize panelBorder legendNRow legendNCol colorLabels colorValues plotly
#' @return A ggplot object.
#' @export
#' @examples
#' plotClusterGeneViolin(pbmcPlot, varFeatures(pbmcPlot)[1:10])
plotClusterGeneViolin <- function(
        object,
        gene,
        groupBy = NULL,
        colorBy = NULL,
        box = FALSE,
        boxAlpha = 0.1,
        yFunc = function(x) log1p(x*1e4),
        showLegend = !is.null(colorBy),
        xlabAngle = 40,
        ...
) {
    groupBy <-  groupBy %||% object@uns$defaultCluster
    gene <- unique(gene)
    featureDF <- retrieveCellFeature(object, gene, slot = "normData")
    # Account for the case where the gene is not detected in any cell
    ngene <- ncol(featureDF)
    geneUse <- gene[gene %in% colnames(featureDF)]
    if (!is.null(yFunc)) featureDF <- as.data.frame(apply(featureDF, 2, yFunc))
    if (!isFALSE(groupBy)) {
        featureDF[[groupBy]] <- .fetchCellMetaVar(
            object = object,
            variables = groupBy,
            checkCategorical = TRUE
        )
    } else {
        groupBy <- "group"
        featureDF[[groupBy]] <- factor("All")
    }
    if (!is.null(colorBy)) {
        colorVar <- .fetchCellMetaVar(
            object = object, variables = colorBy, checkCategorical = TRUE
        )
        featureDF[[colorBy]] <- factor(colorVar)
    }

    featureDF <- featureDF %>%
        .pivot_longer(
            cols = seq(ngene),
            names_to = "gene",
            values_to = "Expression"
        ) %>%
        dplyr::mutate(gene = factor(.data[["gene"]], levels = geneUse))
    if (!is.null(colorBy)) {
        p <- ggplot2::ggplot(
            data = featureDF,
            mapping = ggplot2::aes(
                x = .data[[groupBy]],
                y = .data[["Expression"]],
                fill = .data[[colorBy]]
            ))
    } else {
        p <- ggplot2::ggplot(
            data = featureDF,
            mapping = ggplot2::aes(
                x = .data[[groupBy]],
                y = .data[["Expression"]],
                fill = .data[[groupBy]]
            ))
    }
    p <- p + ggplot2::geom_violin()
    if (isTRUE(box)) {
        p <- p + ggplot2::geom_boxplot(alpha = boxAlpha)
    }
    p <- p +
        ggplot2::facet_wrap(
            stats::formula("~gene"),
            scales = "free_y",
            nrow = ngene,
            strip.position = "left"
        )
    .ggplotLigerTheme(
        p,
        showLegend = showLegend,
        xlabAngle = xlabAngle,
        ...
    )
}




#' Create barcode-rank plot for each dataset
#' @description
#' This function ranks the total count of each cell within each dataset and make
#' line plot. This function is simply for examining the input raw count data
#' and does not infer any recommended cutoff for removing non-cell barcodes.
#' @param object A \linkS4class{liger} object.
#' @inheritDotParams .ggScatter dotSize dotAlpha raster
#' @inheritDotParams .ggplotLigerTheme title subtitle xlab ylab baseSize titleSize subtitleSize xTextSize xTitleSize yTextSize yTitleSize panelBorder plotly
#' @export
#' @return A list object of ggplot for each dataset
#' @examples
#' plotBarcodeRank(pbmc)
plotBarcodeRank <- function(
        object,
        ...
) {
    pl <- list()
    for (d in names(object)) {
        df <- data.frame(
            rank = seq_len(lengths(object)[d]),
            nUMI = sort(
                cellMeta(object, columns = "nUMI", useDatasets = d),
                decreasing = TRUE
            )
        )
        pl[[d]] <- .ggScatter(df, x = "rank", y = "nUMI", ...) +
            ggplot2::scale_y_log10() +
            ggplot2::scale_x_log10() +
            ggplot2::geom_line()
    }
    pl
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
#' annotation on the pie chart. Run \code{options(ggrepel.max.overlaps = n)}
#' before plotting to set allowed label overlaps.
#' @param object A \linkS4class{liger} object.
#' @param class1,class2 Each should be a single name of a categorical variable
#' available in \code{cellMeta} slot. Number of cells in each categories in
#' \code{class2} will be served as the denominator when calculating proportions.
#' By default \code{class1 = NULL} and uses default clusters and \code{class2 =
#' "dataset"}.
#' @param method For bar plot, choose whether to draw \code{"stack"} or
#' \code{"group"} bar plot. Default \code{"stack"}.
#' @param showLegend Whether to show the legend. Default \code{TRUE}.
#' @param panelBorder Whether to show rectangle border of the panel instead of
#' using ggplot classic bottom and left axis lines. Default \code{FALSE}.
#' @param inclRev Logical, for barplot, whether to reverse the specification for
#' \code{class1} and \code{class2} and produce two plots. Default \code{FALSE}.
#' @param combinePlot Logical, whether to combine the two plots with
#' \code{\link[cowplot]{plot_grid}} when two plots are created. Default
#' \code{TRUE}.
#' @param useCluster For \code{plotClusterProportions}. Same as \code{class1}
#' while \code{class2} is hardcoded with \code{"dataset"}.
#' @param labelSize,labelColor Settings on pie chart percentage label. Default
#' \code{4} and \code{"white"}.
#' @inheritDotParams .ggplotLigerTheme title subtitle xlab ylab legendFillTitle showLegend legendPosition baseSize titleSize subtitleSize xTextSize xTitleSize yTextSize yTitleSize legendTextSize legendTitleSize panelBorder legendNRow legendNCol colorLabels colorValues colorPalette colorDirection naColor colorLow colorMid colorHigh colorMidPoint plotly
#' @param return.plot `r lifecycle::badge("defunct")`
#' @return ggplot or list of ggplot
#' @rdname plotProportion
#' @export
#' @examples
#' plotProportion(pbmcPlot)
#' plotProportionBar(pbmcPlot, method = "group")
#' plotProportionPie(pbmcPlot)
plotProportion <- function(
        object,
        class1 = NULL,
        class2 = "dataset",
        method = c("stack", "group", "pie"),
        ...
) {
    class1 <- class1 %||% object@uns$defaultCluster
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
        class1 = NULL,
        class2 = "dataset",
        showLegend = FALSE,
        panelBorder = TRUE,
        ...
) {
    class1 <- class1 %||% object@uns$defaultCluster
    if (length(class1) != 1 ||
        length(class2) != 1)
        cli::cli_abort(
            "{.var class1} and {.var class2} must be name of one categorical variable in {.code cellMeta(object)}"
        )
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
        class1 = NULL,
        class2 = "dataset",
        method = c("stack", "group"),
        inclRev = FALSE,
        panelBorder = TRUE,
        combinePlot = TRUE,
        ...
) {
    class1 <- class1 %||% object@uns$defaultCluster
    if (length(class1) != 1 ||
        length(class2) != 1)
        cli::cli_abort("{.var class1} and {.var class2} must be name of one categorical variable in {.code cellMeta(object)}")
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
        useCluster = NULL,
        return.plot = FALSE,
        ...
) {
    .deprecateArgs(defunct = "return.plot")
    lifecycle::deprecate_warn(
        when = "1.99.0", what = "plotClusterProportions()",
        with = "plotProportionDot()",
        details = "See help(\"plotProportion\") for more new options."
    )
    useCluster <- useCluster %||% object@uns$defaultCluster
    plotProportionDot(object, class1 = useCluster, class2 = "dataset", ...)
}

#' @rdname plotProportion
#' @param circleColors Character vector of colors. \code{plotProportionPie}
#' parameter for setting the colors of circles, i.e. categorical variable
#' controlled by \code{class2}. Default \code{NULL} uses ggplot default hues.
#' @export
plotProportionPie <- function(
        object,
        class1 = NULL,
        class2 = "dataset",
        labelSize = 4,
        labelColor = "black",
        circleColors = NULL,
        ...
) {
    class1 <- class1 %||% object@uns$defaultCluster
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
        ggplot2::geom_tile(colour = "white", linewidth = 0.3, height = 0.9) +
        ggplot2::coord_polar()
    # Add class2 annotation
    class2fct <- droplevels(df[[class2]])
    class2Uniq <- levels(class2fct)
    lgdDF <- data.frame(slope = 0, intercept = seq_along(class2Uniq) + 0.5,
                        class2 = class2Uniq)
    names(lgdDF)[3] <- class2
    p <- p +
        ggplot2::geom_abline(
            data = lgdDF,
            mapping = ggplot2::aes(slope = .data[["slope"]],
                                   intercept = .data[["intercept"]],
                                   colour = .data[[class2]]),
            linewidth = 2
        )
    if (!is.null(circleColors)) {
        if (length(circleColors) < length(class2Uniq)) {
            cli::cli_alert_warning("Less {.field circleColors} ({.val {length(circleColors)}}) specified than required from {.val {class2}} ({.val {length(class2Uniq)}}). Using default.")
        } else {
            p <- p + ggplot2::scale_color_manual(values = circleColors)
        }
    }

    if (!requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggplot2::geom_text(
            size = labelSize, color = labelColor,
            position = ggplot2::position_nudge(y = 0.25)
        )
    } else {
        p <- p + ggrepel::geom_text_repel(
            size = labelSize, color = labelColor, force = 1,
            nudge_y = 0.25, bg.colour = "white"
        )
    }
    .ggplotLigerTheme(p, ...) +
        ggplot2::theme(
            axis.line = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank()
        )
}

#' Box plot of cluster proportion in each dataset, grouped by condition
#' @description
#' This function calculate the proportion of each category (e.g. cluster, cell
#' type) within each dataset, and then make box plot grouped by condition. The
#' proportion of all categories within one dataset sums up to 1. The condition
#' variable must be a variable of dataset, i.e. each dataset must belong to only
#' one condition.
#'
#' @param object A \linkS4class{liger} object.
#' @param useCluster Name of variable in \code{cellMeta(object)}. Default
#' \code{NULL} uses default cluster.
#' @param conditionBy Name of the variable in \code{cellMeta(object)} that
#' represents the condition. Must be a high level variable of the
#' \code{sampleBy} variable, i.e. each sample must belong to only one condition.
#' Default \code{NULL} does not group samples by condition.
#' @param sampleBy Name of the variable in \code{cellMeta(object)} that
#' represents individual samples. Default \code{"dataset"}.
#' @param splitByCluster Logical, whether to split the wide grouped box plot by
#' cluster, into a list of boxplots for each cluster. Default \code{FALSE}.
#' @param dot Logical, whether to add dot plot on top of the box plot. Default
#' \code{FALSE}.
#' @param dotSize Size of the dot. Default uses user option "ligerDotSize", or
#' \code{1} if not set.
#' @param dotJitter Logical, whether to jitter the dot to avoid overlapping
#' within a box when many dots are presented. Default \code{FALSE}.
#' @inheritDotParams .ggplotLigerTheme title subtitle xlab ylab legendFillTitle showLegend legendPosition baseSize titleSize subtitleSize xTextSize xTitleSize yTextSize yTitleSize legendTextSize legendTitleSize panelBorder legendNRow legendNCol colorLabels colorValues colorPalette colorDirection naColor colorLow colorMid colorHigh colorMidPoint plotly
#' @export
#' @return A ggplot object or a list of ggplot objects if
#' \code{splitByCluster = TRUE}.
#' @examples
#' # "boxes" are expected to appear as horizontal lines, because there's no
#' # "condition" variable that groups the datasets in the example object, and
#' # thus only one value exists for each "box".
#' plotProportionBox(pbmcPlot, conditionBy = "dataset")
plotProportionBox <- function(
        object,
        useCluster = NULL,
        conditionBy = NULL,
        sampleBy = "dataset",
        splitByCluster = FALSE,
        dot = FALSE,
        dotSize = getOption("ligerDotSize", 1),
        dotJitter = FALSE,
        ...
) {
    useCluster <- useCluster %||% object@uns$defaultCluster
    if (is.null(useCluster)) {
        cli::cli_abort("No cluster specified nor default set.")
    }
    clusterVar <- .fetchCellMetaVar(object, useCluster, checkCategorical = TRUE)
    datasetVar <- .fetchCellMetaVar(object, sampleBy, checkCategorical = TRUE)
    compositionTable <- table(datasetVar, clusterVar)
    dfLong <- data.frame(compositionTable)
    names(dfLong) <- c(sampleBy, useCluster, "Count")

    if (!is.null(conditionBy)) {
        conditionVar <- .fetchCellMetaVar(
            object = object, variables = conditionBy, checkCategorical = TRUE
        )
        # Check that condition variable is strictly a high level variable of dataset
        if (!all(rowSums(table(datasetVar, conditionVar) > 0) == 1)) {
            cli::cli_abort("Condition variable must be a high level variable of the datasets, i.e. each dataset must belong to only one condition.")
        }

        conditionTable <- table(datasetVar, conditionVar)
        conditionMap <- apply(
            conditionTable,
            MARGIN = 1,
            function(row) colnames(conditionTable)[row > 0]
        )

        dfLong[[conditionBy]] <- factor(
            conditionMap[dfLong[[sampleBy]]],
            levels = levels(conditionVar)
        )
    }

    dfLong %<>%
        dplyr::group_by(dataset) %>%
        dplyr::mutate(
            Proportion = .data[["Count"]] / sum(.data[["Count"]]),
        )
    if (isTRUE(splitByCluster)) {
        plist <- lapply(levels(clusterVar), function(cluster) {
            p <- ggplot2::ggplot(
                data = dfLong[dfLong[[useCluster]] == cluster, ],
                mapping = (
                    if (!is.null(conditionBy))
                        ggplot2::aes(
                            x = .data[[conditionBy]],
                            y = .data[["Proportion"]],
                            fill = .data[[conditionBy]]
                        )
                    else ggplot2::aes(
                        x = .data[[useCluster]],
                        y = .data[["Proportion"]]
                    )
                )
            ) + (
                if (isTRUE(dot))
                    ggplot2::geom_point(
                        size = dotSize,
                        color = "black",
                        position =
                            if (isTRUE(dotJitter)) ggplot2::position_jitter()
                        else "identity"
                    )
                else
                    NULL
            ) +
                ggplot2::geom_boxplot() +
                ggplot2::ggtitle(paste0(useCluster, ": ", cluster))
            return(.ggplotLigerTheme(p, ...))
        })
        names(plist) <- levels(clusterVar)
        return(plist)
    } else {
        p <- ggplot2::ggplot(
            data = dfLong,
            mapping = (
                if (!is.null(conditionBy))
                    ggplot2::aes(
                        x = .data[[useCluster]],
                        y = .data[["Proportion"]],
                        fill = .data[[conditionBy]]
                    )
                else ggplot2::aes(
                    x = .data[[useCluster]],
                    y = .data[["Proportion"]]
                )
            )
        ) +
            (if (isTRUE(dot))
                ggplot2::geom_point(
                    size = dotSize,
                    color = "black",
                    position =
                        if (isTRUE(dotJitter)) ggplot2::position_jitter()
                    else "identity"
                )
             else
                 NULL) +
            ggplot2::geom_boxplot()
        return(.ggplotLigerTheme(p, ...))
    }

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Volcano plot ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create volcano plot for Wilcoxon test result
#' @description
#' \code{plotVolcano} is a simple implementation and shares most of arguments
#' with other rliger plotting functions. \code{plotEnhancedVolcano} is a
#' wrapper function of \code{EnhancedVolcano::EnhancedVolcano()}, which has
#' provides substantial amount of arguments for graphical control. However, that
#' requires the installation of package "EnhancedVolcano".
#'
#' \code{highlight} and \code{labelTopN} both controls the feature name
#' labeling, whereas \code{highlight} is considered first. If both are as
#' default (\code{NULL}), all significant features will be labeled.
#' @param result Data frame table returned by \code{\link{runMarkerDEG}} or
#' \code{\link{runPairwiseDEG}}.
#' @param group Selection of one group available from \code{result$group}. If
#' only one group is available from \code{result}, default \code{NULL} uses it.
#' @param logFCThresh Number for the threshold on the absolute value of the log2
#' fold change statistics. Default \code{1}.
#' @param padjThresh Number for the threshold on the adjusted p-value
#' statistics. Default \code{0.01}.
#' @param highlight A character vector of feature names to be highlighted.
#' Default \code{NULL}.
#' @param labelTopN Number of top differential expressed features to be labeled
#' on the top of the dots. Ranked by adjusted p-value first and absolute value
#' of logFC next. Default \code{NULL}.
#' @param dotSize,dotAlpha Numbers for universal aesthetics control of dots.
#' Default \code{2} and \code{0.8}.
#' @param legendPosition Text indicating where to place the legend. Choose from
#' \code{"top"}, \code{"bottom"}, \code{"left"} or \code{"right"}. Default
#' \code{"top"}.
#' @param labelSize Size of labeled top features and line annotations. Default
#' \code{4}.
#' @inheritDotParams .ggScatter dotOrder raster labelText labelTextSize seed
#' @inheritDotParams .ggplotLigerTheme title subtitle legendColorTitle showLegend baseSize titleSize subtitleSize xTextSize xTitleSize yTextSize yTitleSize legendTextSize legendTitleSize panelBorder
#' @return ggplot
#' @export
#' @examples
#' plotVolcano(deg.pw, "0.stim")
plotVolcano <- function(
        result,
        group = NULL,
        logFCThresh = 1,
        padjThresh = 0.01,
        highlight = NULL,
        labelTopN = NULL,
        dotSize = 2,
        dotAlpha = 0.8,
        legendPosition = "top",
        labelSize = 4,
        ...
) {
    if (!is.factor(result$group)) result$group <- factor(result$group)
    result$group <- droplevels(result$group)
    if (is.null(group)) {
        if (nlevels(result$group) == 1) groupUse <- levels(result$group)
        else {
            cli::cli_abort(
                c("Please specify one group to visualize",
                  i = "Available ones: {.val {levels(result$group)}}")
            )
        }
    } else if (length(group) != 1 ||
        !group %in% result$group) {
        cli::cli_abort(
            c("Please specify one available group to visualize",
              i = "Available ones: {.val {unique(result$group)}}")
        )
    }
    # Avoid dplyr namespace weird bug
    # That when "group" is a column in the data, and `group` exists as a
    # variable out of the pipe, doing `.data[['group']] == group` will be
    # interpreted as `.data[['group']] == .data[['group']]`
    groupUse <- group
    minPosPadj <- min(result$padj[result$padj > 0], na.rm = TRUE) / 10

    result <- result %>%
        dplyr::filter(.data[['group']] == groupUse, !is.na(.data[['padj']])) %>%
        dplyr::mutate(Significance = dplyr::case_when(
            abs(.data[['logFC']]) > logFCThresh &
                .data[['padj']] < padjThresh ~ "padj & logFC",
            abs(.data[['logFC']]) > logFCThresh ~ "logFC",
            .data[['padj']] < padjThresh ~ "padj",
            .default = "Not significant"
        )) %>%
        dplyr::mutate(Significance = factor(.data[["Significance"]],
                                            levels = c("Not significant",
                                                       "logFC", "padj", "padj & logFC"))) %>%
        dplyr::mutate(padj = ifelse(.data[['padj']] == 0, minPosPadj, .data[['padj']])) %>%
        dplyr::mutate(padj = -log10(.data[['padj']])) %>%
        dplyr::arrange(dplyr::desc(.data[['padj']]),
                       dplyr::desc(.data[['logFC']]))

    # Decide which genes to label
    if (!is.null(highlight)) {
        notFound <- highlight[!highlight %in% result$feature]
        if (length(notFound) > 0) {
            cli::cli_alert_warning(
                "The following features are not found in the result or got NA p-value: {.val {notFound}}"
            )
        }
        result <- result %>%
            dplyr::mutate(
                label = dplyr::case_when(
                    .data[['feature']] %in% highlight ~ .data[['feature']],
                    TRUE ~ NA_character_
                )
            )
    } else {
        if (is.null(labelTopN)) {
            labelTopN <- sum(result$Significance == "padj & logFC")
        } else {
            labelTopN <- min(labelTopN, sum(result$Significance == "padj & logFC"))
        }
        labelIdx <- which(result$Significance == "padj & logFC")[seq(labelTopN)]
        result$label <- NA_character_
        result$label[labelIdx] <- result$feature[labelIdx]
    }

    # Prepare for lines that mark the cutoffs
    vlineLab <- data.frame(
        X = c(-logFCThresh, logFCThresh)
    )
    hlineLab <- data.frame(
        Y = c(-log10(padjThresh))
    )
    halfWidth <- max(abs(result$logFC), na.rm = TRUE)
    height <- max(result$padj, na.rm = TRUE)
    lfcVertDrift <- height * 0.03
    lfcHorizDrift <- halfWidth * 0.02
    padjVertDrift <- height * 0.02
    padjHorizDrift <- halfWidth * 0.001
    p <- .ggScatter(result, x = "logFC", y = "padj",
                    colorBy = "Significance", zeroAsNA = FALSE,
                    labelBy = "label",
                    xlab = "Log2 Fold Change",
                    ylab = "-log10 Adjusted P-value",
                    colorValues = c("black", "#ef2301", "#416ae1", "#238b22"),
                    legendPosition = legendPosition,
                    dotSize = dotSize, dotAlpha = dotAlpha,
                    ggrepelLabelTick = TRUE, ...) +
        ggplot2::xlim(-max(abs(result$logFC)), max(abs(result$logFC))) +
        ggplot2::geom_vline(data = vlineLab,
                            mapping = ggplot2::aes(xintercept = .data[["X"]]),
                            linetype = "longdash") +
        ggplot2::geom_hline(data = hlineLab,
                            ggplot2::aes(yintercept = .data[["Y"]]),
                            linetype = "longdash") +
        ggplot2::annotate("text",
                          x = -logFCThresh - lfcHorizDrift, y = -lfcVertDrift,
                          label = paste0("lower log2FC cutoff: ", -logFCThresh),
                          size = labelSize, hjust = 1, vjust = 1) +
        ggplot2::annotate("text",
                          x = logFCThresh + lfcHorizDrift, y = -lfcVertDrift,
                          label = paste0("higher log2FC cutoff: ", logFCThresh),
                          size = labelSize, hjust = 0, vjust = 1) +
        ggplot2::annotate("text", angle = 90,
                          x = - halfWidth + padjHorizDrift,
                          y = -log10(padjThresh) + padjVertDrift,
                          hjust = 0, vjust = 1,
                          label = paste("p-adj cutoff:", padjThresh),
                          size = labelSize)
    return(p)
}

#' Create volcano plot with EnhancedVolcano
#' @inheritParams plotVolcano
#' @param ... Arguments passed to EnhancedVolcano::EnhancedVolcano(), except
#' that \code{toptable}, \code{lab}, \code{x} and \code{y} are prefilled by this
#' wrapper.
#' @returns ggplot
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
#'     defaultCluster(pbmc) <- pbmcPlot$leiden_cluster
#'     # Test the DEG between "stim" and "ctrl", within each cluster
#'     result <- runPairwiseDEG(
#'         pbmc,
#'         groupTest = "stim",
#'         groupCtrl = "ctrl",
#'         variable1 = "dataset",
#'         splitBy = "defaultCluster"
#'     )
#'     plotEnhancedVolcano(result, "0.stim")
#' }
#' }
plotEnhancedVolcano <- function(
        result,
        group,
        ...
) {
    if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) { # nocov start
        cli::cli_abort(
            "Package {.pkg EnhancedVolcano} is needed for this function to work.
            Please install it by command:
            {.code BiocManager::install('EnhancedVolcano')}")
    } # nocov end
    result <- result[result$group == group, ]
    EnhancedVolcano::EnhancedVolcano(
        toptable = result,
        lab = result$feature, x = "logFC", y = "padj",
        ...
    )
}



#' Create density plot basing on specified coordinates
#' @description This function shows the cell density presented in a 2D
#' dimensionality reduction coordinates. Density is shown with coloring and
#' contour lines. A scatter plot of the dimensionality reduction is added as
#' well. The density plot can be splitted by categorical variables (e.g.
#' \code{"dataset"}), while the scatter plot will always be shown for all cells
#' in subplots as a reference of the global structure.
#' @param object A \linkS4class{liger} object
#' @param useDimRed Name of the variable storing dimensionality reduction result
#' in the \code{cellMeta} slot. Default uses default dimension reduction.
#' @param splitBy Character vector of categorical variable names in
#' \code{cellMeta} slot. Split all cells by groupings on this/these variable(s)
#' to produce a density plot containing only the cells in each group. Default
#' \code{NULL}.
#' @param combinePlot Logical, whether to utilize
#' \code{\link[cowplot]{plot_grid}} to combine multiple plots into one. Default
#' \code{TRUE} returns combined ggplot. \code{FALSE} returns a list of ggplot
#' or a single ggplot when only one plot is requested.
#' @param minDensity A positive number to filter out low density region colored
#' on plot. Default \code{8}. Setting zero will show density on the whole panel.
#' @param contour Logical, whether to draw the contour line. Default
#' \code{TRUE}.
#' @param contourLineWidth Numeric, the width of the contour line. Default
#' \code{0.3}.
#' @param contourBins Number of contour bins. Higher value generates more
#' contour lines. Default \code{5}.
#' @param dot Logical, whether to add scatter plot of all cells, even when
#' density plot is splitted with \code{splitBy}. Default \code{TRUE}.
#' @param dotColor,dotSize,dotAlpha Numeric, controls the appearance of all
#' dots. Default \code{"grey"}, \code{0.6} and \code{0.3}, respectively.
#' @param dotRaster Logical, whether to rasterize the scatter plot. Default
#' \code{NULL} automatically rasterizes the dots when number of total cells to
#' be plotted exceeds 100,000.
#' @param title Text of main title of the plots. Default \code{NULL}. Length
#' of character vector input should match with number of plots generated.
#' @param legendFillTitle Text of legend title. Default \code{"Density"}.
#' @param colorPalette Name of the option for
#' \code{\link[ggplot2]{scale_fill_viridis_c}}. Default \code{"magma"}.
#' @param colorDirection Color gradient direction for
#' \code{\link[ggplot2]{scale_fill_viridis_c}}. Default \code{-1}.
#' @inheritDotParams .ggplotLigerTheme title subtitle xlab ylab showLegend legendPosition baseSize titleSize xTitleSize yTitleSize legendTitleSize subtitleSize xTextSize yTextSize legendTextSize panelBorder plotly
#' @return A ggplot object when only one plot is generated, A ggplot object
#' combined with \code{\link[cowplot]{plot_grid}} when multiple plots and
#' \code{combinePlot = TRUE}. A list of ggplot when multiple plots and
#' \code{combinePlot = FALSE}.
#' @export
#' @examples
#' # Example dataset has small number of cells, thus cutoff adjusted.
#' plotDensityDimRed(pbmcPlot, minDensity = 1)
plotDensityDimRed <- function(
        object,
        useDimRed = NULL,
        splitBy = NULL,
        combinePlot = TRUE,
        minDensity = 8,
        contour = TRUE,
        contourLineWidth = 0.3,
        contourBins = 5,
        dot = TRUE,
        dotColor = "grey",
        dotSize = 0.6,
        dotAlpha = 0.3,
        dotRaster = NULL,
        title = NULL,
        legendFillTitle = "Density",
        colorPalette = "magma",
        colorDirection = -1,
        ...
) {
    dr <- as.data.frame(dimRed(object, useDimRed))
    splitVar <- .fetchCellMetaVar(object, splitBy, checkCategorical = TRUE,
                                  drop = FALSE)

    if (!is.null(splitVar) && ncol(splitVar) > 0) {
        # Will be splitting into subplots
        drList <- split(dr, splitVar)
        title <- .checkArgLen(title, length(drList), .stop = FALSE)
    } else {
        # Will return a single ggplot
        if (length(title) > 1) {
            cli::cli_alert_warning("{.var title} has length greater than 1 while only a single plot is generated. Using the first value only.")
            title <- title[1]
        }
        drList <- list(dr)
    }
    plotList <- list()
    if (length(drList) == 1) {
        return(.ggDensity(drList[[1]], dotCoordDF = drList[[1]],
                          title = title, minDensity = minDensity,
                          contour = contour,
                          contourLineWidth = contourLineWidth,
                          contourBins = contourBins, dot = dot,
                          dotColor = dotColor, dotSize = dotSize,
                          dotAlpha = dotAlpha,
                          legendFillTitle = legendFillTitle,
                          colorPalette = colorPalette,
                          colorDirection = colorDirection,
                          dotRaster = dotRaster, ...))
    } else {
        if (is.null(title)) title <- names(drList)
        for (i in seq_along(drList))
            plotList[[i]] <- .ggDensity(drList[[i]], dotCoordDF = dr,
                                        title = title[i],
                                        minDensity = minDensity,
                                        contour = contour,
                                        contourLineWidth = contourLineWidth,
                                        contourBins = contourBins,
                                        dotColor = dotColor, dot = dot,
                                        dotSize = dotSize, dotAlpha = dotAlpha,
                                        legendFillTitle = legendFillTitle,
                                        colorPalette = colorPalette,
                                        colorDirection = colorDirection,
                                        dotRaster = dotRaster,
                                        ...)
        names(plotList) <- names(drList)
        if (isTRUE(combinePlot))
            return(cowplot::plot_grid(plotlist = plotList,
                                      align = "hv", axis = "tblr"))
        else return(plotList)
    }
}

.ggDensity <- function(
        coordDF,
        dotCoordDF,
        minDensity = 8,
        contour = TRUE,
        contourLineWidth = 0.3,
        contourBins = 5,
        dot = TRUE,
        dotColor = "grey",
        dotSize = 0.6,
        dotAlpha = 0.3,
        dotRaster = NULL,
        legendFillTitle = "Density",
        colorPalette = "magma",
        colorDirection = -1,
        ...
) {
    dotRaster <- .checkRaster(nrow(dotCoordDF), dotRaster)
    zeroAsNA <- function(x) {
        x[x < minDensity] <- NA
        x
    }
    x <- colnames(coordDF)[1]
    y <- colnames(coordDF)[2]
    xLim <- c(min(coordDF[,1]) - 2, max(coordDF[,1]) + 2)
    yLim <- c(min(coordDF[,2]) - 2, max(coordDF[,2]) + 2)
    p <- ggplot2::ggplot(coordDF,
                         ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
        ggplot2::stat_density_2d_filled(
            ggplot2::aes(
                fill = zeroAsNA(ggplot2::after_stat(.data[["count"]]))
            ),
            geom = "raster", contour = FALSE, na.rm = TRUE) +
        ggplot2::scale_fill_viridis_c(
            option = colorPalette,
            direction = colorDirection,
            na.value = "white",
            guide = ggplot2::guide_colorbar(title = legendFillTitle)) +
        ggplot2::xlim(xLim) +
        ggplot2::ylim(yLim)
    if (isTRUE(contour)) {
        p <- p +
            ggplot2::geom_density_2d(color = "black",
                                     linewidth = contourLineWidth,
                                     bins = contourBins)
    }
    if (isTRUE(dot)) {
        if (!isTRUE(dotRaster)) {
            p <- p +
                ggplot2::geom_point(
                    data = dotCoordDF, size = dotSize, stroke = 0,
                    colour = "grey", alpha = dotAlpha
                )
        } else {
            p <- p +
                scattermore::geom_scattermore(
                    data = dotCoordDF, pointsize = dotSize, stroke = 0,
                    color = "grey", alpha = dotAlpha
                )
        }
    }

    .ggplotLigerTheme(p, ...)
}


#' Visualize factor expression and gene loading
#' @rdname plotGeneLoadings
#' @param object A \linkS4class{liger} object with valid factorization result.
#' @param markerTable Returned result of \code{\link{getFactorMarkers}}.
#' @param useFactor Integer index for which factor to visualize.
#' @param useDimRed Name of the variable storing dimensionality reduction result
#' in the \code{cellMeta} slot. Default \code{"UMAP"}.
#' @param nLabel Integer, number of top genes to be shown with text labels.
#' Default \code{15}.
#' @param nPlot Integer, number of top genes to be shown in the loading rank
#' plot. Default \code{30}.
#' @inheritDotParams plotDimRed colorByFunc cellIdx shapeBy titles
#' @inheritDotParams .ggScatter dotSize dotAlpha trimHigh trimLow raster
#' @inheritDotParams .ggplotLigerTheme xlab ylab legendColorTitle legendShapeTitle showLegend legendPosition baseSize titleSize subtitleSize xTextSize xTitleSize yTextSize yTitleSize legendTextSize legendTitleSize colorPalette colorDirection naColor panelBorder
#' @export
#' @examples
#' result <- getFactorMarkers(pbmcPlot, "ctrl", "stim")
#' plotGeneLoadings(pbmcPlot, result, useFactor = 2)
plotGeneLoadings <- function(
        object,
        markerTable,
        useFactor,
        useDimRed = NULL,
        nLabel = 15,
        nPlot = 30,
        ...
) {
    allDotArgs <- list(...)
    pdrArgs <- names(as.list(args(plotDimRed)))
    pglrArgs <- names(as.list(args(plotGeneLoadingRank)))
    pdrDotArgs <- intersect(names(allDotArgs), pdrArgs)
    pglrDotArgs <- intersect(names(allDotArgs), pglrArgs)
    otherDotArgs <- setdiff(names(allDotArgs), c(pdrDotArgs, pglrDotArgs))
    p1 <- do.call(
        what = plotDimRed,
        args = c(
            list(
                object = object,
                colorBy = useFactor,
                useDimRed = useDimRed,
                slot = "H.norm",
                zeroAsNA = TRUE,
                dotOrder = "asc",
                splitBy = NULL
            ),
            allDotArgs[pdrDotArgs],
            allDotArgs[otherDotArgs]
        )
    )
    bottom <- do.call(
        what = plotGeneLoadingRank,
        args = c(
            list(
                object = object,
                markerTable = markerTable,
                useFactor = useFactor,
                nLabel = nLabel,
                nPlot = nPlot
            ),
            allDotArgs[pglrDotArgs],
            allDotArgs[otherDotArgs]
        )
    )
    bottom <- bottom[c(1, 3, 2)]
    bottomComb <- cowplot::plot_grid(plotlist = bottom, nrow = 1)
    cowplot::plot_grid(p1, bottomComb, nrow = 2)
}

#' @rdname plotGeneLoadings
#' @export
plotGeneLoadingRank <- function(
        object,
        markerTable,
        useFactor,
        nLabel = 15,
        nPlot = 30,
        ...
) {
    # Table-object matching check
    dataset1 <- names(markerTable)[1]
    dataset1 <- .checkUseDatasets(object, useDatasets = dataset1)
    dataset2 <- names(markerTable)[3]
    dataset2 <- .checkUseDatasets(object, useDatasets = dataset2)

    geneList <- list()
    # subset to specific factor and sort by p-value
    geneList$V1 <- markerTable[[1]][markerTable[[1]]$factor_num == useFactor, ]
    geneList$V1 <- geneList$V1[order(geneList$V1$pval), ]$feature
    geneList$V2 <- markerTable[[3]][markerTable[[3]]$factor_num == useFactor, ]
    geneList$V2 <- geneList$V2[order(geneList$V2$pval), ]$feature
    # don't sort for W
    geneList$W <- markerTable[[2]][markerTable[[2]]$factor_num == useFactor,]$feature
    geneList <- lapply(geneList, function(g) g[seq(min(nLabel, length(g)))])

    loadingList <- getMatrix(object, "V", dataset = c(dataset1, dataset2))
    W <- getMatrix(object, "W")
    loadingList$W <- pmin(W + loadingList[[1]], W + loadingList[[2]])
    titles <- c(dataset1, dataset2, "Shared")
    plotList <- list()
    for (i in seq_along(loadingList)) {
        # loading: V or W matrix, g x k
        loading <- loadingList[[i]]
        topGenes <- geneList[[i]]
        sorted <- sort(loading[,useFactor])
        topLoaded <- rev(sorted)[seq(nPlot)]
        topLoaded <- names(topLoaded)
        topGenes <- topLoaded[topLoaded %in% topGenes]
        if (length(topGenes) == 0) topGenes <- "no genes"

        geneDF <- data.frame(loadings = sorted,
                             # Fix the xlimit within 0 - 1
                             # So together with ggplot2::coord_cartesian
                             # The top gene annotation can be shown at a fixed
                             # position
                             xpos = seq(0, 1, length.out = length(sorted)),
                             top_k = factor(names(sorted) %in% topGenes,
                                            levels = c(TRUE, FALSE)))

        ylimTxt <- max(geneDF$loadings)
        plotList[[i]] <- .ggScatter(geneDF, "xpos", "loadings", "top_k",
                                    title = titles[i], labelText = FALSE,
                                    colorValues = c("#8227A0", "black"), ...) +
            ggplot2::annotate(
                "text",
                x = 1.1,
                y = seq(ylimTxt, 0, length.out = nLabel)[seq_along(topGenes)],
                label = topGenes, hjust = 0, col = "#8227A0"
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                axis.ticks.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank(),
                axis.title.x = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                panel.grid.major.x = ggplot2::element_blank(),
                panel.grid.minor.x = ggplot2::element_blank(),
                legend.position = "none"
            ) +
            ggplot2::coord_cartesian(
                xlim = c(0, 1),
                clip = "off"
            ) +
            ggplot2::theme(plot.margin = grid::unit(c(1, 4, 1, 1), "lines"))
    }
    plotList
}


#' Make Riverplot/Sankey diagram that shows label mapping across datasets
#' @description
#' Creates a riverplot/Sankey diagram to show how independent cluster
#' assignments from two datasets map onto a joint clustering. Prior knowledge of
#' cell annotation for the given datasets is required to make sense from the
#' visualization. Dataset original annotation can be added with the syntax shown
#' in example code in this manual. The joint clustering could be generated with
#' \code{\link{runCluster}} or set by any other metadata annotation.
#'
#' Dataset original annotation can be inserted before running this function
#' using \code{cellMeta<-} method. Please see example below.
#'
#' This function depends on CRAN available package "sankey" and it has to be
#' installed in order to make this function work.
#' @export
#' @rdname plotSankey
#' @param object A \linkS4class{liger} object with all three clustering
#' variables available.
#' @param cluster1,cluster2 Name of the variables in \code{cellMeta(object)} for
#' the cluster assignments of dataset 1 and 2, respectively.
#' @param clusterConsensus Name of the joint cluster variable to use. Default
#' uses the default clustering of the object. Can select a variable name in
#' \code{cellMeta(object)}.
#' @param minFrac Numeric. Minimum fraction of cluster for an edge to be shown.
#' Default \code{0.05}.
#' @param minCell Numeric. Minimum number of cells for an edge to be shown.
#' Default \code{10}.
#' @param titles Character vector of three. Customizes the column title text
#' shown. Default uses the variable names \code{cluster1},
#' \code{clusterConsensus} and \code{cluster2}.
#' @param prefixes Character vector of three. Cluster names have to be unique
#' across all three variables, so this is provided to deduplicate the clusters
#' by adding \code{"prefixes[i]-"} before the actual label. This will not be
#' applied when no duplicate is found. Default \code{NULL} uses variable names.
#' An NA value or a string with no character (i.e. \code{""}) does not add the
#' prefix to the corresponding variable.
#' @param labelCex Numeric. Amount by which node label text should be magnified
#' relative to the default. Default \code{1}.
#' @param titleCex Numeric. Amount by which node label text should be magnified
#' relative to the default. Default \code{1.1}.
#' @param colorValues Character vector of color codes to set color for each
#' level in the consensus clustering. Default \code{scPalette}.
#' @param mar Numeric vector of the form \code{c(bottom, left, top, right)}
#' which gives the number of lines of margin to be specified on the four sides
#' of the plot. Increasing the 2nd and 4th values can be helpful when cluster
#' labels are long and extend out side of the plotting region. Default
#' \code{c(2, 2, 4, 2)}.
#' @note
#' This function works as a replacement of the function \code{makeRiverplot}
#' in rliger <1.99. We decide to make a new function because the dependency
#' adopted by the older version is archived on CRAN and will be no longer
#' available.
#' @returns No returned value. The sankey diagram will be displayed instead.
#' @examples
#' # Make fake dataset specific labels from joint clustering result
#' cellMeta(pbmcPlot, "ctrl_cluster", "ctrl") <-
#'     cellMeta(pbmcPlot, "leiden_cluster", "ctrl")
#' cellMeta(pbmcPlot, "stim_cluster", "stim") <-
#'     cellMeta(pbmcPlot, "leiden_cluster", "stim")
#' if (requireNamespace("sankey", quietly = TRUE)) {
#'     plotSankey(pbmcPlot, "ctrl_cluster", "stim_cluster",
#'                titles = c("control", "LIGER", "stim"),
#'                prefixes = c("c", NA, "s"))
#' }
plotSankey <- function(
        object,
        cluster1,
        cluster2,
        clusterConsensus = NULL,
        minFrac = 0.01,
        minCell = 10,
        titles = NULL,
        prefixes = NULL,
        labelCex = 1,
        titleCex = 1.1,
        colorValues = scPalette,
        mar = c(2, 2, 4, 2)
) {
    if (!requireNamespace("sankey", quietly = TRUE)) # nocov start
        cli::cli_abort(
            "Package {.pkg sankey} is needed for this function to work.
            Please install it by command:
            {.code install.packages('sankey')}") # nocov end

    clusterConsensus <- clusterConsensus %||% object@uns$defaultCluster
    clusterDF <- .fetchCellMetaVar(object,
                                   c(cluster1, clusterConsensus, cluster2),
                                   checkCategorical = TRUE, droplevels = TRUE)

    titles <- titles %||% c(cluster1, clusterConsensus, cluster2)
    titles <- .checkArgLen(titles, 3, repN = FALSE, class = "character")
    # Prepare for networkD3 input: Links, Nodes
    cluster1Fct <- droplevels(clusterDF[[1]])
    clusterCFct <- droplevels(clusterDF[[2]])
    cluster2Fct <- droplevels(clusterDF[[3]])
    # Have it separated so that the node matching don't have to be affected by
    # duplicated cluster names across variables
    nodes1 <- levels(cluster1Fct)
    nodesC <- levels(clusterCFct)
    nodes2 <- levels(cluster2Fct)
    .addPrefix <- function(p, n) {
        if (is.na(p) || nchar(p) == 0) n
        else paste0(p, '-', n)
    }

    if (any(duplicated(c(nodes1, nodesC, nodes2)))) {
        prefixes <- prefixes %||% c(cluster1, clusterConsensus, cluster2)
        prefixes <- .checkArgLen(prefixes, 3, repN = FALSE, class = "character")
        nodes1 <- .addPrefix(prefixes[1], nodes1)
        nodesC <- .addPrefix(prefixes[2], nodesC)
        nodes2 <- .addPrefix(prefixes[3], nodes2)
        levels(cluster1Fct) <- nodes1
        levels(clusterCFct) <- nodesC
        levels(cluster2Fct) <- nodes2
    }

    # Organize and filter the edges
    edges1 <- table(cluster1Fct, clusterCFct) %>%
        as.data.frame() %>%
        `colnames<-`(c("source", "target", "weight")) %>%
        # dplyr::filter(.data[["count"]] > 0)
        dplyr::group_by(.data[["source"]]) %>%
        dplyr::filter(.data[["weight"]] > minCell,
                      (.data[["weight"]] / sum(.data[["weight"]])) > minFrac) %>%
        dplyr::mutate(col = colorValues[as.integer(.data[["target"]])])
    edges2 <- table(clusterCFct, cluster2Fct) %>%
        as.data.frame() %>%
        `colnames<-`(c("source", "target", "weight")) %>%
        dplyr::group_by(.data[["target"]]) %>%
        dplyr::filter(.data[["weight"]] > minCell,
                      (.data[["weight"]] / sum(.data[["weight"]])) > minFrac) %>%
        dplyr::mutate(col = colorValues[as.integer(.data[["source"]])])

    # Remove unused nodes according to cleaned edges
    nodes1 <- nodes1[nodes1 %in% edges1$source]
    nodesC <- nodesC[nodesC %in% edges1$target | nodesC %in% edges2$source]
    nodes2 <- nodes2[nodes2 %in% edges2$target]

    edges <- rbind(edges1, edges2) %>%
        as.data.frame %>%
        dplyr::mutate(colorstyle = "col")
    edges[[1]] <- as.character(edges[[1]])
    edges[[2]] <- as.character(edges[[2]])
    nodes <- data.frame(id = c(nodes1, nodesC, nodes2),
                        x = c(rep(1, length(nodes1)),
                              rep(2, length(nodesC)),
                              rep(3, length(nodes2))),
                        col = "grey", cex = labelCex,
                        boxw = 0.05)

    pkgsnap_sankey <- sankey::make_sankey(nodes = nodes, edges = edges)
    sankey::sankey(pkgsnap_sankey, mar = mar) # mar order: bltr

    graphics::mtext(titles[1], side = 3, adj = 0.05, cex = titleCex, font = 2)
    graphics::mtext(titles[2], side = 3, adj = 0.5, cex = titleCex, font = 2)
    graphics::mtext(titles[3], side = 3, adj = 0.95, cex = titleCex, font = 2)
}

#' `r lifecycle::badge("deprecated")` Generate a river (Sankey) plot
#' @description
#' Creates a riverplot to show how separate cluster assignments from two
#' datasets map onto a joint clustering. The joint clustering is by default the
#' object clustering, but an external one can also be passed in. Uses the
#' riverplot package to construct riverplot object and then plot.
#' @param object \code{liger} object. Should run quantileAlignSNF before calling.
#' @param cluster1 Cluster assignments for dataset 1. Note that cluster names
#' should be distinct across datasets.
#' @param cluster2 Cluster assignments for dataset 2. Note that cluster names
#' should be distinct across datasets.
#' @param cluster_consensus Optional external consensus clustering (to use
#' instead of object clusters)
#' @param min.frac Minimum fraction of cluster for edge to be shown (default
#' 0.05).
#' @param min.cells Minumum number of cells for edge to be shown (default 10).
#' @param river.yscale y-scale to pass to riverplot -- scales the edge with
#' values by this factor, can be used to squeeze vertically (default 1).
#' @param river.lty Line style to pass to riverplot (default 0).
#' @param river.node_margin Node_margin to pass to riverplot -- how much
#' vertical space to keep between the nodes (default 0.1).
#' @param label.cex Size of text labels (default 1).
#' @param label.col Color of text labels (defualt "black").
#' @param lab.srt Angle of text labels (default 0).
#' @param river.usr Coordinates at which to draw the plot in form (x0, x1, y0,
#' y1).
#' @param node.order Order of clusters in each set (list with three vectors of
#' ordinal numbers). By default will try to automatically order them
#' appropriately.
#' @return \code{object} with refined cluster assignment updated in
#' \code{"louvain_cluster"} variable in \code{cellMeta} slot. Can be fetched
#' with \code{object$louvain_cluster}
#' @name makeRiverplot-deprecated
#' @seealso \code{\link{rliger-deprecated}}
NULL

#' @rdname rliger-deprecated
#' @section \code{makeRiverplot}:
#' For \code{makeRiverplot}, use \code{\link{plotSankey}} as the replacement.
#' @export
makeRiverplot <- function(object, cluster1, cluster2, cluster_consensus = NULL,
                          min.frac = 0.05, min.cells = 10, river.yscale = 1,
                          river.lty = 0, river.node_margin = 0.1, label.cex = 1,
                          label.col = "black", lab.srt = 0, river.usr = NULL,
                          node.order = "auto") {
    lifecycle::deprecate_stop("1.99.0", "makeRiverplot()",
                              "plotSankey()")
}


#' Visualize a spatial dataset
#' @description
#' Simple visualization of spatial coordinates. See example code for how to have
#' information preset in the object. Arguments to the liger object method are
#' passed down to ligerDataset method.
#' @export
#' @rdname plotSpatial
#' @param object Either a \linkS4class{liger} object containing a spatial
#' dataset or a \linkS4class{ligerSpatialDataset} object.
#' @inheritDotParams .ggScatter dotOrder dotSize dotAlpha raster labelTextSize seed
#' @inheritDotParams .ggplotLigerTheme title subtitle showLegend legendPosition baseSize titleSize xTitleSize yTitleSize legendTitleSize subtitleSize xTextSize yTextSize legendTextSize legendDotSize legendNRow legendNCol colorLabels colorValues naColor
#' @return A ggplot object
#' @examples
#' ctrl.fake.spatial <- as.ligerDataset(dataset(pbmc, "ctrl"), modal = "spatial")
#' fake.coords <- matrix(rnorm(2 * ncol(ctrl.fake.spatial)), ncol = 2)
#' coordinate(ctrl.fake.spatial) <- fake.coords
#' dataset(pbmc, "ctrl") <- ctrl.fake.spatial
#' defaultCluster(pbmc) <- pbmcPlot$leiden_cluster
#' plotSpatial2D(pbmc, dataset = "ctrl")
plotSpatial2D <- function(object, ...) {
    UseMethod("plotSpatial2D", object)
}

#' @export
#' @rdname plotSpatial
#' @method plotSpatial2D liger
#' @param dataset Name of one spatial dataset.
#' @param useCluster Either the name of one variable in \code{cellMeta(object)}
#' or a factor object with annotation that matches with all cells in the
#' specified dataset. Default \code{NULL} uses default clusters.
#' @param legendColorTitle Alternative title text in the legend. Default
#' \code{NULL} uses the variable name set by \code{useCluster}, or
#' \code{"Annotation"} is \code{useCluster} is a customized factor object.
plotSpatial2D.liger <- function(
        object,
        dataset,
        useCluster = NULL,
        legendColorTitle = NULL,
        ...) {
    dataset <- .checkUseDatasets(object, useDatasets = dataset,
                                 modal = "spatial")
    .checkArgLen(dataset, 1, class = "character")
    ld <- dataset(object, dataset)
    useCluster <- useCluster %||%
        defaultCluster(object)[object$dataset == dataset]
    if (length(useCluster) == 1) {
        legendColorTitle <- legendColorTitle %||% useCluster
        useCluster <- cellMeta(object, useCluster, useDatasets = dataset)
    } else {
        useCluster <- .checkArgLen(useCluster, ncol(ld), repN = FALSE, class = "factor")
        legendColorTitle <- legendColorTitle %||% "Annotation"
    }
    plotSpatial2D.ligerSpatialDataset(
        ld, useCluster = useCluster, legendColorTitle = legendColorTitle, ...
    )
}

#' @export
#' @rdname plotSpatial
#' @method plotSpatial2D ligerSpatialDataset
#' @param useDims Numeric vector of two, choosing the coordinates to be drawn
#' on 2D space. (STARmap data could have 3 dimensions.) Default \code{c(1, 2)}.
#' @param xlab,ylab Text label on x-/y-axis. Default \code{NULL} does not show
#' it.
#' @param labelText Logical, whether to label annotation onto the scatter plot.
#' Default \code{FALSE}.
#' @param panelBorder Whether to show rectangle border of the panel instead of
#' using ggplot classic bottom and left axis lines. Default \code{TRUE}.
plotSpatial2D.ligerSpatialDataset <- function(
        object,
        useCluster = NULL,
        legendColorTitle = NULL,
        useDims = c(1, 2),
        xlab = NULL,
        ylab = NULL,
        labelText = FALSE,
        panelBorder = TRUE,
        ...)
{
    .checkArgLen(useCluster, ncol(object), repN = FALSE, class = "factor")
    legendColorTitle <- legendColorTitle %||% "Annotation"

    coord <- coordinate(object)
    .checkArgLen(useDims, 2, repN = FALSE, class = "numeric")
    coord <- coord[, useDims]
    plotDF <- as.data.frame(coord)
    colnames(plotDF) <- c("x", "y")
    xRange <- c(min(plotDF$x), max(plotDF$x))
    yRange <- c(min(plotDF$y), max(plotDF$y))

    if (is.null(useCluster)) {
        .ggScatter(plotDF = plotDF, x = "x", y = "y",
                   xlab = xlab, ylab = ylab, panelBorder = panelBorder, ...) +
            ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank()) +
            ggplot2::coord_fixed(xlim = xRange, ylim = yRange)
    } else {
        plotDF[[legendColorTitle]] <- factor(useCluster)
        .ggScatter(plotDF = plotDF, x = "x", y = "y", colorBy = legendColorTitle,
                   xlab = xlab, ylab = ylab, labelText = labelText,
                   legendColorTitle = legendColorTitle,
                   panelBorder = panelBorder, ...) +
            ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank()) +
            ggplot2::coord_fixed(xlim = xRange, ylim = yRange)
    }
}
