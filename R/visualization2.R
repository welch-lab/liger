# Scatter Plots of DimRed ######################################################
#' Generate Dimensionality Reduction Plot with Coloring
#' @description some text
#' @param object A \linkS4class{liger} object.
#' @param useCluster Name of variable in \code{cell.meta} slot. Default
#' \code{"louvain_cluster"}.
#' @param useDimRed Name of type of the coordinate matrix in \code{cell.meta}
#' slot. Currently, choose from \code{"UMAP"} or \code{"TSNE"}. Default
#' \code{"UMAP"}.
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
#' \code{\link{.ggCellScatter}}, \code{\link{.ggplotLigerTheme}} for additional
#' graphic setting
#' @rdname plotClusterDimRed
#' @export
plotClusterDimRed <- function(
        object,
        useCluster = "louvain_cluster",
        useDimRed = c("UMAP", "TSNE"),
        ...) {
    useDimRed <- match.arg(useDimRed)
    xVar <- paste0(useDimRed, ".1")
    yVar <- paste0(useDimRed, ".2")
    plotCellScatter(object, x = xVar, y = yVar, colorBy = useCluster,
                    slot = "cell.meta", dotOrder = "shuffle", ...)
}

#' @rdname plotClusterDimRed
#' @export
plotDatasetDimRed <- function(
        object,
        useDimRed = c("UMAP", "TSNE"),
        ...) {
    useDimRed <- match.arg(useDimRed)
    xVar <- paste0(useDimRed, ".1")
    yVar <- paste0(useDimRed, ".2")
    plotCellScatter(object, x = xVar, y = yVar, colorBy = "dataset",
                    slot = "cell.meta", labelText = FALSE,
                    dotOrder = "shuffle", ...)
}

#' @rdname plotClusterDimRed
#' @export
plotByDatasetAndCluster <- function(
        object,
        useDimRed = c("UMAP", "TSNE"),
        useCluster = "louvain_cluster",
        combinePlots = TRUE,
        ...
) {
    useDimRed <- match.arg(useDimRed)
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

#' @rdname plotClusterDimRed
#' @export
plotGeneDimRed <- function(
        object,
        features,
        useDimRed = c("UMAP", "TSNE"),
        log = TRUE,
        scaleFactor = 1e4,
        zeroAsNA = TRUE,
        colorPalette = "C",
        ...
) {
    useDimRed <- match.arg(useDimRed)
    xVar <- paste0(useDimRed, ".1")
    yVar <- paste0(useDimRed, ".2")
    scaleFunc <- function(x) {
        if (!is.null(scaleFactor)) x <- scaleFactor*x
        if (isTRUE(log)) x <- log2(x + 1)
        x
    }
    plotCellScatter(object, x = xVar, y = yVar, colorBy = features,
                    slot = "norm.data", colorByFunc = scaleFunc,
                    dotOrder = "ascending", zeroAsNA = zeroAsNA,
                    colorPalette = colorPalette, ...)
}

#' @rdname plotClusterDimRed
#' @export
plotFactorDimRed <- function(
        object,
        factors,
        useDimRed = c("UMAP", "TSNE"),
        trimHigh = 0.03,
        zeroAsNA = TRUE,
        colorPalette = "D",
        ...
) {
    useDimRed <- match.arg(useDimRed)
    xVar <- paste0(useDimRed, ".1")
    yVar <- paste0(useDimRed, ".2")
    plotCellScatter(object, x = xVar, y = yVar, colorBy = factors,
                    slot = "H.norm", dotOrder = "ascending",
                    trimHigh = trimHigh, zeroAsNA = zeroAsNA,
                    colorPalette = colorPalette, ...)
}


# Violin Plots #################################################################
#' @export
plotGeneViolin <- function(
        object,
        gene,
        by.dataset = TRUE,
        groupBy = NULL,
        ...
) {
    splitBy <- NULL
    if (isTRUE(by.dataset)) splitBy <- "dataset"

    if (is.null(groupBy)) {
        if ("leiden_cluster" %in% names(cell.meta(object)))
            groupBy <- "leiden_cluster"
        else if ("louvain_cluster" %in% names(cell.meta(object)))
            groupBy <- "louvain_cluster"
        else if ("H.norm_cluster" %in% names(cell.meta(object)))
            groupBy <- "H.norm_cluster"
    } else if (isFALSE(groupBy)) {
        groupBy <- NULL
    }

    plotList <- plotCellViolin(
        object,
        y = gene,
        slot = "norm.data",
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

# Others #######################################################################

#' @export
plotProportionDot <- function(
        object,
        class1 = "louvain_cluster",
        class2 = "dataset",
        legend = FALSE,
        panelBorder = TRUE,
        ...
) {
    if (length(class1) != 1 ||
        length(class2) != 1)
        stop("`class1` and `class2` must be name of one categorical variable ",
             "in `cell.meta` slot.")
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
    .ggplotLigerTheme(p, legend = legend, panelBorder = panelBorder, ...)
}

#' @export
plotProportionBar <- function(
        object,
        class1 = "louvain_cluster",
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
             "in `cell.meta` slot.")
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

#' @export
plotClusterProportions <- function(
        object,
        useCluster = "louvain_cluster",
        return.plot = FALSE,
        ...
) {
    .deprecateArgs(defunct = "return.plot")
    plotProportionDot(object, class1 = useCluster, class2 = "dataset", ...)
}

#' @export
plotProportion <- function(
        object,
        class1 = "louvain_cluster",
        class2 = "dataset",
        method = c("stack", "group"),
        ...
) {
    method <- match.arg(method)
    p1 <- plotProportionDot(object, class1 = class1, class2 = class2, ...)
    p2 <- plotProportionBar(object, class1 = class1, class2 = class2,
                            inclRev = TRUE, combinePlot = TRUE,
                            method = method, ...)
    cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(1, 2))
}
