#' Generate scatter plot(s) using liger object
#' @description This function allows for using available cell metadata to build
#' the x-/y-axis. Available per-cell data can be used to form the color/shape
#' annotation, including cell metadata, raw or processed gene expression, and
#' unnormalized or aligned factor loading. Multiple coloring variable is allowed
#' from the same specification of \code{slot}, and this returns a list of plots
#' with different coloring values. Users can further split the plot(s) by
#' grouping on cells (e.g. datasets).
#' @details Available option for \code{slot} include: \code{"cell.meta"},
#' \code{"raw.data"}, \code{"norm.data"}, \code{"scale.data"}, \code{"H.norm"}
#' and \code{"H"}. When \code{"raw.data"}, \code{"norm.data"} or
#' \code{"scale.data"}, \code{colorBy} has to be a character vector of feature
#' names. When \code{"H.norm"} or \code{"H"}, \code{colorBy} can be any valid
#' index to select one factor of interests. Note that character index follows
#' \code{"factor_[k]"} format, with replacing \code{[k]} with an integer.
#'
#' When \code{"cell.meta"}, \code{colorBy} has to be an available column name in
#' the table. Note that, for \code{colorBy} as well as \code{x}, \code{y},
#' \code{shapeBy} and \code{splitBy}, since a matrix object is feasible in
#' \code{cell.meta} table, using a column (e.g. named as \code{"column1"} in a
#' certain matrix (e.g. named as \code{"matrixVar"}) should follow the syntax of
#' \code{"matrixVar.column1"}. When the matrix does not have a "colname"
#' attribute, the subscription goes with \code{"matrixVar.V1"},
#' \code{"matrixVar.V2"} and etc. Use \code{"UMAP.1"}, \code{"UMAP.2"},
#' \code{"TSNE.1"} or \code{"TSNE.2"} for the 2D embeddings generated with
#' rliger package. These are based on the nature of \code{as.data.frame} method
#' on a \code{\link[S4Vectors]{DataFrame}} object.
#' @param object \linkS4class{liger} object
#' @param x,y Available variable name in \code{cell.meta} slot to look for
#' the dot coordinates. See details.
#' @param colorBy Available variable name in specified \code{slot} to look for
#' color annotation information. See details. Default \code{NULL} generates
#' all-black dots.
#' @param slot Choose the slot to find the \code{colorBy} variable. See details.
#' Default \code{"cell.meta"}.
#' @param colorByFunc Default \code{NULL}. A function object that expects a
#' vector/factor/data.frame retrieved by \code{colorBy} as the only input, and
#' returns an object of the same size, so that the all color "aes" are replaced
#' by this output. Useful when, for example, users need to scale the gene
#' expression shown on plot.
#' @param splitBy Character vector of categorical variable names in
#' \code{cell.meta} slot. Split all cells by groupings on this/these variable(s)
#' to produce a scatter plot containing only the cells in each group. Default
#' \code{NULL}.
#' @param shapeBy Available variable name in \code{cell.meta} slot to look for
#' categorical annotation to be reflected by dot shapes. Default \code{NULL}.
#' @param ... More plot setting arguments. See \code{\link{.ggCellScatter}} and
#' \code{\link{.ggplotLigerTheme}}.
#' @return A ggplot object when a single plot is intended. A list of ggplot
#' objects, when multiple \code{colorBy} variables and/or \code{splitBy} are
#' set. When \code{plotly = TRUE}, all ggplot objects become plotly (htmlwidget)
#' objects.
#' @export
plotCellScatter <- function(
        object,
        x,
        y,
        colorBy = NULL,
        slot = c("cell.meta", "raw.data", "norm.data",
                 "scale.data", "H.norm", "H"),
        colorByFunc = NULL,
        splitBy = NULL,
        shapeBy = NULL,
        ...
) {
    slot <- match.arg(slot)
    plotDF <- as.data.frame(cell.meta(object))[, c(x, y, shapeBy, splitBy)]

    # Create copies of `plotDF` in `plotDFList`, where each `plotDF` has only
    # one `colorBy` variable
    plotDFList <- list()
    colorByParam <- list()
    if (!is.null(colorBy)) {
        colorDF <- retrieveCellFeature(object, colorBy, slot)
        # When retrieving H/H.norm, exact colname might not be what `colorBy` is
        colorBy <- colnames(colorDF)
        if (!is.null(colorByFunc))
            colorDF[, colorBy] <- colorByFunc(colorDF[, colorBy])
        for (i in seq_along(colorBy)) {
            plotDFList[[colorBy[i]]] <- cbind(plotDF,
                                              colorDF[, i, drop = FALSE])
            colorByParam[[colorBy[i]]] <- colorBy[i]
        }
    } else {
        plotDFList[[1]] <- plotDF
        colorByParam <- list(NULL)
    }

    # Split each `colorBy` specific `plotDF` in `plotDFList` by `splitBy`, so
    # `plotDFList` first becomes a nested list, where the first level is by
    # `colorBy`, and the second level is by `splitBy`
    if (!is.null(splitBy)) {
        for (i in seq_along(plotDFList)) {
            # Might be just one plotDF when 0/1 colorBy set, or more when >1
            # colorBy set
            plotDFList[[i]] <- split(plotDFList[[i]], plotDFList[[i]][,splitBy])
            names(plotDFList[[i]]) <- paste(names(plotDFList)[i],
                                            names(plotDFList[[i]]),
                                            sep = ".")
            colorByParam[[i]] <- rep(colorByParam[[i]],
                                     length(names(plotDFList[[i]])))
        }
        # Then flatten the nested list by concatenating each colorBy-list.
        plotDFList <- Reduce(c, plotDFList)
        colorByParam <- Reduce(c, colorByParam)
        names(colorByParam) <- names(plotDFList)
    }

    plotList <- list()
    for (i in seq_along(plotDFList)) {
        plotList[[i]] <- .ggCellScatter(plotDF = plotDFList[[i]], x = x, y = y,
                       colorBy = colorByParam[[i]], shapeBy = shapeBy, ...)
    }
    names(plotList) <- names(plotDFList)

    if (length(plotList) == 1) {
        return(plotList[[1]])
    } else {
        return(plotList)
    }
}

#' Produce single scatter plot with data frame passed from upstream
#' @param plotDF Data frame like object (fortifiable) that contains all
#' necessary information to make the plot.
#' @param x,y,colorBy,shapeBy See \code{\link{plotCellScatter}}.
#' @param dotOrder Controls the order that each dot is added to the plot. Choose
#' from \code{"shuffle"}, \code{"ascending"}, or \code{"descending"}. Default
#' \code{"shuffle"}, useful when coloring by categories that overlaps (e.g.
#' "dataset"), \code{"ascending"} can be useful when coloring by a continuous
#' variable (e.g. gene expression) where high values needs more
#' highlight. \code{NULL} use default order.
#' @param dotSize,dotAlpha Numeric, controls the size or transparency of all
#' dots. Default \code{0.6} and \code{0.9}.
#' @param raster Logical, whether to rasterize the plot. Default \code{NULL}
#' automatically rasterize the plot when number of total cells to be plotted
#' exceeds 100,000.
#' @param labelText Logical, whether to show text label at the median position
#' of each categorical group specified by \code{colorBy}. Default \code{TRUE}.
#' Does not work when continuous coloring is specified.
#' @param labelTextSize Numeric, controls the size of label size when
#' \code{labelText = TRUE}. Default \code{4}.
#' @param seed Random seed for reproducibility. Default \code{1}.
#' @param ... More theme setting arguments passed to
#' \code{\link{.ggplotLigerTheme}}.
#' @return ggplot object by default. When \code{plotly = TRUE}, returns
#' plotly (htmlwidget) object.
#' @importFrom rlang .data
.ggCellScatter <- function(
        plotDF,
        x,
        y,
        colorBy = NULL,
        shapeBy = NULL,
        dotOrder = c("shuffle", "ascending", "descending"),
        dotSize = 1,
        dotAlpha = 0.9,
        raster = NULL,
        labelText = TRUE,
        labelTextSize = 4,
        seed = 1,
        ...
) {
    dotOrder <- match.arg(dotOrder)
    set.seed(seed)
    if (is.null(raster)) {
        # Automatically decide whether to rasterize depending on number of cells
        if (nrow(plotDF) > 1e5) {
            raster <- TRUE
            .log("NOTE: Points are rasterized as number of cells/nuclei ",
                 "plotted exceeds 100,000.\n",
                 "Use `raster = FALSE` or `raster = TRUE` to force plot in ",
                 "vector form or not.")
        } else raster <- FALSE
    }

    if (!is.null(dotOrder)) {
        if (dotOrder == "shuffle") {
            idx <- sample(nrow(plotDF))
            plotDF <- plotDF[idx, ]
        } else if (dotOrder == "ascending") {
            plotDF <- plotDF[order(plotDF[[colorBy]], decreasing = FALSE),]
        } else {
            plotDF <- plotDF[order(plotDF[[colorBy]], decreasing = TRUE),]
        }
    }

    # TODO: Any way to avoid using these many conditions?
    if (!is.null(colorBy)) {
        if (!is.null(shapeBy))
            p <- ggplot2::ggplot(plotDF,
                                 ggplot2::aes(x = .data[[x]],
                                              y = .data[[y]],
                                              color = .data[[colorBy]],
                                              shape = .data[[shapeBy]]))
        else
            p <- ggplot2::ggplot(plotDF,
                                 ggplot2::aes(x = .data[[x]],
                                              y = .data[[y]],
                                              color = .data[[colorBy]]))
    } else {
        if (!is.null(shapeBy))
            p <- ggplot2::ggplot(plotDF,
                                 ggplot2::aes(x = .data[[x]],
                                              y = .data[[y]],
                                              shape = .data[[shapeBy]]))
        else
            p <- ggplot2::ggplot(plotDF,
                                 ggplot2::aes(x = .data[[x]],
                                              y = .data[[y]]))
    }
    if (isTRUE(raster)) {
        p <- p + scattermore::geom_scattermore(pointsize = dotSize,
                                               alpha = dotAlpha)
    } else {
        if (!is.logical(raster))
            warning("Please use TRUE/FALSE for `raster`. Not rasterizing.")
        p <- p + ggplot2::geom_point(size = dotSize, stroke = 0,
                                     alpha = dotAlpha)
    }

    # For categorical grouping
    if (!is.null(colorBy) && is.factor(plotDF[[colorBy]])) {
        if (isTRUE(labelText) & !is.null(colorBy)) {
            textData <- dplyr::group_by(plotDF, .data[[colorBy]])
            textData <- dplyr::summarise(textData,
                                         x = stats::median(.data[[x]]),
                                         y = stats::median(.data[[y]]))
            p <- p + ggrepel::geom_text_repel(
                data = textData,
                mapping = ggplot2::aes(x = .data[["x"]],
                                       y = .data[["y"]],
                                       label = .data[[colorBy]]),
                color = "black", size = labelTextSize, inherit.aes = FALSE
            )
            # Important to have `inherit.aes = F` above, otherwise
            # `geom_text_repel` looks for "shapeBy" setting which this newly
            # generated label coordinate table just doesn't have.
        }
    }

    p <- .ggplotLigerTheme(p, ...)

    return(p)
}

#' Generate violin/box plot(s) using liger object
#' @description This function allows for using available cell metadata, feature
#' expression or factor loading to generate violin plot, and grouping the data
#' with available categorical cell metadata. Available categorical cell metadata
#' can be used to form the color annotation. When it is different from the
#' grouping, it forms a nested grouping. Multiple y-axis variables are allowed
#' from the same specification of \code{slot}, and this returns a list of violin
#' plot for each. Users can further split the plot(s) by grouping on cells (e.g.
#' datasets).
#' @details Available option for \code{slot} include: \code{"cell.meta"},
#' \code{"raw.data"}, \code{"norm.data"}, \code{"scale.data"}, \code{"H.norm"}
#' and \code{"H"}. When \code{"raw.data"}, \code{"norm.data"} or
#' \code{"scale.data"}, \code{y} has to be a character vector of feature names.
#' When \code{"H.norm"} or \code{"H"}, \code{colorBy} can be any valid index to
#' select one factor of interests. Note that character index follows
#' \code{"factor_[k]"} format, with replacing \code{[k]} with an integer.
#'
#' When \code{"cell.meta"}, \code{y} has to be an available column name in
#' the table. Note that, for \code{y} as well as \code{groupBy}, \code{colorBy}
#' and \code{splitBy} since a matrix object is feasible in \code{cell.meta}
#' table, using a column (e.g. named as \code{"column1"} in a certain matrix
#' (e.g. named as \code{"matrixVar"}) should follow the syntax of
#' \code{"matrixVar.column1"}. When the matrix does not have a "colname"
#' attribute, the subscription goes with \code{"matrixVar.V1"},
#' \code{"matrixVar.V2"} and etc. These are based on the nature of
#' \code{as.data.frame} method on a \code{\link[S4Vectors]{DataFrame}} object.
#' @param object \linkS4class{liger} object
#' @param y Available variable name in \code{slot} to look for the value to
#' visualize.
#' @param groupBy,colorBy Available variable name in \code{cell.meta} slot to
#' look for categorical grouping. See details. Default \code{NULL} produces no
#' grouping and all-black graphic elements.
#' @param slot Choose the slot to find the \code{y} variable. See Details.
#' Default \code{"cell.meta"}.
#' @param yFunc A function object that expects a vector/factor/data.frame
#' retrieved by \code{y} as the only input, and returns an object of the same
#' size, so that the y-axis is replaced by this output. Useful when, for
#' example, users need to scale the gene expression shown on plot.
#' @param splitBy Character vector of categorical variable names in
#' \code{cell.meta} slot. Split all cells by groupings on this/these variable(s)
#' to produce a violin plot containing only the cells in each group. Default
#' \code{NULL}.
#' @param ... More plot setting arguments. See \code{\link{.ggCellViolin}} and
#' \code{\link{.ggplotLigerTheme}}.
#' @return A ggplot object when a single plot is intended. A list of ggplot
#' objects, when multiple \code{y} variables and/or \code{splitBy} are set. When
#' \code{plotly = TRUE}, all ggplot objects become plotly (htmlwidget) objects.
#' @export
plotCellViolin <- function(
        object,
        y,
        groupBy = NULL,
        slot = c("cell.meta", "raw.data", "norm.data",
                 "scale.data", "H.norm", "H"),
        yFunc = NULL,
        colorBy = NULL,
        splitBy = NULL,
        ...
) {
    slot <- match.arg(slot)

    # `groupBy` can be NULL, if so plotDF has ncell x 0 dimension
    if (is.null(groupBy)) {
        groupBy <- "All Cells"
        plotDF <- data.frame(factor(rep(NA, ncol(object))),
                             row.names = colnames(object))
        colnames(plotDF) <- groupBy
    } else {
        plotDF <- as.data.frame(cell.meta(object))[, groupBy, drop = FALSE]
    }
    plotDF[, splitBy] <- as.data.frame(cell.meta(object))[, splitBy,
                                                          drop = FALSE]
    if (!is.null(colorBy)) {
        colorDF <- as.data.frame(cell.meta(object))[, c(colorBy), drop = FALSE]
        plotDF[, colorBy] <- colorDF
    }

    plotDFList <- list()
    yParam <- list()

    # Create copies of `plotDF` in `plotDFList`, where each `plotDF` has only
    # one `y` variable
    yDF <- retrieveCellFeature(object, y, slot)
    # When retrieving H/H.norm, exact colname might not be what `colorBy` is
    y <- colnames(yDF)
    if (!is.null(yFunc))
        yDF[, y] <- yFunc(yDF[, y])
    for (i in seq_along(y)) {
        if (!y[i] %in% colnames(plotDF)) {
            plotDFList[[y[i]]] <- cbind(plotDF, yDF[, i, drop = FALSE])
        }
        yParam[[y[i]]] <- y[i]
    }

    # Split each `y` specific `plotDF` in `plotDFList` by `splitBy`, so
    # `plotDFList` first becomes a nested list, where the first level is by `y`,
    # and the secod level is by `splitBy`
    if (!is.null(splitBy)) {
        for (i in seq_along(plotDFList)) {
            plotDFList[[i]] <- split(plotDFList[[i]], plotDFList[[i]][,splitBy])
            names(plotDFList[[i]]) <- paste(names(plotDFList)[i],
                                            names(plotDFList[[i]]),
                                            sep = ".")
            yParam[[i]] <- rep(yParam[[i]], length(names(plotDFList[[i]])))
        }
        # Then flatten the nested list by concatenating each y-list.
        plotDFList <- Reduce(c, plotDFList)
        yParam <- Reduce(c, yParam)
        names(yParam) <- names(plotDFList)
    }

    plotList <- list()
    for (i in seq_along(plotDFList)) {
        plotList[[i]] <- .ggCellViolin(plotDF = plotDFList[[i]],
                                       y = yParam[[i]], groupBy = groupBy,
                                       colorBy = colorBy, ...)
    }
    names(plotList) <- names(plotDFList)

    if (length(plotList) == 1) {
        return(plotList[[1]])
    } else {
        return(plotList)
    }
}

#' Produce single violin plot with data frame passed from upstream
#' @param plotDF Data frame like object (fortifiable) that contains all
#' necessary information to make the plot.
#' @param y,groupBy,colorBy See \code{\link{plotCellViolin}}.
#' @param violin,box,dot Logical, whether to add violin plot, box plot or dot
#' (scatter) plot, respectively. Layers are added in the order of dot, violin,
#' and violin on the top surface. By default, only violin plot is generated.
#' @param violinAlpha,boxAlpha Numeric, controls the transparency of layers.
#' Default \code{0.8}, \code{0.6}, respectively.
#' @param violinWidth,boxWidth Numeric, controls the width of violin/box
#' bounding box. Default \code{0.9} and \code{0.4}.
#' @param dotColor,dotSize Numeric, globally controls the appearance of all
#' dots. Default \code{"black"} and \code{0.6}.
#' @param raster Logical, whether to rasterize the dot plot. Default \code{NULL}
#' automatically rasterizes the dot plot when number of total cells to be
#' plotted exceeds 100,000.
#' @param seed Random seed for reproducibility. Default \code{1}.
#' @param ... More theme setting arguments passed to
#' \code{\link{.ggplotLigerTheme}}.
#' @return ggplot object by default. When \code{plotly = TRUE}, returns
#' plotly (htmlwidget) object.
.ggCellViolin <- function(
        plotDF,
        y,
        groupBy = NULL,
        colorBy = NULL,
        violin = TRUE,
        violinAlpha = 0.8,
        violinWidth = 0.9,
        box = FALSE,
        boxAlpha = 0.6,
        boxWidth = 0.4,
        dot = FALSE,
        dotColor = "black",
        dotSize = 0.6,
        raster = NULL,
        seed = 1,
        ...
) {
    if (is.null(raster)) {
        # Automatically decide whether to rasterize depending on number of cells
        if (ncol(plotDF) > 1e5) {
            raster <- TRUE
            .log("NOTE: Points are rasterized as number of cells/nuclei ",
                 "plotted exceeds 100,000.\n",
                 "Use `raster = FALSE` or `raster = TRUE` to force plot in ",
                 "vector form or not.")
        } else raster <- FALSE
    }
    if (is.null(colorBy)) {
        plot <- ggplot2::ggplot(plotDF,
                                ggplot2::aes(x = .data[[groupBy]],
                                             y = .data[[y]]))
    } else {
        plot <- ggplot2::ggplot(plotDF,
                                ggplot2::aes(x = .data[[groupBy]],
                                             y = .data[[y]],
                                             colour = .data[[colorBy]],
                                             fill = .data[[colorBy]]))
    }

    if (isTRUE(dot)) {
        if (!is.null(dotColor)) {
            if (isTRUE(raster))
                plot <- plot + scattermore::geom_scattermore(
                    size = dotSize, color = dotColor, stroke = 0,
                    position = "jitter"
                )
            else
                plot <- plot + ggplot2::geom_jitter(size = dotSize,
                                                    color = dotColor,
                                                    stroke = 0, height = 0)
        } else {
            if (isTRUE(raster))
                plot <- plot + scattermore::geom_scattermore(
                    size = dotSize, stroke = 0, position = "jitter"
                )
            else
                plot <- plot + ggplot2::geom_jitter(size = dotSize,
                                                    stroke = 0, height = 0)
        }
    }
    if (isTRUE(violin))
        plot <- plot + ggplot2::geom_violin(alpha = violinAlpha,
                                            position = "dodge",
                                            width = violinWidth)
    if (isTRUE(box))
        plot <- plot + ggplot2::geom_boxplot(alpha = boxAlpha, fill = "white",
                                             position = "dodge",
                                             width = boxWidth)

    plot <- .ggplotLigerTheme(plot, ...)
    if (groupBy == "All Cells") {
        plot <- plot + ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                      axis.ticks.x = ggplot2::element_blank(),
                                      axis.title.x = ggplot2::element_blank())
    }
    plot
}

#' Generic ggplot theme setting for rliger package
#' @description Controls content and size of all peripheral texts.
#' @param plot ggplot object passed from wrapper plotting functions
#' @param title,subtitle,xlab,ylab Main title, subtitle or X/Y axis title text.
#' By default, no main title or subtitle will be set, and X/Y axis title will be
#' the names of variables used for plotting. Use \code{NULL} to hide elements.
#' \code{TRUE} for \code{xlab} or \code{ylab} shows default values.
#' @param legendColorTitle,legendShapeTitle Set alternative title text for
#' legend on color or shape. Default \code{NULL} shows the original variable
#' name.
#' @param legend Whether to show the legend. Default \code{TRUE}.
#' @param baseSize One-parameter control of all text sizes. Individual text
#' element sizes can be controlled by other size arguments. "Title" sizes are
#' 2 points larger than "text" sizes when being controlled by this.
#' @param titleSize,xTitleSize,yTitleSize,legendTitleSize Size of main title,
#' axis titles and legend title. Default \code{NULL} controls by
#' \code{baseSize + 2}.
#' @param subtitleSize,xTextSize,yTextSize,legendTextSize Size of subtitle text,
#' axis texts and legend text. Default \code{NULL} controls by \code{baseSize}.
#' @param legendDotSize Allow dots in legend region to be large enough to see
#' the colors/shapes clearly. Default \code{4}.
#' @param legendNRow,legendNCol Integer, when too many categories in one
#' variable, arranges number of rows or columns. Default \code{NULL},
#' automatically split to \code{ceiling(levels(variable)/10)} columns.
#' @param colorPalette For continuous coloring, a palette name or an index to
#' select from available options from ggplot
#' \code{\link[ggplot2]{scale_brewer}}. Additionally, \code{"viridis"} TODO are
#' also supported. Default \code{"viridis"}.
#' @param colorLow,colorMid,colorHigh,colorMidPoint All four of these must be
#' specified to customize palette with
#' \code{\link[ggplot2]{scale_colour_gradient2}}. Default \code{NULL}.
#' @param plotly Whether to use plotly to enable web based interactive browsing
#' for the plot. Requires installation of package "plotly". Default
#' \code{FALSE}.
#' @return Updated ggplot object by default. When \code{plotly = TRUE}, returns
#' plotly (htmlwidget) object.
.ggplotLigerTheme <- function(
        plot,
        # All text content
        title = NULL,
        subtitle = NULL,
        xlab = TRUE,
        ylab = TRUE,
        legendColorTitle = NULL,
        legendShapeTitle = NULL,
        legendSizeTitle = NULL,
        legend = TRUE,
        # All sizes
        baseSize = 10,
        titleSize = NULL,
        subtitleSize = NULL,
        xTextSize = NULL,
        xTitleSize = NULL,
        yTextSize = NULL,
        yTitleSize = NULL,
        legendTextSize = NULL,
        legendTitleSize = NULL,
        legendDotSize = 4,
        # Other
        legendNRow = NULL,
        legendNCol = NULL,
        colorPalette = "plasma",
        colorDirection = -1,
        colorLow = NULL,
        colorMid = NULL,
        colorHigh = NULL,
        colorMidPoint = NULL,
        plotly = FALSE
) {
    if (!is.null(title))
        plot <- plot + ggplot2::ggtitle(title, subtitle = subtitle)

    # Broadcast one-param setting to each
    Subtitle <- xText <- yText <- legendText <- baseSize
    Title <- xTitle <- yTitle <- legendTitle <- baseSize + 2
    # And set specific ones if specified
    if (!is.null(titleSize)) Title <- titleSize
    if (!is.null(subtitleSize)) Subtitle <- subtitleSize
    if (!is.null(xTextSize)) xText <- xTextSize
    if (!is.null(xTitleSize)) xTitle <- xTitleSize
    if (!is.null(yTextSize)) yText <- yTextSize
    if (!is.null(yTitleSize)) yTitle <- yTitleSize
    if (!is.null(legendTextSize)) legendText <- legendTextSize
    if (!is.null(legendTitleSize)) legendTitle <- legendTitleSize

    # Set x/y axis titles
    if (!isTRUE(xlab)) {
        if (is.null(xlab) || isFALSE(xlab)) {
            xlab <- NULL
        }
        plot <- plot + ggplot2::xlab(xlab)
    }
    if (!isTRUE(ylab)) {
        if (is.null(ylab) || isFALSE(ylab)) {
            ylab <- NULL
        }
        plot <- plot + ggplot2::ylab(ylab)
    }

    # Set sizes
    plot <- plot +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = Title),
            plot.subtitle = ggplot2::element_text(size = Subtitle),
            axis.text.x = ggplot2::element_text(size = xText),
            axis.title.x = ggplot2::element_text(size = xTitle),
            axis.text.y = ggplot2::element_text(size = yText),
            axis.title.y = ggplot2::element_text(size = yTitle),
            legend.text = ggplot2::element_text(size = legendText),
            legend.title = ggplot2::element_text(size = legendTitle)
        )

    # Whether to show legend
    if (!isTRUE(legend))
        plot <- plot + ggplot2::theme(legend.position = "none")
    else {
        # legend region settings. Need to prepare a list so we call
        # `guides()` once, otherwise any previous calls will be overwritten.
        guide <- list(colour = list(), shape = list(), size = list())
        guideFunc <- list(colour = NULL, shape = NULL, size = NULL)
        legendTitle <- list(colour = legendColorTitle,
                            shape = legendShapeTitle,
                            size = legendSizeTitle)
        for (a in names(guide)) {
            varName <- rlang::as_label(plot$mapping[[a]])
            if (varName != "NULL") {
                # If colour or shape is set
                if (is.factor(plot$data[[varName]])) {
                    # Categorical setting
                    guideFunc[[a]] <- ggplot2::guide_legend
                    # Set dot size in legend
                    guide[[a]]$override.aes <- list(size = legendDotSize)
                    # Set nrow/ncol to arrange the legend categories
                    if (is.null(legendNRow) & is.null(legendNCol)) {
                        # When nothing set, ggplot automatically makes it one
                        # column, which might be too long, so I add a auto-
                        # arrangement to limit max nrow to 10 but still evenly
                        # distribute the columns.
                        nCategory <- length(levels(plot$data[[varName]]))
                        if (nCategory > 10)
                            legendNCol <- ceiling(nCategory/10)
                    }
                    guide[[a]]$nrow <- legendNRow
                    guide[[a]]$ncol <- legendNCol
                } else {
                    # continuous setting
                    if (a == "colour") {
                        guideFunc[[a]] <- ggplot2::guide_colourbar
                        # Set continuous palette
                        plot <- plot +
                            .setColorBarPalette(palette = colorPalette,
                                                direction = colorDirection,
                                                low = colorLow, mid = colorMid,
                                                high = colorHigh,
                                                midPoint = colorMidPoint)
                    }
                    else
                        guideFunc[[a]] <- ggplot2::guide_legend

                }
                # Set title for the variable
                if (!is.null(legendTitle[[a]]))
                    guide[[a]]$title <- legendTitle[[a]]
                # Finalise the setting
                guide[[a]] <- do.call(guideFunc[[a]], guide[[a]])
            }
        }
        plot <- plot + ggplot2::guides(
            colour = guide$colour,
            shape = guide$shape
        )
    }

    if (isTRUE(plotly)) {
        if (requireNamespace("plotly", quietly = FALSE)) {
            plot <- plotly::ggplotly(plot)
        } else {
            warning('Run `install.packages("plotly")` to enable web based ',
                    "interactive browsing. Returning original ggplot.")
        }
    }
    return(plot)
}

.setColorBarPalette <- function(
        palette = "magma",
        direction = 1,
        low = NULL,
        mid = NULL,
        high = NULL,
        midPoint = NULL
) {
    viridisOptions <- c(
        "magma", "A", "inferno", "B", "plasma", "C", "viridis", "D",
        "cividis", "E", "rocket", "F", "mako", "G", "turbo", "H"
    )
    # TODO: discrete color palette
    layer <- NULL
    if (!is.null(low) & !is.null(mid) &
        !is.null(high) & !is.null(midPoint)) {
        # Only start to build customized color bar if all four arguments
        # are specified
        layer <- ggplot2::scale_colour_gradient2(low = low, mid = mid,
                                                 high = high,
                                                 midpoint = midPoint)
    } else if (!is.null(palette)) {
        # Otherwise, choose a palette based on non-NULL name
        if (palette %in% viridisOptions) {
            layer <- ggplot2::scale_colour_viridis_c(option = palette,
                                                     direction = direction)
        }
        else
            layer <- ggplot2::scale_colour_distiller(palette = palette)
    }
    # When nothing set, return NULL. "+ NULL" on a ggplot object doesn't
    # change anything
    return(layer)
}
