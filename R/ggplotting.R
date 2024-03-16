#' Generate scatter plot(s) using liger object
#' @description This function allows for using available cell metadata to build
#' the x-/y-axis. Available per-cell data can be used to form the color/shape
#' annotation, including cell metadata, raw or processed gene expression, and
#' unnormalized or aligned factor loading. Multiple coloring variable is allowed
#' from the same specification of \code{slot}, and this returns a list of plots
#' with different coloring values. Users can further split the plot(s) by
#' grouping on cells (e.g. datasets).
#' @details Available option for \code{slot} include: \code{"cellMeta"},
#' \code{"rawData"}, \code{"normData"}, \code{"scaleData"}, \code{"H.norm"}
#' and \code{"H"}. When \code{"rawData"}, \code{"normData"} or
#' \code{"scaleData"}, \code{colorBy} has to be a character vector of feature
#' names. When \code{"H.norm"} or \code{"H"}, \code{colorBy} can be any valid
#' index to select one factor of interests. Note that character index follows
#' \code{"Factor_[k]"} format, with replacing \code{[k]} with an integer.
#'
#' When \code{"cellMeta"}, \code{colorBy} has to be an available column name in
#' the table. Note that, for \code{colorBy} as well as \code{x}, \code{y},
#' \code{shapeBy} and \code{splitBy}, since a matrix object is feasible in
#' \code{cellMeta} table, using a column (e.g. named as \code{"column1"} in a
#' certain matrix (e.g. named as \code{"matrixVar"}) should follow the syntax of
#' \code{"matrixVar.column1"}. When the matrix does not have a "colname"
#' attribute, the subscription goes with \code{"matrixVar.V1"},
#' \code{"matrixVar.V2"} and etc. Use \code{"UMAP.1"}, \code{"UMAP.2"},
#' \code{"TSNE.1"} or \code{"TSNE.2"} for the 2D embeddings generated with
#' rliger package. These are based on the nature of \code{as.data.frame} method
#' on a \code{\link[S4Vectors]{DataFrame}} object.
#' @param object A \linkS4class{liger} object
#' @param colorBy Available variable name in specified \code{slot} to look for
#' color annotation information. See details. Default \code{NULL} generates
#' all-black dots.
#' @param slot Choose the slot to find the \code{colorBy} variable. See details.
#' Default \code{"cellMeta"}.
#' @param colorByFunc Default \code{NULL}. A function object that expects a
#' vector/factor/data.frame retrieved by \code{colorBy} as the only input, and
#' returns an object of the same size, so that the all color "aes" are replaced
#' by this output. Useful when, for example, users need to scale the gene
#' expression shown on plot.
#' @param cellIdx Character, logical or numeric index that can subscribe cells.
#' Missing or \code{NULL} for all cells.
#' @param splitBy Character vector of categorical variable names in
#' \code{cellMeta} slot. Split all cells by groupings on this/these variable(s)
#' to produce a scatter plot containing only the cells in each group. Default
#' \code{NULL}.
#' @param shapeBy Available variable name in \code{cellMeta} slot to look for
#' categorical annotation to be reflected by dot shapes. Default \code{NULL}.
#' @param titles Title text. A character scalar or a character vector with as
#' many elements as multiple plots are supposed to be generated. Default
#' \code{NULL}.
#' @param ... More plot setting arguments. See \code{\link{.ggScatter}} and
#' \code{\link{.ggplotLigerTheme}}.
#' @return A ggplot object when a single plot is intended. A list of ggplot
#' objects, when multiple \code{colorBy} variables and/or \code{splitBy} are
#' set. When \code{plotly = TRUE}, all ggplot objects become plotly (htmlwidget)
#' objects.
#' @export
#' @examples
#' plotDimRed(pbmcPlot, colorBy = "dataset", slot = "cellMeta",
#'            labelText = FALSE)
#' plotDimRed(pbmcPlot, colorBy = "S100A8", slot = "normData",
#'            dotOrder = "ascending", dotSize = 2)
#' plotDimRed(pbmcPlot, colorBy = 2, slot = "H.norm",
#'            dotOrder = "ascending", dotSize = 2, colorPalette = "viridis")
plotDimRed <- function(
        object,
        colorBy = NULL,
        useDimRed = NULL,
        slot = c("cellMeta", "rawData", "normData",
                 "scaleData", "H.norm", "H",
                 "normPeak", "rawPeak"),
        colorByFunc = NULL,
        cellIdx = NULL,
        splitBy = NULL,
        shapeBy = NULL,
        titles = NULL,
        ...
) {
    slot <- match.arg(slot)
    # useDimRed <- useDimRed %||% object@uns$defaultDimRed
    # useDimRed <- .findDimRedName(object, useDimRed, stopOnNull = TRUE, returnFirst = TRUE)
    plotDF <- as.data.frame(dimRed(object, useDimRed, cellIdx = cellIdx))
    x <- colnames(plotDF)[1]
    y <- colnames(plotDF)[2]

    ann <- .fetchCellMetaVar(object, variables = c(shapeBy, splitBy),
                             checkCategorical = TRUE, cellIdx = cellIdx,
                             drop = FALSE, droplevels = TRUE)
    if (!is.null(ann)) plotDF <- cbind(plotDF, ann)
    cellIdx <- .idxCheck(object, cellIdx, orient = "cell")
    # Create copies of `plotDF` in `plotDFList`, where each `plotDF` has only
    # one `colorBy` variable
    plotDFList <- list()
    colorByParam <- list()
    if (!is.null(colorBy)) {
        colorDF <- retrieveCellFeature(object, feature = colorBy,
                                       slot = slot, cellIdx = cellIdx,
                                       verbose = FALSE)
        # When retrieving H/H.norm, exact colname might not be what `colorBy` is
        colorBy <- colnames(colorDF)
        if (!is.null(colorByFunc))
            colorDF[, colorBy] <- colorByFunc(colorDF[, colorBy])
        for (i in seq_along(colorBy)) {
            if (!colorBy[i] %in% colnames(plotDF)) {
                #subDF <- cbind(plotDF, colorDF[, i, drop = FALSE])
                plotDFList[[colorBy[i]]] <- cbind(plotDF,
                                                  colorDF[, i, drop = FALSE])
            } else {
                plotDFList[[colorBy[i]]] <- plotDF
                plotDFList[[colorBy[i]]][[colorBy[i]]] <- colorDF[, i]
            }
            # plotDFList[[colorBy[i]]] <- subDF[cellIdx, , drop = FALSE]
            colorByParam[[colorBy[i]]] <- colorBy[i]
        }

    } else {
        plotDFList[[1]] <- plotDF#[cellIdx, , drop = FALSE]
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
                                     length(plotDFList[[i]]))
        }
        # Then flatten the nested list by concatenating each colorBy-list.
        plotDFList <- Reduce(c, plotDFList)
        colorByParam <- Reduce(c, colorByParam)
        if (!is.null(colorByParam)) {
            names(colorByParam) <- names(plotDFList)
        }
    }

    plotList <- list()
    titles <- .checkArgLen(titles, n = length(plotDFList), class = "ANY", .stop = FALSE)
    for (i in seq_along(plotDFList)) {
        cliID <- cli::cli_process_start("Plotting feature {.val {names(plotDFList)[i]}} on {.val {nrow(plotDFList[[i]])}} cells")
        plotList[[i]] <- .ggScatter(plotDF = plotDFList[[i]], x = x, y = y,
                                    colorBy = colorByParam[[i]],
                                    shapeBy = shapeBy, title = titles[i], ...)
        cli::cli_process_done(cliID)
    }
    names(plotList) <- names(plotDFList)

    if (length(plotList) == 1) {
        return(plotList[[1]])
    } else {
        return(plotList)
    }
}

#' Produce single scatter plot with data frame passed from upstream
#' @details Having package "ggrepel" installed can help adding tidier text
#' labels on the scatter plot.
#' @param plotDF Data frame like object (fortifiable) that contains all
#' necessary information to make the plot.
#' @param x,y Available variable name in \code{cellMeta} slot to look for
#' the dot coordinates. See details.
#' @param colorBy,shapeBy See \code{\link{plotDimRed}}.
#' @param dotOrder Controls the order that each dot is added to the plot. Choose
#' from \code{"shuffle"}, \code{"ascending"}, or \code{"descending"}. Default
#' \code{"shuffle"}, useful when coloring by categories that overlaps (e.g.
#' "dataset"), \code{"ascending"} can be useful when coloring by a continuous
#' variable (e.g. gene expression) where high values needs more
#' highlight. \code{NULL} use default order.
#' @param dotSize,dotAlpha Numeric, controls the size or transparency of all
#' dots. Default \code{getOption("ligerDotSize")} (1) and \code{0.9}.
#' @param trimHigh,trimLow Numeric, limit the largest or smallest value of
#' continuous \code{colorBy} variable. Default \code{NULL}.
#' @param zeroAsNA Logical, whether to set zero values in continuous
#' \code{colorBy} variable to \code{NA} so the color of these value.
#' @param raster Logical, whether to rasterize the plot. Default \code{NULL}
#' automatically rasterize the plot when number of total cells to be plotted
#' exceeds 100,000.
#' @param labelBy A variable name available in \code{plotDF}. If the variable is
#' categorical (a factor), the label position will be the median coordinates of
#' all dots within the same group. Unique labeling in character vector for each
#' dot is also acceptable. Default \code{colorBy}.
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
.ggScatter <- function(
        plotDF,
        x,
        y,
        colorBy = NULL,
        shapeBy = NULL,
        dotOrder = c("shuffle", "ascending", "descending"),
        dotSize = getOption("ligerDotSize"),
        dotAlpha = 0.9,
        trimHigh = NULL,
        trimLow = NULL,
        zeroAsNA = TRUE,
        raster = NULL,
        labelBy = colorBy,
        labelText = TRUE,
        labelTextSize = 4,
        seed = 1,
        ...
) {
    dotOrder <- match.arg(dotOrder)
    set.seed(seed)
    raster <- .checkRaster(nrow(plotDF), raster)
    if (!is.null(colorBy)) {
        if (dotOrder == "shuffle") {
            # Always put NA at bottom layer, i.e. plot them first
            isNA <- which(is.na(plotDF[[colorBy]]))
            nonNA <- which(!is.na(plotDF[[colorBy]]))
            idx <- sample(nonNA)
            plotDF <- plotDF[c(isNA, idx), ]
        } else if (dotOrder == "ascending") {
            plotDF <- plotDF[order(plotDF[[colorBy]], decreasing = FALSE,
                                   na.last = FALSE),]
        } else {
            plotDF <- plotDF[order(plotDF[[colorBy]], decreasing = TRUE,
                                   na.last = FALSE),]
        }
        if (!is.factor(plotDF[[colorBy]])) {
            if (!is.null(trimHigh))
                plotDF[[colorBy]][plotDF[[colorBy]] > trimHigh] <- trimHigh
            if (!is.null(trimLow))
                plotDF[[colorBy]][plotDF[[colorBy]] < trimLow] <- trimLow
            if (isTRUE(zeroAsNA))
                plotDF[[colorBy]][plotDF[[colorBy]] == 0] <- NA
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
        p <- p + ggplot2::geom_point(size = dotSize, stroke = 0,
                                     alpha = dotAlpha)
    }

    # For categorical grouping
    if (!is.null(labelBy) &&
        (is.factor(plotDF[[labelBy]]) || is.character(plotDF[[labelBy]]))) {
        if (isTRUE(labelText)) {
            textData <- dplyr::group_by(plotDF, .data[[labelBy]])
            textData <- dplyr::summarise(textData,
                                         x = stats::median(.data[[x]]),
                                         y = stats::median(.data[[y]]))
            if (!requireNamespace("ggrepel", quietly = TRUE)) {
                p <- p + ggplot2::annotate(
                    "text", x = textData$x, y = textData$y,
                    label = textData[[labelBy]], color = "black",
                    size = labelTextSize
                )
            } else {
                p <- p + ggrepel::geom_text_repel(
                    data = textData,
                    mapping = ggplot2::aes(x = .data[["x"]],
                                           y = .data[["y"]],
                                           label = .data[[labelBy]]),
                    color = "black", size = labelTextSize, inherit.aes = FALSE,
                    bg.colour = "white", bg.r = .2
                )
            }

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
#' @details Available option for \code{slot} include: \code{"cellMeta"},
#' \code{"rawData"}, \code{"normData"}, \code{"scaleData"}, \code{"H.norm"}
#' and \code{"H"}. When \code{"rawData"}, \code{"normData"} or
#' \code{"scaleData"}, \code{y} has to be a character vector of feature names.
#' When \code{"H.norm"} or \code{"H"}, \code{colorBy} can be any valid index to
#' select one factor of interests. Note that character index follows
#' \code{"Factor_[k]"} format, with replacing \code{[k]} with an integer.
#'
#' When \code{"cellMeta"}, \code{y} has to be an available column name in
#' the table. Note that, for \code{y} as well as \code{groupBy}, \code{colorBy}
#' and \code{splitBy} since a matrix object is feasible in \code{cellMeta}
#' table, using a column (e.g. named as \code{"column1"} in a certain matrix
#' (e.g. named as \code{"matrixVar"}) should follow the syntax of
#' \code{"matrixVar.column1"}. When the matrix does not have a "colname"
#' attribute, the subscription goes with \code{"matrixVar.V1"},
#' \code{"matrixVar.V2"} and etc. These are based on the nature of
#' \code{as.data.frame} method on a \code{\link[S4Vectors]{DataFrame}} object.
#'
#' \code{groupBy} is basically send to \code{ggplot2::aes(x)}, while
#' \code{colorBy} is for the "colour" aesthetics. Specifying \code{colorBy}
#' without \code{groupBy} visually creates grouping but there will not be
#' varying values on the x-axis, so \code{boxWidth} will be forced to the same
#' value as \code{violinWidth} under this situation.
#' @param object \linkS4class{liger} object
#' @param y Available variable name in \code{slot} to look for the value to
#' visualize.
#' @param groupBy,colorBy Available variable name in \code{cellMeta} slot to
#' look for categorical grouping. See details. Default \code{NULL} produces no
#' grouping and all-black graphic elements.
#' @param slot Choose the slot to find the \code{y} variable. See Details.
#' Default \code{"cellMeta"}.
#' @param yFunc A function object that expects a vector/factor/data.frame
#' retrieved by \code{y} as the only input, and returns an object of the same
#' size, so that the y-axis is replaced by this output. Useful when, for
#' example, users need to scale the gene expression shown on plot.
#' @param cellIdx Character, logical or numeric index that can subscribe cells.
#' Missing or \code{NULL} for all cells.
#' @param splitBy Character vector of categorical variable names in
#' \code{cellMeta} slot. Split all cells by groupings on this/these variable(s)
#' to produce a violin plot containing only the cells in each group. Default
#' \code{NULL}.
#' @param titles Title text. A character scalar or a character vector with as
#' many elements as multiple plots are supposed to be generated. Default
#' \code{NULL}.
#' @param ... More plot setting arguments. See \code{\link{.ggCellViolin}} and
#' \code{\link{.ggplotLigerTheme}}.
#' @return A ggplot object when a single plot is intended. A list of ggplot
#' objects, when multiple \code{y} variables and/or \code{splitBy} are set. When
#' \code{plotly = TRUE}, all ggplot objects become plotly (htmlwidget) objects.
#' @export
#' @examples
#' plotCellViolin(pbmcPlot, y = "nUMI", groupBy = "dataset", slot = "cellMeta")
#' plotCellViolin(pbmcPlot, y = "nUMI", groupBy = "leiden_cluster",
#'                slot = "cellMeta", splitBy = "dataset",
#'                colorBy = "leiden_cluster",
#'                box = TRUE, dot = TRUE,
#'                ylab = "Total counts per cell",
#'                colorValues = RColorBrewer::brewer.pal(8, "Set1"))
#' plotCellViolin(pbmcPlot, y = "S100A8", slot = "normData",
#'                yFunc = function(x) log2(10000*x + 1),
#'                groupBy = "dataset", colorBy = "leiden_cluster",
#'                box = TRUE, ylab = "S100A8 Expression")
plotCellViolin <- function(
        object,
        y,
        groupBy = NULL,
        slot = c("cellMeta", "rawData", "normData",
                 "scaleData", "H.norm", "H"),
        yFunc = NULL,
        cellIdx = NULL,
        colorBy = NULL,
        splitBy = NULL,
        titles = NULL,
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
        plotDF <- .fetchCellMetaVar(object, groupBy, checkCategorical = TRUE,
                                    drop = FALSE, cellIdx = cellIdx,
                                    droplevels = TRUE)
    }
    plotDF[,splitBy] <- .fetchCellMetaVar(object, splitBy, cellIdx = cellIdx,
                                          checkCategorical = TRUE, drop = FALSE,
                                          droplevels = TRUE)
    if (!is.null(colorBy))
        plotDF[,colorBy] <- .fetchCellMetaVar(object, colorBy,
                                              cellIdx = cellIdx,
                                              checkCategorical = TRUE)

    plotDFList <- list()
    yParam <- list()

    # Create copies of `plotDF` in `plotDFList`, where each `plotDF` has only
    # one `y` variable
    yDF <- retrieveCellFeature(object, y, slot, cellIdx = cellIdx,
                               verbose = FALSE)

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
    titles <- .checkArgLen(titles, n = length(plotDFList), class = "ANY", .stop = FALSE)
    for (i in seq_along(plotDFList)) {
        plotList[[i]] <- .ggCellViolin(plotDF = plotDFList[[i]],
                                       y = yParam[[i]], groupBy = groupBy,
                                       colorBy = colorBy, title = titles[i],
                                       ...)
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
#' dots. Default \code{"black"} and \code{getOption("ligerDotSize")} (1).
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
        dotSize = getOption("ligerDotSize"),
        raster = NULL,
        seed = 1,
        ...
) {
    raster <- .checkRaster(nrow(plotDF), raster)
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
        if (!identical(colorBy, groupBy)) boxWidth <- violinWidth
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
#' @param legendColorTitle,legendFillTitle,legendShapeTitle,legendSizeTitle Set
#' alternative title text for legend on aes of color, fill, shape and size,
#' respectively. Default \code{NULL} shows the original variable name.
#' @param showLegend Whether to show the legend. Default \code{TRUE}.
#' @param legendPosition Text indicating where to place the legend. Choose from
#' \code{"top"}, \code{"bottom"}, \code{"left"} or \code{"right"}. Default
#' \code{"right"}.
#' @param baseSize One-parameter control of all text sizes. Individual text
#' element sizes can be controlled by other size arguments. "Title" sizes are
#' 2 points larger than "text" sizes when being controlled by this.
#' @param titleSize,xTitleSize,yTitleSize,legendTitleSize Size of main title,
#' axis titles and legend title. Default \code{NULL} controls by
#' \code{baseSize + 2}.
#' @param subtitleSize,xTextSize,yTextSize,legendTextSize Size of subtitle text,
#' axis texts and legend text. Default \code{NULL} controls by \code{baseSize}.
#' @param xFacetSize,yFacetSize Size of facet label text. Default \code{NULL}
#' controls by \code{baseSize - 2}.
#' @param legendDotSize Allow dots in legend region to be large enough to see
#' the colors/shapes clearly. Default \code{4}.
#' @param panelBorder Whether to show rectangle border of the panel instead of
#' using ggplot classic bottom and left axis lines. Default \code{FALSE}.
#' @param colorLabels,colorValues Each a vector with as many values as the
#' number of categories for the categorical coloring aesthetics. Labels will be
#' the shown text and values will be the color code. These are passed to
#' \code{\link[ggplot2]{scale_color_manual}}. Default uses an internal selected
#' palette if there are <= 26 colors needed, or ggplot hues otherwise, and plot
#' original labels (levels of the factor).
#' @param legendNRow,legendNCol Integer, when too many categories in one
#' variable, arranges number of rows or columns. Default \code{NULL},
#' automatically split to \code{ceiling(levels(variable)/10)} columns.
#' @param colorPalette For continuous coloring, an index or a palette name to
#' select from available options from ggplot
#' \code{\link[ggplot2]{scale_brewer}} or \code{\link[viridisLite]{viridis}}.
#' Default \code{"magma"}.
#' @param colorDirection Choose \code{1} or \code{-1}. Applied when
#' \code{colorPalette} is from Viridis options. Default \code{-1} use darker
#' color for higher value, while \code{1} reverses this direction.
#' @param colorLow,colorMid,colorHigh,colorMidPoint All four of these must be
#' specified to customize palette with
#' @param naColor The color code for \code{NA} values. Default \code{"#DEDEDE"}.
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
        legendFillTitle = NULL,
        legendShapeTitle = NULL,
        legendSizeTitle = NULL,
        showLegend = TRUE,
        legendPosition = "right",
        # All sizes
        baseSize = getOption("ligerBaseSize"),
        titleSize = NULL,
        subtitleSize = NULL,
        xTextSize = NULL,
        xFacetSize = NULL,
        xTitleSize = NULL,
        yTextSize = NULL,
        yFacetSize = NULL,
        yTitleSize = NULL,
        legendTextSize = NULL,
        legendTitleSize = NULL,
        legendDotSize = 4,
        # Other
        panelBorder = FALSE,
        legendNRow = NULL,
        legendNCol = NULL,
        # Coloring
        colorLabels = NULL,
        colorValues = NULL,
        colorPalette = "magma",
        colorDirection = -1,
        naColor = "#DEDEDE",
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
    xFacet <- yFacet <- baseSize - 2
    # And set specific ones if specified
    if (!is.null(titleSize)) Title <- titleSize
    if (!is.null(subtitleSize)) Subtitle <- subtitleSize
    if (!is.null(xTextSize)) xText <- xTextSize
    if (!is.null(xFacetSize)) xFacet <- xFacetSize
    if (!is.null(xTitleSize)) xTitle <- xTitleSize
    if (!is.null(yTextSize)) yText <- yTextSize
    if (!is.null(yFacetSize)) yFacet <- yFacetSize
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
            strip.text.x = ggplot2::element_text(size = xFacet),
            strip.text.y = ggplot2::element_text(size = yFacet),
            legend.text = ggplot2::element_text(size = legendText),
            legend.title = ggplot2::element_text(size = legendTitle)
        )

    if (isTRUE(panelBorder)) {
        plot <- plot + ggplot2::theme(
            axis.line = ggplot2::element_line(linewidth = 0),
            panel.border = ggplot2::element_rect(fill = NA, colour = "black",
                                                 linewidth = 0.7)
        )
    }

    # legend region settings. Need to prepare a list so we call
    # `guides()` once, otherwise any previous calls will be overwritten.
    guide <- list(colour = list(), fill = list(),
                  shape = list(), size = list())
    guideFunc <- list(colour = NULL, shape = NULL, size = NULL)
    legendTitle <- list(colour = legendColorTitle,
                        fill = legendFillTitle,
                        shape = legendShapeTitle,
                        size = legendSizeTitle)
    for (a in names(guide)) {
        varName <- rlang::as_label(plot$mapping[[a]])
        if (varName == "NULL") next
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
            if (is.null(colorLabels)) {
                colorLabels <- levels(plot$data[[varName]])
            }
            if (is.null(colorValues)) {
                if (nlevels(plot$data[[varName]]) <= length(scPalette))
                    colorValues <- scPalette[seq_len(nlevels(plot$data[[varName]]))]
                else {
                    colorValues <- scales::hue_pal()(
                        length(levels(plot$data[[varName]]))
                    )
                }
            }
            if (a %in% c("colour", "fill")) {
                plot <- plot +
                    .setColorLegendPalette(plot$data[[varName]],
                                           aesType = a,
                                           labels = colorLabels,
                                           values = colorValues,
                                           naColor = naColor)
            }
        } else {
            # continuous setting
            if (a %in% c("colour", "fill")) {
                guideFunc[[a]] <- ggplot2::guide_colourbar
                # Set continuous palette
                plot <- plot +
                    .setColorBarPalette(aesType = a,
                                        palette = colorPalette,
                                        direction = colorDirection,
                                        naColor = naColor,
                                        low = colorLow, mid = colorMid,
                                        high = colorHigh,
                                        midPoint = colorMidPoint)
            } else {
                guideFunc[[a]] <- ggplot2::guide_legend
            }

            if (a == "size") {
                guide[[a]]$fill <- "black"
            }

        }
        # Set title for the variable
        if (!is.null(legendTitle[[a]]))
            guide[[a]]$title <- legendTitle[[a]]
        # Finalise the setting
        guide[[a]] <- do.call(guideFunc[[a]], guide[[a]])
    }
    guide <- lapply(guide, function(x) if (!identical(x, list())) x else NULL)
    plot <- plot +
        ggplot2::guides(
            colour = guide$colour,
            shape = guide$shape,
            size = guide$size,
            fill = guide$fill
        ) +
        ggplot2::theme(
            legend.position = legendPosition
        )

    if (isFALSE(showLegend)) {
        plot <- plot + ggplot2::theme(legend.position = "none")
    }

    if (isTRUE(plotly)) {
        if (requireNamespace("plotly", quietly = TRUE)) {
            plot <- plotly::ggplotly(plot)
        } else {
            cli::cli_alert_danger(
                "Package {.pkg plotly} is needed for interactive browsing."
            )
            cli::cli_alert_info("Please run {.code install.packages('plotly')} to enable it.")
            cli::cli_alert_info("Returning the original {.cls ggplot}.")
        }
    }
    return(plot)
}

.setColorLegendPalette <- function(
        fct,
        aesType = c("colour", "fill"),
        labels = NULL,
        values = NULL,
        naColor = "#DEDEDE"
) {
    aesType <- match.arg(aesType)
    layer <- NULL
    lvlUse <- table(fct) > 0
    labels <- labels[lvlUse]
    values <- values[lvlUse]
    if (!is.null(labels) && !is.null(values)) {
        if (aesType == "colour") {
            layer <- ggplot2::scale_color_manual(values = values,
                                                 labels = labels,
                                                 na.value = naColor)
        } else {
            layer <- ggplot2::scale_fill_manual(values = values,
                                                labels = labels,
                                                na.value = naColor)
        }
    }
    return(layer)
}

.setColorBarPalette <- function(
        aesType = c("colour", "fill"),
        palette = "magma",
        direction = 1,
        naColor = "#DEDEDE",
        low = NULL,
        mid = NULL,
        high = NULL,
        midPoint = NULL
) {
    aesType <- match.arg(aesType)
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
                                                 midpoint = midPoint,
                                                 na.value = naColor)
    } else if (!is.null(palette)) {
        # Otherwise, choose a palette based on non-NULL name
        if (palette %in% viridisOptions) {
            if (aesType == "colour")
                layer <- ggplot2::scale_colour_viridis_c(option = palette,
                                                         direction = direction,
                                                         na.value = naColor)
            else
                layer <- ggplot2::scale_fill_viridis_c(option = palette,
                                                       direction = direction,
                                                       na.value = naColor)
        }
        else
            if (aesType == "colour")
                layer <- ggplot2::scale_colour_distiller(palette = palette,
                                                         direction = direction,
                                                         na.value = naColor)
            else
                layer <- ggplot2::scale_fill_distiller(palette = palette,
                                                       direction = direction,
                                                       na.value = naColor)
    }
    # When nothing set, return NULL. "+ NULL" on a ggplot object doesn't
    # change anything
    return(layer)
}
