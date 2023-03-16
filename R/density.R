#' Create density plot basing on specified coordinates
#' @description This function shows the cell density presented in a 2D
#' dimensionality reduction coordinates. Density is shown with coloring and
#' contour lines. A scatter plot of the dimensionality reduction is added as
#' well. The density plot can be splitted by categorical variables (e.g.
#' \code{"dataset"}), while the scatter plot will always be shown for all cells
#' in subplots as a reference of the global structure.
#' @param object A \linkS4class{liger} object
#' @param useDimRed Name of the variable storing dimensionality reduction result
#' in the \code{cellMeta} slot. Default \code{"UMAP"}.
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
#' @param ... More theme setting arguments passed to
#' \code{\link{.ggplotLigerTheme}}.
#' @return A ggplot object when only one plot is generated, A ggplot object
#' combined with \code{\link[cowplot]{plot_grid}} when multiple plots and
#' \code{combinePlot = TRUE}. A list of ggplot when multiple plots and
#' \code{combinePlot = FALSE}.
#' @export
#' @examples
#' data("pbmcPlot", package = "rliger")
#' # Example dataset has small number of cells, thus cutoff adjusted.
#' plotDensityDimRed(pbmcPlot, minDensity = 1)
plotDensityDimRed <- function(
        object,
        useDimRed = "UMAP",
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
    xVar <- paste0(useDimRed, ".1")
    yVar <- paste0(useDimRed, ".2")
    dr <- .fetchCellMetaVar(object, c(xVar, yVar))
    splitVar <- .fetchCellMetaVar(object, splitBy, checkCategorical = TRUE,
                                  drop = FALSE)

    if (!is.null(splitVar) && ncol(splitVar) > 0) {
        # Will be splitting into subplots
        drList <- split(dr, splitVar)
        if (!is.null(title) && length(title) != length(drList)) {
            warning("Length of `title` does not match to number of categories ",
                    "identified from `splitBy`. Using default titles. ")
            title <- NULL
        }
    } else {
        # Will return a single ggplot
        if (length(title) > 1) {
            warning("`title` has length greater than 1 while only a single ",
                    "plot is generated. Using the first value only. ")
            title <- title[1]
        }
        drList <- list(dr)
    }
    plotList <- list()
    if (length(drList) == 0) stop("No plot could be generated")
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
