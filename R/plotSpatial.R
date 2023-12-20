#' Visualize a spatial dataset
#' @export
#' @rdname plotSpatial
#' @param object Either a \linkS4class{liger} object containing a spatial
#' dataset or a \linkS4class{ligerSpatialDataset} object.
#' @param ... Arguments passed to other methods. \code{.liger} method passes
#' everything to \code{.ligerSpatialDataset} method, and the latter passes
#' everything to \code{\link{.ggScatter}} and then
#' \code{\link{.ggplotLigerTheme}}.
#' @return A ggplot object
plotSpatial2D <- function(object, ...) {
    UseMethod("plotSpatial2D", object)
}

#' @export
#' @rdname plotSpatial
#' @method plotSpatial2D liger
#' @param dataset Name of one spatial dataset.
#' @param useCluster Either the name of one variable in \code{cellMeta(object)}
#' or a factor object with annotation that matches with all cells in the
#' specified dataset. Default \code{NULL} does not add any annotation.
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
    .checkArgLen(dataset, 1)
    ld <- dataset(object, dataset)
    useCluster <- useCluster %||%
        defaultCluster(object)[object$dataset == dataset]
    if (length(useCluster) == 1) {
        legendColorTitle <- legendColorTitle %||% useCluster
        useCluster <- cellMeta(object, useCluster, useDataset = dataset)
    } else {
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
plotSpatial2D.ligerSpatialDataset <- function(
        object,
        useCluster = NULL,
        legendColorTitle = NULL,
        useDims = c(1, 2),
        xlab = NULL,
        ylab = NULL,
        labelText = FALSE,
        ...)
{
    .checkArgLen(useCluster, ncol(object))
    legendColorTitle <- legendColorTitle %||% "Annotation"

    coord <- coordinate(object)
    .checkArgLen(useDims, 2)
    coord <- coord[, useDims]
    plotDF <- as.data.frame(coord)
    colnames(plotDF) <- c("x", "y")
    xRange <- c(min(plotDF$x), max(plotDF$x))
    yRange <- c(min(plotDF$y), max(plotDF$y))

    plotDF[[legendColorTitle]] <- factor(useCluster)
    .ggScatter(plotDF = plotDF, x = "x", y = "y", colorBy = legendColorTitle,
               xlab = xlab, ylab = ylab, labelText = labelText,
               legendColorTitle = legendColorTitle, ...) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.text = ggplot2::element_blank()) +
        ggplot2::coord_fixed(xlim = xRange, ylim = yRange)
}
