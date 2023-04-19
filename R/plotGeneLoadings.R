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
#' @param ... Additional plot theme setting arguments passed to
#' \code{\link{.ggScatter}} and \code{\link{.ggplotLigerTheme}}.
#' @export
#' @examples
#' result <- getFactorMarkers(pbmcPlot, "ctrl", "stim")
#' plotGeneLoadings(pbmcPlot, result, useFactor = 2)
plotGeneLoadings <- function(
        object,
        markerTable,
        useFactor,
        useDimRed = "UMAP",
        nLabel = 15,
        nPlot = 30,
        ...
) {
    p1 <- plotCellScatter(
        object,
        x = paste0(useDimRed, ".1"), y = paste0(useDimRed, ".2"),
        colorBy = useFactor, slot = "H.norm", zeroAsNA = TRUE, dotOrder = "asc",
        splitBy = NULL, shapeBy = NULL, colorByFunc = NULL, ...
    )
    bottom <- plotGeneLoadingRank(object, markerTable, useFactor,
                                  nLabel, nPlot, ...)
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
