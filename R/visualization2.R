#' @export
plotByDatasetAndCluster <- function(
        object,
        useDimRed = c("UMAP", "TSNE"),
        cluster = "louvain_cluster",
        combinePlots = TRUE,
        labelTextDataset = FALSE,
        labelTextCluster = TRUE,
        ...
) {
    useDimRed <- match.arg(useDimRed)
    xVar <- paste0(useDimRed, ".1")
    yVar <- paste0(useDimRed, ".2")

    p1 <- plotCellScatter(object, x = xVar, y = yVar, colorBy = "dataset",
                          labelText = labelTextDataset, ...)
    p2 <- plotCellScatter(object, x = xVar, y = yVar, colorBy = cluster,
                          labelText = labelTextCluster,...)
    plot <- list(dataset = p1, cluster = p2)
    if (isTRUE(combinePlots)) {
        plot <- cowplot::plot_grid(plotlist = plot, nrow = 1,
                                   align = "h", axis = "tblr")
    }
    return(plot)
}

plotGene <- function(
        object,
        gene,
        ...
) {

}

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
        yFun = function(x) log2(10000*x + 1),
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
