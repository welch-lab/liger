data("pbmcPlot", package = "rliger2")
expect_gg <- function(...) {
    objects <- list(...)
    for (i in seq_along(objects)) {
        testthat::expect_is(objects[[i]], "ggplot")
    }

}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Scatter plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Scatter plots")
test_that("scatter plots", {
    # Wrappers
    expect_gg(
        plotClusterDimRed(pbmcPlot),
        plotDatasetDimRed(pbmcPlot),
        plotByDatasetAndCluster(pbmcPlot),
        plotGeneDimRed(pbmcPlot, "S100A8"),
        plotFactorDimRed(pbmcPlot, 1)
    )

    # General
    expect_is(
        plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2", splitBy = "dataset"),
        "list"
    )
    expect_is(
        plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2", colorBy = "dataset",
                        splitBy = "dataset"),
        "list"
    )
    for (o in c("shuffle", "ascending", "descending")) {
        expect_gg(plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2",
                                  colorBy = "dataset", dotOrder = o))
    }
    expect_gg(
        plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2", colorBy = "S100A8",
                        slot = "normData", trimHigh = 5, trimLow = 0),
        plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2",
                        colorBy = "leiden_cluster", shapeBy = "dataset"),
        plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2", colorBy = NULL,
                        shapeBy = "dataset"),
        plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2",
                        colorBy = "leiden_cluster", raster = TRUE)
    )
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Violin plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Violin plots")
test_that("Violin plots", {
    # Wrappers
    expect_is(plotGeneViolin(pbmcPlot, "S100A8"), "list")
    expect_is(plotGeneViolin(pbmcPlot, "S100A8", groupBy = FALSE, dot = TRUE,
                             dotColor = NULL), "list")
    expect_gg(
        plotGeneViolin(pbmcPlot, "S100A8", byDataset = FALSE),
        plotTotalCountViolin(pbmc, dot = TRUE),
        plotGeneDetectedViolin(pbmc, dot = TRUE, dotColor = NULL,
                               raster = TRUE)
    )

    # General
    expect_gg(plotCellViolin(pbmcPlot, "nUMI", groupBy = NULL,
                             colorBy = "leiden_cluster", box = TRUE, dot = TRUE,
                             raster = TRUE))
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ggplot themes
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("ggplot themes")
test_that("ggplot themes", {
    nCluster <- length(levels(pbmcPlot$leiden_cluster))
    expect_gg(
        plotClusterDimRed(pbmcPlot,
                          titleSize = 10,
                          subtitleSize = 8,
                          xTextSize = 8,
                          xFacetSize = 8,
                          xTitleSize = 10,
                          yTextSize = 8,
                          yFacetSize = 8,
                          yTitleSize = 10,
                          legendTextSize = 8,
                          legendTitleSize = 10,
                          xlab = NULL, ylab = NULL,
                          panelBorder = TRUE,
                          colorLabels = letters[1:nCluster],
                          colorValues = RColorBrewer::brewer.pal(nCluster,
                                                                 "Set1"),
                          showLegend = FALSE),
        plotTotalCountViolin(pbmcPlot, colorBy = "dataset",
                             colorLabel = letters[1:2],
                             colorValues = c("green", "purple"),
                             legendFillTitle = "hi", legendColorTitle = "hi"),
        plotGeneDimRed(pbmcPlot, "S100A8", colorLow = "red",
                       colorMid = "white",
                       colorHigh = "blue",
                       colorMidPoint = 5),
        plotGeneDimRed(pbmcPlot, "S100A8", colorPalette = "Pastel1")
    )
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Density plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Density plot")
test_that("Density plot", {
    expect_gg(
        expect_warning(plotDensityDimRed(pbmcPlot, splitBy = "dataset",
                                         title = "one"),
                       "Length of `title` does"),
        expect_warning(plotDensityDimRed(pbmcPlot, title = letters[1:3],
                                         dotRaster = TRUE),
                       "`title` has length greater than")
    )
    expect_is(plotDensityDimRed(pbmcPlot, "UMAP", splitBy = "dataset",
                                title = names(pbmcPlot), combinePlot = FALSE),
              "list")
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Proportion plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Proportion plots")
test_that("Proportion plots", {
    expect_gg(
        plotProportion(pbmcPlot),
        plotProportion(pbmcPlot, method = "pie"),
        plotProportionBar(pbmcPlot, method = "group"),
        plotClusterProportions(pbmcPlot)
    )
    expect_is(plotProportionBar(pbmcPlot, inclRev = TRUE, combinePlot = FALSE),
              "list")
    expect_error(plotProportionDot(pbmcPlot, letters),
                 "`class1` and `class2` must be")
    expect_error(plotProportionBar(pbmcPlot, letters),
                 "`class1` and `class2` must be")
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Heatmap
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Heatmap")
test_that("Heatmap", {
    expect_is(
        plotFactorHeatmap(pbmcPlot),
        "HeatmapList"
    )
    expect_error(
        plotFactorHeatmap(pbmcPlot, viridisOption = letters),
        "`viridisOption` has to be one value from the available choices"
    )
    expect_is(
        plotGeneHeatmap(pbmcPlot, varFeatures(pbmcPlot),
                        useCellMeta = "leiden_cluster", scale = TRUE),
        "HeatmapList"
    )
    extAnn <- cellMeta(pbmcPlot, c("nUMI", "dataset"), as.data.frame = TRUE)
    expect_is(
        plotGeneHeatmap(pbmcPlot, varFeatures(pbmcPlot),
                        cellAnnotation = extAnn, transpose = TRUE,
                        cellTextSize = 8,
                        featureTextSize = 8,
                        cellTitleSize = 10,
                        featureTitleSize = 10,
                        legendTextSize = 8,
                        legendTitleSize = 10),
        "HeatmapList"
    )
    expect_warning(
        plotGeneHeatmap(pbmcPlot, varFeatures(pbmcPlot), scale = TRUE,
                        trim = c(-5, 0, 5)),
        "`trim` must be a numeric vector of two values"
    )
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dot Plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Dot Plot")
test_that("Dot Plot", {
    expect_error(plotClusterGeneDot("hello"), "Please use a liger object")
    expect_error(plotClusterFactorDot(pbmcPlot, viridisOption = letters),
                 "`viridisOption` has to be one")
    expect_is(plotClusterGeneDot(pbmcPlot, varFeatures(pbmcPlot)[1:5],
                                 cellTextSize = 8,
                                 featureTextSize = 8,
                                 cellTitleSize = 10,
                                 featureTitleSize = 10,
                                 legendTextSize = 8,
                                 legendTitleSize = 10),
              "HeatmapList")
    expect_is(
        plotClusterGeneDot(pbmcPlot,
                           features = data.frame(varFeatures(pbmcPlot)[1:10],
                                                 factor(rep(letters[1:5], 2))),
                           transpose = TRUE),
        "HeatmapList"
    )

    expect_error(plotClusterFactorDot("hello"), "Please use a liger object")
    expect_is(plotClusterFactorDot(pbmcPlot, factorScaleFunc = function(x) x),
              "HeatmapList")
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Gene loading
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Gene loading")
test_that("Gene loading", {
    res <- getFactorMarkers(pbmcPlot, "ctrl", "stim")
    expect_gg(plotGeneLoadings(pbmcPlot, res, 1))
})

