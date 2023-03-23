data("pbmcPlot", package = "rliger2")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Scatter plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Scatter plots")
test_that("scatter plots", {
    # Wrappers
    expect_is(plotClusterDimRed(pbmcPlot), "ggplot")
    expect_is(plotDatasetDimRed(pbmcPlot), "ggplot")
    expect_is(plotByDatasetAndCluster(pbmcPlot), "ggplot")
    expect_is(plotGeneDimRed(pbmcPlot, "S100A8"), "ggplot")
    expect_is(plotFactorDimRed(pbmcPlot, 1), "ggplot")

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
        expect_is(plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2",
                                  colorBy = "dataset", dotOrder = o), "ggplot")
    }
    expect_is(plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2",
                              colorBy = "S100A8", slot = "normData",
                              trimHigh = 5, trimLow = 0), "ggplot")

    expect_is(
        plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2",
                        colorBy = "leiden_cluster", shapeBy = "dataset"),
        "ggplot"
    )
    expect_is(
        plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2",
                        colorBy = NULL, shapeBy = "dataset"),
        "ggplot"
    )
    expect_is(
        plotCellScatter(pbmcPlot, "UMAP.1", "UMAP.2",
                        colorBy = "leiden_cluster", raster = TRUE),
        "ggplot"
    )
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Violin plots
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Violin plots")
test_that("Violin plots", {

})
