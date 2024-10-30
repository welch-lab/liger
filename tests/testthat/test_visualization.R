data("pbmcPlot", package = "rliger")

withNewH5Copy <- function(fun) {
    ctrlpath.orig <- system.file("extdata/ctrl.h5", package = "rliger")
    stimpath.orig <- system.file("extdata/stim.h5", package = "rliger")
    if (!file.exists(ctrlpath.orig))
        stop("Cannot find original h5 file at: ", ctrlpath.orig)
    # if (file.exists("ctrltest.h5")) file.remove("ctrltest.h5")
    # if (file.exists("stimtest.h5")) file.remove("stimtest.h5")
    # pwd <- getwd()
    # # Temp setting for GitHub Actions
    # fsep <- ifelse(Sys.info()["sysname"] == "Windows", "\\", "/")
    # if (Sys.info()["sysname"] == "Windows") {
    #     pwd <- file.path("C:\\Users", Sys.info()["user"], "Documents", fsep = fsep)
    # }

    # ctrlpath <- file.path(pwd, "ctrltest.h5", fsep = fsep)
    # stimpath <- file.path(pwd, "stimtest.h5", fsep = fsep)
    ctrlpath <- tempfile(pattern = "ctrltest_", fileext = ".h5")
    stimpath <- tempfile(pattern = "stimtest_", fileext = ".h5")
    cat("Working ctrl H5 file path: ", ctrlpath, "\n")
    cat("Working stim H5 file path: ", stimpath, "\n")
    file.copy(ctrlpath.orig, ctrlpath, copy.mode = TRUE)
    file.copy(stimpath.orig, stimpath, copy.mode = TRUE)
    if (!file.exists(ctrlpath))
        stop("Cannot find copied h5 file at: ", ctrlpath)
    if (!file.exists(stimpath))
        stop("Cannot find copied h5 file at: ", stimpath)

    fun(list(ctrl = ctrlpath, stim = stimpath))

    if (file.exists(ctrlpath)) unlink(ctrlpath)
    if (file.exists(stimpath)) unlink(stimpath)
}

closeH5Liger <- function(object) {
    for (d in names(object)) {
        if (isH5Liger(object, d)) {
            h5file <- getH5File(object, d)
            h5file$close()
        }
    }
}

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
        plotDimRed(pbmcPlot, splitBy = "dataset"),
        "list"
    )
    expect_is(
        plotDimRed(pbmcPlot, colorBy = "dataset",
                        splitBy = "dataset"),
        "list"
    )
    for (o in c("shuffle", "ascending", "descending")) {
        expect_gg(plotDimRed(pbmcPlot, colorBy = "dataset", dotOrder = o))
    }
    expect_gg(
        plotDimRed(pbmcPlot, colorBy = "S100A8", slot = "normData",
                   trimHigh = 5, trimLow = 0),
        plotDimRed(pbmcPlot, colorBy = "leiden_cluster", shapeBy = "dataset"),
        plotDimRed(pbmcPlot, colorBy = NULL, shapeBy = "dataset")
    )
    skip_if_not_installed("scattermore")
    expect_gg(
        plotDimRed(pbmcPlot, colorBy = "leiden_cluster", raster = TRUE)
    )

    expect_gg(plotGroupClusterDimRed(pbmcPlot))
    do.call(expect_gg, plotGroupClusterDimRed(pbmcPlot, combinePlot = FALSE))

    do.call(expect_gg, plotBarcodeRank(pbmc))
    # Fake operation to create ATAC datasets
    pbmcPlot@datasets$ctrl <- as.ligerDataset(dataset(pbmcPlot, "ctrl"), "atac")
    pbmcPlot@datasets$stim <- as.ligerDataset(dataset(pbmcPlot, "stim"), "atac")
    normPeak(pbmcPlot, "ctrl") <- normData(pbmcPlot, "ctrl")
    normPeak(pbmcPlot, "stim") <- normData(pbmcPlot, "stim")
    expect_gg(plotPeakDimRed(pbmcPlot, "ISG15"))

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
        plotTotalCountViolin(pbmc, dot = TRUE)
    )

    expect_gg(
        plotClusterGeneViolin(pbmcPlot, "S100A8", box = TRUE, colorBy = "dataset"),
        plotClusterGeneViolin(pbmcPlot, "S100A8", groupBy = FALSE)
    )
    # General
    skip_if_not_installed("scattermore")
    expect_gg(
        plotGeneDetectedViolin(pbmc, dot = TRUE, dotColor = NULL, raster = TRUE),
        plotCellViolin(pbmcPlot, "nUMI", groupBy = NULL,
                       colorBy = "leiden_cluster", box = TRUE, dot = TRUE,
                       raster = TRUE)
    )
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
        expect_no_warning(plotDensityDimRed(pbmcPlot, splitBy = "dataset",
                                            title = "one"))
    )
    expect_is(plotDensityDimRed(pbmcPlot, "UMAP", splitBy = "dataset",
                                title = names(pbmcPlot), combinePlot = FALSE),
              "list")
    skip_if_not(requireNamespace("scattermore", quietly = TRUE))
    expect_gg(
        expect_message(plotDensityDimRed(pbmcPlot, title = letters[1:3],
                                         dotRaster = TRUE),
                       "`title` has length greater than")
    )
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
        plotProportionDot(pbmcPlot)
    )
    expect_is(plotProportionBar(pbmcPlot, inclRev = TRUE, combinePlot = FALSE),
              "list")
    expect_error(plotProportionDot(pbmcPlot, letters),
                 "`class1` and `class2` must be")
    expect_error(plotProportionBar(pbmcPlot, letters),
                 "`class1` and `class2` must be")

    defaultCluster(pbmcPlot) <- NULL
    expect_error(plotProportionBox(pbmcPlot), "No cluster specified nor default set")
    defaultCluster(pbmcPlot) <- "leiden_cluster"
    expect_error(plotProportionBox(pbmcPlot, conditionBy = "leiden_cluster"),
                 "Condition variable must be a high level variable of the datasets")
    expect_gg(
        plotProportionBox(pbmcPlot, dot = TRUE),
        plotProportionBox(pbmcPlot, conditionBy = "dataset")
    )
    do.call(expect_gg, plotProportionBox(pbmcPlot, splitByCluster = TRUE, dot = TRUE))
    do.call(expect_gg, plotProportionBox(pbmcPlot, splitByCluster = TRUE, conditionBy = "dataset"))

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

context("spatial coordinates")
test_that("Plot spatial coordinates", {
    ctrl.fake.spatial <- as.ligerDataset(dataset(pbmc, "ctrl"), modal = "spatial")
    fake.coords <- matrix(rnorm(2 * ncol(ctrl.fake.spatial)), ncol = 2)
    dimnames(fake.coords) <- list(colnames(ctrl.fake.spatial), c("x", "y"))
    coordinate(ctrl.fake.spatial) <- fake.coords
    dataset(pbmc, "ctrl") <- ctrl.fake.spatial
    expect_gg(plotSpatial2D(pbmc, dataset = "ctrl"))
    expect_gg(plotSpatial2D(pbmc, dataset = "ctrl", useCluster = "dataset"))
})

context("Sankey")
test_that("PlotSankey", {
    skip_if_not_installed("sankey")
    cellMeta(pbmcPlot, "ctrl_cluster", "ctrl") <-
        cellMeta(pbmcPlot, "leiden_cluster", "ctrl")
    cellMeta(pbmcPlot, "stim_cluster", "stim") <-
        cellMeta(pbmcPlot, "leiden_cluster", "stim")
    pdfName <- tempfile(pattern = "fig_", fileext = ".pdf")
    grDevices::pdf(file = pdfName)
    expect_no_error({
        plotSankey(pbmcPlot, "ctrl_cluster", "stim_cluster",
                   titles = c("control", "LIGER", "stim"),
                   prefixes = c("c", NA, "s"))
    })
    grDevices::dev.off()
    if (file.exists(pdfName)) unlink(pdfName)
})
