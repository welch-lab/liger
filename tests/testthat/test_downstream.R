data("pbmc", package = "rliger2")

withNewH5Copy <- function(fun) {
    ctrlpath.orig <- system.file("extdata/ctrl.h5", package = "rliger2")
    stimpath.orig <- system.file("extdata/stim.h5", package = "rliger2")
    if (!file.exists(ctrlpath.orig))
        stop("Cannot find original h5 file at: ", ctrlpath.orig)
    if (file.exists("ctrltest.h5")) file.remove("ctrltest.h5")
    if (file.exists("stimtest.h5")) file.remove("stimtest.h5")
    pwd <- getwd()
    # Temp setting for GitHub Actions
    fsep <- ifelse(Sys.info()["sysname"] == "Windows", "\\", "/")
    if (Sys.info()["sysname"] == "Windows") {
        pwd <- file.path("C:\\Users", Sys.info()["user"], "Documents", fsep = fsep)
    }

    ctrlpath <- file.path(pwd, "ctrltest.h5", fsep = fsep)
    stimpath <- file.path(pwd, "stimtest.h5", fsep = fsep)
    cat("Working ctrl H5 file path: ", ctrlpath, "\n")
    cat("Working stim H5 file path: ", stimpath, "\n")
    file.copy(ctrlpath.orig, ctrlpath, copy.mode = TRUE)
    file.copy(stimpath.orig, stimpath, copy.mode = TRUE)
    if (!file.exists(ctrlpath))
        stop("Cannot find copied h5 file at: ", ctrlpath)

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

process <- function(object, f = TRUE, q = TRUE) {
    object <- normalize(object)
    object <- selectGenes(object)
    object <- scaleNotCenter(object)
    if (f) object <- runOnlineINMF(object, k = 20, minibatchSize = 100)
    if (q) object <- quantileNorm(object)
    return(object)
}

is_online <- function() {
    tryCatch({
        readLines("https://cran.r-project.org/", n = 1)
        TRUE
    },
    warning = function(w) invokeRestart("muffleWarning"),
    error = function(e) FALSE)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Clustering
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Clustering")
test_that("clustering", {
    pbmc <- process(pbmc, f = FALSE, q = FALSE)
    expect_error(runCluster(pbmc), "No factor loading ")

    pbmc <- runOnlineINMF(pbmc, k = 20, minibatchSize = 100)
    expect_message(runCluster(pbmc, nRandomStarts = 1),
                   "leiden clustering on unnormalized")
    expect_message(runCluster(pbmc, nRandomStarts = 1, method = "louvain"),
                   "louvain clustering on unnormalized")

    pbmc <- quantileNorm(pbmc)
    expect_message(runCluster(pbmc, nRandomStarts = 1),
                   "leiden clustering on quantile normalized")
    expect_message(runCluster(pbmc, nRandomStarts = 1, method = "louvain"),
                   "louvain clustering on quantile normalized")

    # Tests for singleton grouping. Need to find the case where there are singletons
    # expect_message(runCluster(pbmc, nRandomStarts = 1,
    #                          partitionType = "CPMVertexPartition"),
    #                "63 singletons identified. 112 final clusters")
    # pbmc <- runCluster(pbmc, nRandomStarts = 1,
    #                          partitionType = "CPMVertexPartition",
    #                          groupSingletons = FALSE)
    # expect_equal(nlevels(pbmc$leiden_cluster), 113)
    # expect_true("singleton" %in% pbmc$leiden_cluster)

    # Test downsampling with the info already here
    pbmc.small1 <- downsample(pbmc, maxCells = 100)
    expect_true(all.equal(sapply(pbmc.small1@datasets, ncol), c(ctrl = 44, stim = 56)))
    idx <- downsample(pbmc, balance = c("dataset"), maxCells = 100,
                              returnIndex = TRUE)
    expect_is(idx, "integer")
    expect_equal(length(idx), 200)

    pbmc <- mapCellMeta(pbmc, from = "dataset", newTo = "pseudo",
                        `ctrl` = "CTRL")
    expect_identical(pbmc$pseudo,
                     setNames(factor(c(rep("CTRL", 300), rep("stim", 300))),
                              colnames(pbmc)))
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dimensionality reduction
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("dimensionality reduction")
test_that("dimensionality reduction", {
    pbmc <- process(pbmc)
    expect_message(runUMAP(pbmc, useRaw = TRUE),
                   "Generating UMAP on unnormalized")
    expect_message(pbmc <- runUMAP(pbmc, useRaw = FALSE),
                   "Generating UMAP on quantile normalized")
    expect_equal(dim(pbmc$UMAP), c(ncol(pbmc), 2))

    expect_message(runTSNE(pbmc, useRaw = TRUE),
                   "Generating TSNE \\(Rtsne\\) on unnormalized")
    expect_message(pbmc <- runTSNE(pbmc, useRaw = FALSE),
                   "Generating TSNE \\(Rtsne\\) on quantile normalized")
    expect_equal(dim(pbmc$TSNE), c(ncol(pbmc), 2))

    expect_error(runTSNE(pbmc, method = "fft"),
                 "Please pass in path to FIt-SNE directory as fitsne.path.")
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Differential Expression
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Differential Expression")
test_that("wilcoxon", {
    expect_error(runMarkerDEG(pbmc),
                 "No `conditionBy` given or default cluster not set")
    pbmc <- process(pbmc)
    pbmc <- runCluster(pbmc, nRandomStarts = 1)
    res0 <- runMarkerDEG(pbmc, conditionBy = "dataset", useDatasets = 1)
    expect_true(all(is.nan(res0$pval)))

    res1 <- runMarkerDEG(pbmc)
    expect_equal(dim(res1), c(249 * nlevels(pbmc$leiden_cluster), 10))
    res2 <- runMarkerDEG(pbmc, conditionBy = "dataset", splitBy = "leiden_cluster")
    expect_is(res2, "list")
    hm1 <- plotMarkerHeatmap(pbmc, res1, dedupBy = "l")
    hm2 <- plotMarkerHeatmap(pbmc, res1, dedupBy = "p")
    expect_is(hm1, "HeatmapList")
    expect_is(hm2, "HeatmapList")
    expect_is(plotVolcano(res1, 0), "ggplot")
    expect_is(plotEnhancedVolcano(res1, 0), "ggplot")

    expect_error(getFactorMarkers(pbmc, "ctrl", "stim", factorShareThresh = 0),
                 "No factor passed the dataset specificity threshold")
    expect_warning(
        expect_message(
            res3 <- getFactorMarkers(pbmc, "ctrl", "stim", printGenes = TRUE)
        )
    )
    expect_is(res3, "list")
    expect_identical(names(res3), c("ctrl", "shared", "stim", "num_factors_V1",
                                   "num_factors_V2"))

    expect_error(runGOEnrich(res1, group = "a"),
                 "Selected groups not available")
    expect_error(runGOEnrich(res1, group = 0, orderBy = c("logFC", "pval")),
                 "Only one `orderBy`")
    expect_error(runGOEnrich(res1, orderBy = "score"),
                 "`orderBy` should be one of")
    if (is_online()) {
        go1 <- runGOEnrich(res1, group = 0, orderBy = "logFC", significant = FALSE)
        expect_is(go1, "list")
        expect_is(go1$`0`$result, "data.frame")
        go2 <- runGOEnrich(res1, group = 0, orderBy = "pval", significant = FALSE)
        expect_is(go2, "list")
        expect_is(go2$`0`$result, "data.frame")
    }
})

# test_that("pseudo bulk - group wise", {
    # expect_error(runPseudoBulkDEG("hey"), "Please use a `liger` object.")
    #
    # rawData(datasets(pbmcPlot)[[1]]) <- rawData(dataset(pbmc, 1))
    # rawData(datasets(pbmcPlot)[[2]]) <- rawData(dataset(pbmc, 2))
    #
    # expect_error(runPseudoBulkDEG(pbmcPlot),
    #              "Either `groups` or `markerBy`")
    # expect_error(runPseudoBulkDEG(pbmcPlot, groups = 1, markerBy = "dataset"),
    #              "Only one of `groups` and `markerBy`")
    #
    # expect_error(runPseudoBulkDEG(pbmcPlot, groups = 1:10),
    #              "Please use a named list for `groups`")
    # expect_error(runPseudoBulkDEG(pbmcPlot,
    #                               groups = list(1)),
    #              "Please use at least 2 elements in `groups` list")
    # # Auto naming when un-named comparison group
    # expect_message(
    #     res <- runPseudoBulkDEG(
    #         pbmcPlot,
    #         groups = list(
    #             pbmcPlot$leiden_cluster == 1,
    #             pbmcPlot$leiden_cluster == 2
    #         )
    #     ),
    #     "Generating pseudo-bulks for condition \"group1\""
    # )
    # expect_is(res, "data.frame")
    #
    # res <- runPseudoBulkDEG(
    #     pbmcPlot,
    #     groups = list(
    #         pbmcPlot$leiden_cluster == 1,
    #         pbmcPlot$leiden_cluster == 2
    #     ),
    #     useCellMetaVar = "dataset"
    # )
    # expect_is(res, "data.frame")
    #
    # expect_error(
    #     res <- runPseudoBulkDEG(
    #         pbmcPlot,
    #         groups = list(
    #             c1 = pbmcPlot$leiden_cluster == 1,
    #             c2 = pbmcPlot$leiden_cluster == 2
    #         ),
    #         replicateAnn = list()
    #     ),
    #     "Please use a `data.frame` or a `factor` to specify replicate "
    # )
    #
    # ann <- cellMeta(pbmcPlot, "dataset", as.data.frame = TRUE, drop = FALSE)
    #
    # ann1 <- data.frame(v1 = ann$dataset)
    # expect_error(
    #     res <- runPseudoBulkDEG(
    #         pbmcPlot,
    #         groups = list(
    #             c1 = pbmcPlot$leiden_cluster == 1,
    #             c2 = pbmcPlot$leiden_cluster == 2
    #         ),
    #         replicateAnn = ann1
    #     ),
    #     "Not all cells involved in `groups` are annotated in "
    # )
    #
    # ann2 <- ann$dataset
    # expect_no_error(
    #     runPseudoBulkDEG(
    #         pbmcPlot,
    #         groups = list(
    #             c1 = pbmcPlot$leiden_cluster == 1,
    #             c2 = pbmcPlot$leiden_cluster == 2
    #         ),
    #         replicateAnn = ann2
    #     )
    # )
    #
    # expect_error(
    #     runPseudoBulkDEG(
    #         pbmcPlot,
    #         groups = list(
    #             c1 = pbmcPlot$leiden_cluster == 1,
    #             c2 = pbmcPlot$leiden_cluster == 2
    #         ),
    #         replicateAnn = ann2[pbmcPlot$leiden_cluster == 1]
    #     ),
    #     "Unable to format replicate annotation with given"
    # )
#
#     ann3 <- ann2
#     names(ann3) <- colnames(pbmcPlot)
#     expect_no_error(
#         runPseudoBulkDEG(
#             pbmcPlot,
#             groups = list(
#                 c1 = pbmcPlot$leiden_cluster == 1,
#                 c2 = pbmcPlot$leiden_cluster == 2
#             ),
#             replicateAnn = ann3[pbmcPlot$leiden_cluster %in% 1:3]
#         )
#     )
#     expect_error(
#         runPseudoBulkDEG(
#             pbmcPlot,
#             groups = list(
#                 c1 = pbmcPlot$leiden_cluster == 1,
#                 c2 = pbmcPlot$leiden_cluster == 2
#             ),
#             replicateAnn = ann3[pbmcPlot$leiden_cluster == 1]
#         ),
#         "Missing cells: ctrl_AAGGCTTGGTTCGA.1"
#     )
#
#     expect_warning(
#         runPseudoBulkDEG(
#             pbmcPlot,
#             groups = list(
#                 c1 = pbmcPlot$leiden_cluster == 1 & pbmcPlot$dataset == "ctrl",
#                 c2 = pbmcPlot$leiden_cluster == 0 & pbmcPlot$dataset == "stim"
#             ),
#             useCellMetaVar = "dataset"
#         ),
#         "will not create pseudo-bulks but test at single cell level"
#     )
#
#     res2 <- runPseudoBulkDEG(
#         pbmcPlot, markerBy = "leiden_cluster"
#     )
#     expect_equal(ncol(res2), 5)
# })

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GSEA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("GSEA")
custom <- list(
    `Immune System` = c("9636", "2633", "6282", "6280", "6279", "2207", "2214",
                        "6402", "91543", "6233", "10578", "3553", "5473",
                        "3627", "51316", "929", "972")
)

if (is_online()) {
    test_that("gsea", {
        expect_warning({
            expect_is(runGSEA(pbmcPlot, genesets = "Immune System"), "list")
            expect_is(runGSEA(pbmcPlot, customGenesets = custom), "list")
        })
    })
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ATAC
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("ATAC")
data("bmmc")
test_that("ATAC", {
    bmmc <- normalize(bmmc)
    bmmc <- selectGenes(bmmc)
    bmmc <- scaleNotCenter(bmmc)
    bmmc <- runOnlineINMF(bmmc, minibatchSize = 80)
    bmmc <- quantileNorm(bmmc)
    bmmc <- normalizePeak(bmmc)
    bmmc <- imputeKNN(bmmc, reference = "atac", queries = "rna")
    expect_is(dataset(bmmc, "rna"), "ligerATACDataset")
    corr <- linkGenesAndPeaks(
        bmmc, useDataset = "rna",
        pathToCoords = system.file("extdata/hg19_genes.bed",
                                   package = "rliger2")
    )
    expect_is(corr, "dgCMatrix")
})

