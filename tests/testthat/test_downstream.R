data("pbmc", package = "rliger")

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
    skip_if_not_installed("RcppPlanc")
    pbmc <- process(pbmc, f = FALSE, q = FALSE)
    expect_error(runCluster(pbmc), "No cell factor loading available")

    pbmc <- runOnlineINMF(pbmc, k = 20, minibatchSize = 100)
    expect_message(runCluster(pbmc, nRandomStarts = 1),
                   "leiden clustering on unaligned")

    expect_message(runCluster(pbmc, nRandomStarts = 1, method = "louvain"),
                   "louvain clustering on unaligned")

    pbmc <- quantileNorm(pbmc)
    expect_message(pbmc <- runCluster(pbmc, nRandomStarts = 1, saveSNN = TRUE),
                   "leiden clustering on aligned")
    expect_is(defaultCluster(pbmc, droplevels = TRUE), "factor")
    expect_is(pbmc@uns$snn, "dgCMatrix")
    expect_message(pbmc <- runCluster(pbmc, nRandomStarts = 1, method = "louvain"),
                   "louvain clustering on aligned")
    expect_message(defaultCluster(pbmc, name = "louvain_cluster") <- "louvain_cluster",
                   "Cannot have")
    expect_error(defaultCluster(pbmc) <- "notexist", "Selected variable does not exist")
    defaultCluster(pbmc) <- pbmc$leiden_cluster
    expect_identical(pbmc$leiden_cluster, pbmc$defaultCluster)
    expect_error(defaultCluster(pbmc) <- factor(letters), "Length of")
    defaultCluster(pbmc) <- NULL
    defaultCluster(pbmc, name = "leiden") <- unname(pbmc$leiden_cluster)
    expect_identical(pbmc$leiden, pbmc$leiden_cluster)

    fakevar <- pbmc$leiden_cluster
    names(fakevar)[1:26] <- letters
    expect_error(defaultCluster(pbmc) <- fakevar, "Not all `names")



    expect_equal(calcPurity(pbmc, "leiden_cluster", "leiden_cluster"), 1)
    expect_error(calcPurity(pbmc, letters, "leiden_cluster"),
                 "Longer/shorter `trueCluster` than cells considered requires")
    expect_message(calcPurity(pbmc, unname(pbmc$leiden_cluster), "leiden_cluster"), "Assuming unnamed")

    expect_equal(calcARI(pbmc, "leiden_cluster", "leiden_cluster"), 1)
    expect_error(calcARI(pbmc, letters, "leiden_cluster"),
                 "Longer/shorter `trueCluster` than cells considered requires")
    expect_message(calcARI(pbmc, unname(pbmc$leiden_cluster), "leiden_cluster"), "Assuming unnamed")

    expect_equal(calcNMI(pbmc, "leiden_cluster", "leiden_cluster"), 1)
    expect_error(calcNMI(pbmc, letters, "leiden_cluster"),
                 "Longer/shorter `trueCluster` than cells considered requires")
    expect_message(calcNMI(pbmc, unname(pbmc$leiden_cluster), "leiden_cluster"), "Assuming unnamed")
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
    skip_if_not_installed("RcppPlanc")
    pbmc <- process(pbmc)
    expect_message(runUMAP(pbmc, useRaw = TRUE),
                   "Generating UMAP on unaligned")
    expect_error(dimRed(pbmc), "available in this")
    expect_message(pbmc <- runUMAP(pbmc, useRaw = FALSE),
                   "Generating UMAP on aligned")
    pbmc@uns$defaultDimRed <- NULL
    expect_message(dimRed(pbmc), "No default")
    defaultDimRed(pbmc) <- "UMAP"
    expect_error(defaultDimRed(pbmc) <- letters, "Can only set one")
    expect_identical(defaultDimRed(pbmc), dimRed(pbmc, "UMAP"))
    expect_equal(dim(dimRed(pbmc)), c(ncol(pbmc), 2))
    expect_no_error(dimRed(pbmc, "UMAP2") <- dimRed(pbmc, "UMAP"))
    expect_equal(nrow(dimRed(pbmc, name = 1, cellIdx = 1:10)), 10)
    expect_equal(nrow(dimRed(pbmc, name = 1, useDatasets = names(pbmc))), ncol(pbmc))
    expect_equal(nrow(dimRed(pbmc, name = "UMAP", cellIdx = 1:10)), 10)
    expect_equal(nrow(dimRed(pbmc, name = "UMAP", useDatasets = names(pbmc))), ncol(pbmc))
    expect_no_error(dimRed(pbmc, 2) <- NULL)

    expect_message(runTSNE(pbmc, useRaw = TRUE),
                   "Generating TSNE \\(Rtsne\\) on unaligned")
    expect_message(pbmc <- runTSNE(pbmc, useRaw = FALSE),
                   "Generating TSNE \\(Rtsne\\) on aligned")
    expect_equal(dim(dimRed(pbmc, "TSNE")), c(ncol(pbmc), 2))

    expect_error(runTSNE(pbmc, method = "fft"),
                 "Please pass in path to FIt-SNE directory as fitsne.path.")
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Differential Expression
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Differential Expression")
test_that("wilcoxon", {
    skip_if_not_installed("RcppPlanc")
    expect_error(runMarkerDEG(pbmc),
                 "No `conditionBy` given or default cluster not set")
    pbmc <- process(pbmc)
    pbmc <- runCluster(pbmc, nRandomStarts = 1)
    res0 <- runMarkerDEG(pbmc, conditionBy = "dataset", useDatasets = 1, method = "wilcox")
    expect_true(all(is.nan(res0$pval)))

    res1 <- runMarkerDEG(pbmc, method = "wilcox")
    expect_equal(dim(res1), c(249 * nlevels(pbmc$leiden_cluster), 10))
    res2 <- runMarkerDEG(pbmc, conditionBy = "dataset", splitBy = "leiden_cluster", method = "wilcox")
    expect_is(res2, "data.frame")
    hm1 <- plotMarkerHeatmap(pbmc, res1, dedupBy = "l")
    hm2 <- plotMarkerHeatmap(pbmc, res1, dedupBy = "p")
    expect_is(hm1, "HeatmapList")
    expect_is(hm2, "HeatmapList")
    expect_is(plotVolcano(res1, 0), "ggplot")
    if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
        expect_is(plotEnhancedVolcano(res1, 0), "ggplot")
    }

    expect_error(getFactorMarkers(pbmc, "ctrl", "stim", factorShareThresh = 0),
                 "No factor passed the dataset specificity threshold")
    expect_message(
        res3 <- getFactorMarkers(pbmc, "ctrl", "stim", printGenes = TRUE)
    )
    expect_is(res3, "list")
    expect_identical(names(res3), c("ctrl", "shared", "stim", "num_factors_V1",
                                   "num_factors_V2"))

    skip_if_not_installed("gprofiler2")
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
        go3 <- runGOEnrich(res1, group = c(0, 1), orderBy = "pval", significant = FALSE)

        expect_is(plotGODot(go1, pvalThresh = 1), "ggplot")
        expect_error(plotGODot(go1, group = "ctrl"), "Specified group not available")
        expect_message(plotGODot(go1, group = '0'), "No enough matching")
        expect_is(plotGODot(go3, pvalThresh = 1), "list")
    }
})


test_that("pseudo bulk", {
    skip_if_not_installed("RcppPlanc")
    skip_if_not_installed("DESeq2")
    nIsecGenes <- length(Reduce(intersect, lapply(rawData(pbmc), rownames)))
    pbmc <- process(pbmc)
    pbmc <- runCluster(pbmc, nRandomStarts = 1)
    res1 <- runPairwiseDEG(pbmc, groupTest = pbmc$leiden_cluster == 1,
                           groupCtrl = pbmc$leiden_cluster == 2,
                           method = "pseudo")
    expect_is(res1, "data.frame")
    expect_true(all.equal(dim(res1), c(nIsecGenes, 7)))
    res2 <- runPairwiseDEG(pbmc, groupTest = 1, groupCtrl = 2,
                           variable1 = "leiden_cluster",
                           method = "pseudo", useReplicate = "dataset")
    expect_is(res2, "data.frame")
    expect_true(all.equal(dim(res2), c(nIsecGenes, 7)))
    res3 <- runPairwiseDEG(pbmc, groupTest = 1, groupCtrl = 2,
                           variable1 = "leiden_cluster",
                           method = "pseudo")
    expect_true(all.equal(res1[,-2], res3[,-2])) # Different in "group" column
    pbmc$leiden2 <- pbmc$leiden_cluster
    res4 <- runPairwiseDEG(pbmc, groupTest = 1, groupCtrl = 2,
                           variable1 = "leiden_cluster", variable2 = "leiden2",
                           method = "pseudo", useReplicate = "dataset")
    expect_true(all.equal(res2[,-2], res4[,-2])) # Different in "group" column

    expect_error(runPairwiseDEG(pbmc, variable2 = "yo"),
                 "Please see")

    expect_message(
        runPairwiseDEG(
            pbmc, groupTest = pbmc$dataset == "ctrl" & pbmc$leiden_cluster == 0,
            groupCtrl = pbmc$dataset == "stim" & pbmc$leiden_cluster == 0,
            method = "pseudo", useReplicate = "dataset", nPsdRep = 1
        ),
        "Too few replicates"
    )

    pbmc@datasets$ctrl@rawData <- NULL
    expect_error(
        runPairwiseDEG(
            pbmc, groupTest = 1, groupCtrl = 2,
            variable1 = "leiden_cluster",　method = "pseudo",
            useReplicate = "dataset"
        ),
        "not all available for involved datasets"
    )
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# GSEA
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("GSEA")
test_that("gsea", {
    skip_if_not(is_online())
    skip_if_not_installed("org.Hs.eg.db")
    skip_if_not_installed("reactome.db")
    skip_if_not_installed("fgsea")
    skip_if_not_installed("AnnotationDbi")
    custom <- list(
        `Immune System` = c("9636", "2633", "6282", "6280", "6279", "2207", "2214",
                            "6402", "91543", "6233", "10578", "3553", "5473",
                            "3627", "51316", "929", "972")
    )
    expect_warning({
        expect_is(runGSEA(pbmcPlot, genesets = "Immune System"), "list")
        expect_is(runGSEA(pbmcPlot, customGenesets = custom), "list")
    })
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ATAC
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("ATAC")
data("bmmc")
test_that("ATAC", {
    skip_if_not_installed("RcppPlanc")
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")
    skip_if_not_installed("psych")
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
                                   package = "rliger")
    )
    expect_is(corr, "dgCMatrix")
})

