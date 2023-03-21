data("pbmc", package = "rliger2")

process <- function(object, f = TRUE, q = TRUE) {
    object <- normalize(object)
    object <- selectGenes(object)
    object <- scaleNotCenter(object)
    if (f) object <- online_iNMF(object, k = 20, miniBatch_size = 100)
    if (q) object <- quantileNorm(object)
    return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Clustering
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Clustering")
test_that("clustering", {
    pbmc <- process(pbmc, f = FALSE, q = FALSE)
    expect_error(runLeidenCluster(pbmc), "No factor loading ")
    expect_error(runLouvainCluster(pbmc), "No factor loading ")

    pbmc <- online_iNMF(pbmc, k = 20, miniBatch_size = 100)
    expect_message(runLeidenCluster(pbmc, nRandomStarts = 1),
                   "Leiden clustering on unnormalized")
    expect_message(runLouvainCluster(pbmc, nRandomStarts = 1),
                   "Louvain clustering on unnormalized")

    pbmc <- quantileNorm(pbmc)
    expect_message(runLeidenCluster(pbmc, nRandomStarts = 1),
                   "Leiden clustering on quantile normalized")
    expect_message(runLouvainCluster(pbmc, nRandomStarts = 1),
                   "Louvain clustering on quantile normalized")

    expect_message(runLeidenCluster(pbmc, nRandomStarts = 1,
                             partitionType = "CPMVertexPartition"),
                   "63 singletons identified. 112 final clusters")
    pbmc <- runLeidenCluster(pbmc, nRandomStarts = 1,
                             partitionType = "CPMVertexPartition",
                             groupSingletons = FALSE)
    expect_equal(nlevels(pbmc$leiden_cluster), 113)
    expect_true("singleton" %in% pbmc$leiden_cluster)

    expect_warning(louvainCluster(pbmc), "was deprecated in")
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
