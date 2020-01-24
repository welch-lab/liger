## Tests for iNMF and quantile alignment

set.seed(1)

# pbmc.file <- system.file('tests', 'testdata', 'small_pbmc_data.RDS', package = 'liger')
pbmc.file <- "../testdata/small_pbmc_data.RDS"
pbmc.small <- readRDS(pbmc.file)

# preprocessing steps (object required for tests)
ligex <- createLiger(raw.data = pbmc.small, make.sparse = T, take.gene.union = F,
                     remove.missing = T)
ligex <- normalize(ligex)
ligex <- selectGenes(ligex, var.thresh = c(0.3, 0.9), do.plot = F)
ligex <- scaleNotCenter(ligex)

# Tests for iNMF factorization 
####################################################################################
context("iNMF factorization")

ligex <- optimizeALS(ligex, k = 15, lambda = 5, rand.seed = 1)

 test_that("Dataset names passed correctly", {
  expect_identical(names(ligex@H), c('tenx', 'seqwell'))
  expect_identical(names(ligex@V), c('tenx', 'seqwell'))
  expect_identical(rownames(ligex@scale.data[[1]]), rownames(ligex@H[[1]]))
  expect_identical(colnames(ligex@scale.data[[1]]), colnames(ligex@W))
})

test_that("Dimensions are correct", {
  expect_equal(dim(ligex@H[[1]]), c(250, 15))
  expect_equal(dim(ligex@H[[2]]), c(244, 15))
  expect_equal(dim(ligex@V[[1]]), c(15, 1907))
  expect_equal(dim(ligex@V[[2]]), c(15, 1907))
  expect_equal(dim(ligex@W), c(15, 1907))
})

test_that("Factorization is correct", {
  expect_equal(ligex@H[[1]][6, 1:5], c(3.762285e-05, 0, 7.740413e-03, 0, 1.889281e-03), 
               tolerance = 1e-6)
  expect_equal(ligex@H[[2]][6, 1:5], c(0, 0.004608672, 0, 0, 0.001792484), 
               tolerance = 1e-6)
  expect_equal(ligex@V[[1]][5, 6:10], c(0.1577549, 0, 0.5709908, 0, 0.9109557), 
               tolerance = 1e-6)
  expect_equal(ligex@V[[2]][2, 1:5], c(0, 0, 0, 1.7642890, 1.7410389), 
               tolerance = 1e-6)
  expect_equal(unname(ligex@W[2, 1:5]), c(0, 22.6375047, 0, 1.034460, 2.619777), 
               tolerance = 1e-6)
})

# TODO: Add tests for different lambda setting
# TODO: Add tests for updating k, or with new data
# Note that this should be done so that total test time still < 60s or skip on cran

# Tests for shared factor neighborhood quantile alignment
####################################################################################
context("Quantile alignment")

ligex <- quantileAlignSNF(ligex, knn_k = 20, k2 = 200, resolution = 1)

test_that("Dimensions and lengths are correct", {
  expect_equal(dim(ligex@H.norm), c(494, 15))
  expect_equal(length(ligex@alignment.clusters), 494)
  expect_equal(length(ligex@clusters), 494)
  expect_equal(levels(ligex@clusters), c("0", "1", "2", "3", "4"))
})

test_that("Alignment and clustering are correct", {
  expect_equal(ligex@H.norm[5, 1:5], c(4.141647e-03, 0, 0.0020737469, 0, 0),
               tolerance = 1e-6)
  expect_equal(ligex@H.norm[405, 1:5], c(0.0016112599, 0.0037065406, 0.007530993, 1.360769e-04, 3.096619e-05),
               tolerance = 1e-6)
  expect_equal(as.character(ligex@alignment.clusters[3]), "1")
  expect_equal(as.character(ligex@alignment.clusters[203]), "0")
})

# TODO: Add tests for saving of SNF (once better parameter setting in place)
# TODO: Add tests for different knn_k and k2 settings 

ligex <- quantileAlignSNF(ligex, knn_k = 20, k2 = 200, resolution = 1.5)

test_that("Number of clusters increases to correct value", {
  expect_equal(levels(ligex@clusters), c("0", "1", "2", "3", "4", "5", "6", "7"))
})

test_that("New alignment and clustering are correct", {
  expect_equal(ligex@H.norm[5, 1:5], c(0.004141647, 0, 0.002073747, 0, 0),
               tolerance = 1e-6)
  expect_equal(ligex@H.norm[405, 1:5], c(1.875362e-03, 1.949115e-02, 7.738711e-03, 6.670009e-05, 6.529448e-04),
               tolerance = 1e-6)
  expect_equal(as.character(ligex@alignment.clusters[3]), "4")
  expect_equal(as.character(ligex@alignment.clusters[203]), "0")
})

# Tests for dimensional reduction
# These are included here because these functions are object dependent,
# to avoid recalculating factorization for fresh object as it's time-consuming
# TODO: Add smaller example object (with H, H.norm included) that could be loaded
# to make tests more modular 
####################################################################################
context("Dimensional reduction")

ligex <- runTSNE(ligex, use.raw = F, rand.seed = 1, method = 'Rtsne')
ligex_rawtsne <- runTSNE(ligex, use.raw = T, rand.seed = 1, method = 'Rtsne')

test_that("Dimensions are correct", {
  expect_equal(dim(ligex@tsne.coords), c(494, 2))
  expect_equal(dim(ligex_rawtsne@tsne.coords), c(494, 2))
})


# Tests for plotting functions
# Again, included here because these functions are object dependent (see above)
####################################################################################
context("Plotting")

geneloadings_plots <- plotGeneLoadings(ligex, return.plots = T)

test_that("plotGeneLoadings returns correct number of assembled ggplot objects", {
  expect_equal(length(geneloadings_plots), ncol(ligex@H[[1]]))
  expect_is(geneloadings_plots[[1]], class = c("ggassemble", "gg", "ggplot"))
})

plotgenes_plots <- plotGene(ligex, gene = 'CD8A', use.raw = T, return.plots = T)
test_that("plotGene raw returns correct ggplot objects", {
  expect_equal(length(plotgenes_plots), length(ligex@raw.data))
  expect_is(plotgenes_plots[[1]], class = c("ggplot"))
  expect_equal(unname(plotgenes_plots[[1]]$data$gene[45:50]), 
               c(NA, NA, 2, 2, NA, NA))
  expect_equal(unname(plotgenes_plots[[2]]$data$gene[45:50]), 
               c(NA, NA, NA, 3.281916639, NA, NA))
})

plotgenes_plots <- plotGene(ligex, gene = 'CD8A', return.plots = T)
test_that("plotGene normalizes correctly", {
  expect_equal(length(plotgenes_plots), length(ligex@raw.data))
  expect_equal(unname(plotgenes_plots[[1]]$data$gene[45:50]), 
               c(NA, NA, 3.083491754, 2.879568640, NA, NA))
  expect_equal(unname(plotgenes_plots[[2]]$data$gene[45:50]), 
               c(NA, NA, NA, 2.098256709, NA, NA))
})

plotgenes_plots <- plotGene(ligex, gene = 'CD8A', use.scaled = T, return.plots = T)
test_that("plotGene scales correctly", {
  expect_equal(length(plotgenes_plots), length(ligex@raw.data))
  expect_equal(unname(plotgenes_plots[[1]]$data$gene[45:50]), 
               c(NA, NA, 14.69475989, 14.46124833, NA, NA))
  expect_equal(unname(plotgenes_plots[[2]]$data$gene[45:50]), 
               c(NA, NA, NA, 13.61593681, NA, NA))
})

plotgenes_plots <- plotGene(ligex, gene = 'CD8A', plot.by = 'none', 
                            return.plots = T)
test_that("plotGene returns correct number of plots", {
  expect_is(plotgenes_plots, class = c("ggplot"))
})

plotfeatures_plots <- plotFeature(ligex, feature = 'nUMI', by.dataset = T, return.plots = T)
test_that("plotFeature returns correct ggplot objects", {
  expect_equal(length(plotfeatures_plots), length(ligex@raw.data))
  expect_is(plotfeatures_plots[[1]], class = c("ggplot"))
  expect_equal(as.character(plotfeatures_plots[[1]]$data['AATGCGTGGCTATG', 'dataset']), 
               'tenx')
  expect_equal(as.character(plotfeatures_plots[[2]]$data['Bcell_233', 'dataset']), 
               'seqwell')
})

# Ensure that riverplot function runs 
# TODO: Perhaps check graphics output with vdiffr to make this a real unit test
seqwell_c <- factor(sapply(colnames(ligex@raw.data[[2]]), function(x) {
  strsplit(x, "_")[[1]][1]
}))
tenx_c <- factor(rep("tenx", ncol(ligex@raw.data[[1]])))
names(tenx_c) <- colnames(ligex@raw.data[[1]])

makeRiverplot(ligex, tenx_c, seqwell_c, min.frac = 0.05)

# Tests for subsetting, object conversion
# Again, included here because these functions are object dependent (see above)
####################################################################################
context("Object subsetting and conversion")

# Subsetting functionality
# First add a new cell.data column
ligex@cell.data[["clusters"]] <- ligex@alignment.clusters
ligex_subset <- subsetLiger(ligex, clusters.use = c(1, 2, 3, 4))

test_that("Returns correct subsetted object", {
  expect_equal(names(ligex_subset@raw.data), c('tenx', 'seqwell'))
  expect_equal(dim(ligex_subset@raw.data[[1]]), c(10982, 119))
  expect_equal(colnames(ligex_subset@raw.data[[1]])[1:3], c("ATGCCAGACAGTCA", 
                                                            "CCAAAGTGTGAGAA", "GACCAAACGACTAC"))
  expect_equal(dim(ligex_subset@raw.data[[2]]), c(6705, 138))
  expect_equal(colnames(ligex_subset@raw.data[[2]])[1:3], c("Myeloid_457", "Myeloid_671",
                                                            "Myeloid_1040"))
  expect_equal(levels(ligex_subset@clusters), c("1", "2", "3", "4"))
  expect_equal(nrow(ligex_subset@cell.data), 257)
  expect_equal(rownames(ligex_subset@cell.data), rownames(ligex_subset@tsne.coords))
})

# TODO: Add tests for ligerToSeurat and seuratToLiger functions 
# after including functionality related to new cell.data slot 

# Reorganization
# Make dummy grouping
ligex@cell.data[["broad_type"]] <- ligex@alignment.clusters
levels(ligex@cell.data[["broad_type"]]) <- c("non-myeloid", "non-myeloid", "myeloid", "non-myeloid",
                                             "myeloid", "non-myeloid", "non-myeloid", "non-myeloid")
ligex_reorg <- reorganizeLiger(ligex, by.feature = "broad_type", new.label = "protocol",
                               remove.missing = F)
gene_union <- length(union(rownames(ligex@raw.data[[1]]), rownames(ligex@raw.data[[2]])))

test_that("Returns correctly organized object", {
  expect_equal(names(ligex_reorg@raw.data), c('non-myeloid', 'myeloid'))
  expect_equal(dim(ligex_reorg@raw.data[[1]]), c(gene_union, 370))
  expect_equal(rownames(ligex@cell.data[ligex@cell.data$broad_type == "myeloid", ]),
               colnames(ligex_reorg@raw.data[["myeloid"]]))
  expect_equal(dim(ligex_reorg@raw.data[[2]]), c(gene_union, 124))
  expect_equal(rownames(ligex@cell.data[ligex@cell.data$broad_type == "non-myeloid", ]),
               colnames(ligex_reorg@raw.data[["non-myeloid"]]))
})

test_that("Returns correct cell.data columns", {
  expect_true("protocol" %in% colnames(ligex_reorg@cell.data))
})
       

# Tests for imputing query datasets
# Since this function depends on the H.norm matrix, optimizeALS and quantileAlignSNF
# should be performed before this test
####################################################################################
context("Imputing query datasets")

ligex <- imputeKNN(ligex, reference = 'seqwell', weight = TRUE)

test_that("List names and dimensions correct", {
  expect_equal(names(ligex@raw.data), c('tenx', 'seqwell'))
  expect_equal(dim(ligex@raw.data[['tenx']]), c(6712, 250))
})

test_that("Imputation results correct", {
  expect_equivalent(ligex@raw.data[['tenx']][1, 1:5], c(0.1010272982, 0.1480327620, 0.2882094083, 0.1549154381, 0.1010290263),
               tolerance = 1e-8)
})

