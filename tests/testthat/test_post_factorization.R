## Tests for iNMF and quantile alignment

set.seed(1)

pbmc.file <- system.file('tests', 'testdata', 'small_pbmc_data.RDS', package = 'liger')
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
# Note that this should be done so that total test time still < 60s

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

test_that("Returns correct number of assembled ggplot objects", {
  expect_equal(length(geneloadings_plots), ncol(ligex@H[[1]]))
  expect_is(geneloadings_plots[[1]], class = c("ggassemble", "gg", "ggplot"))
})