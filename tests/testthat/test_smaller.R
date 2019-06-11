## Tests for object creation and preprocessing

set.seed(1)

pbmc.file <- system.file('tests', 'testdata', 'smaller_pbmc_data.RDS', package = 'liger')
load(pbmc.file)
pbmc.small <- smaller_pbmc

# Tests for object creation 
####################################################################################
context("Object creation")

ligex <- createLiger(raw.data = pbmc.small, make.sparse = T, take.gene.union = F,
                     remove.missing = T)
test_that("Object instantiation creates liger object", {
  expect_is(ligex, "liger")
})

test_that("Dataset names passed correctly", {
  expect_identical(names(ligex@raw.data), c('tenx', 'seqwell'))
})

test_that("Sparse matrices created", {
  expect_is(ligex@raw.data[[1]], "CsparseMatrix")
})

test_that("cell.data created correctly", {
  expect_is(ligex@cell.data, "data.frame")
  expect_equal(rownames(ligex@cell.data)[1:10], colnames(ligex@raw.data[[1]])[1:10])
  expect_equal(unname(ligex@cell.data[["nUMI"]][3]), 307)
  expect_equal(unname(ligex@cell.data[["nGene"]][253]), 461)
  expect_equal(as.character(ligex@cell.data[["dataset"]][3]), "tenx")
  expect_equal(as.character(ligex@cell.data[["dataset"]][253]), "seqwell")
})

ligex.nofil <- createLiger(raw.data = pbmc.small, make.sparse = T, take.gene.union = F,
                           remove.missing = F)
ligex.union <- createLiger(raw.data = pbmc.small, make.sparse = T, take.gene.union = T,
                           remove.missing = T)
ligex.union.nofil <- createLiger(raw.data = pbmc.small, make.sparse = T, take.gene.union = T,
                                 remove.missing = F)

test_that("Dimensions correct for filtered and not filtered", {
  expect_equal(dim(ligex@raw.data[[1]]), c(1517, 250))
  expect_equal(dim(ligex@raw.data[[2]]), c(2014, 244))
  expect_equal(dim(ligex.nofil@raw.data[[1]]), c(4000, 250))
  expect_equal(dim(ligex.nofil@raw.data[[2]]), c(2014, 244))
  expect_equal(dim(ligex.union@raw.data[[1]]), c(3294, 250))
  expect_equal(dim(ligex.union@raw.data[[2]]), c(3294, 244))
  expect_equal(dim(ligex.union.nofil@raw.data[[1]]), c(5767, 250))
  expect_equal(dim(ligex.union.nofil@raw.data[[2]]), c(5767, 244))
})

rm(ligex.nofil, ligex.union, ligex.union.nofil)

# Tests for data merging
##########################################################################################
context('Sparse merging')

# create fake datasets
dataset1 <- as(matrix(0, nrow = 6, ncol = 5), 'CsparseMatrix')
dataset1[c(1, 5, 14, 18, 21, 28)] <- 1:6
rownames(dataset1) <- paste0('gene', 11:16)
colnames(dataset1) <- paste0('cell', 1:5)

dataset2 <- as(matrix(0, nrow = 6, ncol = 6), 'CsparseMatrix')
dataset2[c(3, 8, 12, 14, 20, 21, 35)] <- 1:7
rownames(dataset2) <- c(paste0('gene', 11:13), paste0('gene', 7:9))
colnames(dataset2) <- paste0('cell', 6:11)

merged <- MergeSparseDataAll(list(dataset1, dataset2))

test_that("Merged entries correct", {
  expect_equal(unname(merged[, 'cell2']), rep(0, 9))
  expect_equal(unname(merged[, 'cell7']), c(0, 2, 0, 0, 0, 0, 0, 0, 3))
  expect_equal(unname(merged['gene12', ]), c(0, 0, 3, 0, 0, 0, 2, 4, 5, 0, 0))
  expect_equal(unname(merged['gene7', ]), rep(0, 11))
  expect_equal(merged['gene13', 'cell9'], 6)
  expect_equal(merged['gene14', 'cell5'], 6)
})

# Tests for normalization
##########################################################################################
context('Normalization')

ligex <- normalize(ligex)

test_that("Dataset names passed correctly", {
  expect_identical(names(ligex@norm.data), c('tenx', 'seqwell'))
})

test_that("Normalization is correct", {
  expect_equal(ligex@norm.data[[1]][5, 2], 0.018518518519, tolerance = 1e-6)
  expect_equal(ligex@norm.data[[1]]['MRPS33', "GCCCATACAGCGTT"], 0.004405286344, tolerance = 1e-6)
  expect_equal(ligex@norm.data[[2]][9, 1], 0.002155172414, tolerance = 1e-6)
  expect_equal(ligex@norm.data[[2]]['SERPINB9', "Bcell_124"], 0.0033921302578, tolerance = 1e-6)
  expect_equal(sum(ligex@norm.data[[1]][, 1]), 1)
  expect_equal(sum(ligex@norm.data[[2]][, 1]), 1)
})

# Tests for gene selection
##########################################################################################
context('Gene selection')

ligex <- selectGenes(ligex, var.thresh = c(0.3, 0.9), do.plot = F)
ligex_higher <- selectGenes(ligex, var.thresh = c(0.5, 0.9), do.plot = F)

# Check for inclusion of significant genes
test_genes <- c('CTBP2', 'CAPN1', 'FCGR1A', 'DEDD', 'RPS18')
gene_inclusion <- sapply(test_genes, function(x) {
  x %in% ligex@var.genes
})
gene_inclusion_higher <- sapply(test_genes, function(x) {
  x %in% ligex_higher@var.genes
})

test_that("Significant genes are selected", {
  expect_true(all(gene_inclusion))
  expect_true(all(gene_inclusion_higher))
})

test_that("Number of genes is correct", {
  expect_equal(length(ligex@var.genes), 69)
  expect_equal(length(ligex_higher@var.genes), 61)
})

ligex_intersect <- selectGenes(ligex, var.thresh = c(0.3, 0.9), do.plot = F, combine = "intersection")

#next two testthats give abnormal output with the current version of the code
test_that("Number of genes is correct for intersection", {
  expect_equal(length(ligex_intersect@var.genes), 10)
})

test_that("Gives warning when no genes selected", {
  expect_error(selectGenes(ligex, var.thresh = c(2.3, 2.3), do.plot = F, 
                           combine = "intersection"))
})

# Tests for gene scaling
##########################################################################################
context('Gene scaling (no centering)')

ligex <- scaleNotCenter(ligex)

test_that("Dataset names passed correctly", {
  expect_identical(names(ligex@scale.data), c('tenx', 'seqwell'))
})

# Genes should now be columns
test_that("Dimensions are correct", {
  expect_equal(dim(ligex@scale.data[[1]]), c(250, 69))
  expect_equal(dim(ligex@scale.data[[2]]), c(244, 69))
})

test_that("Scaling is correct", {
  expect_equal(ligex@scale.data[[1]][3, "S100A6"], 1.2624833128, tolerance = 1e-6)
  expect_equal(ligex@scale.data[[1]][50, "C14orf166"], 1.4257104609, tolerance = 1e-6)
  expect_equal(ligex@scale.data[[2]][5, "RPLP1"], 0.06835049067, tolerance = 1e-6)
  expect_equal(ligex@scale.data[[2]][61,"CD8A"], 3.5040385517, tolerance = 1e-6)
})



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
  expect_equal(dim(ligex@V[[1]]), c(15, 69))
  expect_equal(dim(ligex@V[[2]]), c(15, 69))
  expect_equal(dim(ligex@W), c(15, 69))
})

test_that("Factorization is correct", {
  expect_equal(ligex@H[[1]][6, 1:5], c(0.019174771654, 0, 0,0.003274994786,0.002616858358), 
               tolerance = 1e-6)
  expect_equal(ligex@H[[2]][6, 1:5], c(0, 0.001869488079, 0, 0, 0), 
               tolerance = 1e-6)
  expect_equal(ligex@V[[1]][5, 6:10], c(0, 0, 0, 0.2139965426, 0), 
               tolerance = 1e-6)
  expect_equal(ligex@V[[2]][2, 1:5], c(0, 0, 0.8185788919, 0, 0), 
               tolerance = 1e-6)
  expect_equal(unname(ligex@W[2, 1:5]), c(0, 0.05821537459, 2.46681953578, 0, 0.42598114733), 
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
  expect_equal(levels(ligex@clusters), c("0", "1", "2"))
})

test_that("Alignment and clustering are correct", {
  expect_equal(ligex@H.norm[5, 1:5], c(0.0080483240699, 0, 0.0002286524507, 0.0008272908936, 0),
               tolerance = 1e-6)
  expect_equal(ligex@H.norm[405, 1:5], c(1.106839925e-02, 1.247105730e-02, 0, 0, 2.596066313e-05),
               tolerance = 1e-6)
  expect_equal(as.character(ligex@alignment.clusters[3]), "2")
  expect_equal(as.character(ligex@alignment.clusters[203]), "1")
})

# TODO: Add tests for saving of SNF (once better parameter setting in place)
# TODO: Add tests for different knn_k and k2 settings 

ligex <- quantileAlignSNF(ligex, knn_k = 20, k2 = 200, resolution = 1.5)

test_that("Number of clusters increases to correct value", {
  expect_equal(levels(ligex@clusters), c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
})

test_that("New alignment and clustering are correct", {
  expect_equal(ligex@H.norm[5, 1:5], c(0.0080483240699, 0, 0.0002286524507, 0.0008272908936, 0),
               tolerance = 1e-6)
  expect_equal(ligex@H.norm[405, 1:5], c(1.026702233e-02, 1.247105730e-02, 0, 0, 3.347543122e-05),
               tolerance = 1e-6)
  expect_equal(as.character(ligex@alignment.clusters[3]), "2")
  expect_equal(as.character(ligex@alignment.clusters[203]), "1")
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
  expect_equal(dim(ligex@dr.coords[["tsne"]]), c(494, 2))
  expect_equal(dim(ligex_rawtsne@dr.coords[["tsne"]]), c(494, 2))
})


# Tests for plotting functions
# Again, included here because these functions are object dependent (see above)
####################################################################################
context("Plotting")

geneloadings_plots <- plotGeneLoadings(ligex, return.plots = T)

#liger in arc connect RStudio doesn't recognize this function (?)
test_that("plotGeneLoadings returns correct number of assembled ggplot objects", {
  expect_equal(length(geneloadings_plots), ncol(ligex@H[[1]]))
  expect_is(geneloadings_plots[[1]], class = c("ggassemble", "gg", "ggplot"))
})

plotgenes_plots <- plotGene(ligex, gene = 'CD8A', return.plots = T)
test_that("plotGene returns correct ggplot objects", {
  expect_equal(length(plotgenes_plots), length(ligex@raw.data))
  expect_is(plotgenes_plots[[1]], class = c("ggplot"))
  expect_equal(unname(plotgenes_plots[[1]]$data$gene[45:50]), 
               c(NA, NA, 5.802051540, 5.515470544, NA, NA))
})

#liger in arc connect RStudio doesn't recognize this function (?)
plotfeatures_plots <- plotFeature(ligex, feature = 'nUMI', by.dataset = T, return.plots = T)
test_that("plotFeature returns correct ggplot objects", {
  expect_equal(length(plotfeatures_plots), length(ligex@raw.data))
  expect_is(plotfeatures_plots[[1]], class = c("ggplot"))
  expect_equal(rownames(plotfeatures_plots[[1]]$data)[1:5], 
               c("AATGCGTGGCTATG", "GAAAGATGATTTCC", "TTCCAAACTCCCAC", "CACTGAGACAGTCA",
                 "GACGGCACACGGGA"))
})





# Tests for subsetting, object conversion
# Again, included here because these functions are object dependent (see above)
####################################################################################
context("Object subsetting and conversion")

# Subsetting functionality
# First add a new cell.data column
ligex@cell.data[["clusters"]] = ligex@alignment.clusters
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
  expect_equal(rownames(ligex_subset@cell.data), rownames(ligex_subset@dr.coords["tsne"]))
})

# TODO: Add tests for ligerToSeurat and seuratToLiger functions 
# after including functionality related to new cell.data slot 

