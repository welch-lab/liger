## Tests for iNMF and quantile alignment

set.seed(1)

# pbmc.file <- system.file('tests', 'testdata', 'small_pbmc_data.RDS', package = 'liger')
pbmc.file <- "../testdata/small_pbmc_data.RDS"
pbmc.small <- readRDS(pbmc.file)

# preprocessing steps (object required for tests)
ligex <- createLiger(raw.data = pbmc.small, take.gene.union = F,
                     remove.missing = T)
ligex <- normalize(ligex)
ligex <- selectGenes(ligex, var.thresh = c(0.3, 0.9), do.plot = F)
ligex <- scaleNotCenter(ligex)

# Tests for iNMF factorization
####################################################################################
context("iNMF factorization")

ligex <- optimizeALS(ligex, k = 15, lambda = 5, rand.seed = 1, max.iters = 100, thresh = 1e-4)

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

# Tests for shared factor neighborhood quantile alignment -- deprecated
####################################################################################
#context("Deprecated -- Quantile alignment")

#ligex <- quantileAlignSNF(ligex, knn_k = 20, k2 = 200, resolution = 1)

#test_that("Dimensions and lengths are correct", {
#  expect_equal(dim(ligex@H.norm), c(494, 15))
  # expect_equal(length(ligex@alignment.clusters), 494)
#  expect_equal(length(ligex@clusters), 494)
#  expect_equal(levels(ligex@clusters), c("1", "10", "11", "12", "13", "14", "15", "2", "3", "4",
#                                         "5", "6", "7", "8", "9"))
#})

#test_that("Alignment and clustering are correct", {
#  expect_equal(ligex@H.norm[5, 1:5], c(0.004141647, 0, 0.002073747, 0, 0),
#               tolerance = 1e-6)
#  expect_equal(ligex@H.norm[405, 1:5], c(0.0022715258, 0.0194911522, 0.0077549767, 0, 0.0003304383),
#               tolerance = 1e-6)
#  expect_equal(as.character(ligex@clusters[3]), "13")
  # expect_equal(as.character(ligex@alignment.clusters[203]), "0")
#})

# TODO: Add tests for saving of SNF (once better parameter setting in place)
# TODO: Add tests for different knn_k and k2 settings

#### Tests for new functions
# ligex <- quantileAlignSNF(ligex, knn_k = 20, k2 = 200, resolution = 1.5)
#
# test_that("Number of clusters increases to correct value", {
#   expect_equal(levels(ligex@clusters), c("0", "1", "2", "3", "4", "5", "6", "7"))
# })
#
# test_that("New alignment and clustering are correct", {
#   expect_equal(ligex@H.norm[5, 1:5], c(0.004141647, 0, 0.002073747, 0, 0),
#                tolerance = 1e-6)
#   expect_equal(ligex@H.norm[405, 1:5], c(1.875362e-03, 1.949115e-02, 7.738711e-03, 6.670009e-05, 6.529448e-04),
#                tolerance = 1e-6)
#   expect_equal(as.character(ligex@alignment.clusters[3]), "4")
#   expect_equal(as.character(ligex@alignment.clusters[203]), "0")
# })

# Tests for shared factor neighborhood quantile alignment
####################################################################################
context("Quantile alignment")

ligex <- quantile_norm(ligex, eps = 3, knn_k = 100)

test_that("New alignment and clustering are correct", {
  expect_equal(levels(ligex@clusters), c("1" ,"10","11","12","13","14","15","2","3","4","5","6","7","8","9" ))
  expect_equal(sum(as.character(ligex@clusters) == "5"), 27)
  expect_equal(ligex@H.norm[5, 1:5], c(0.004141647, 0, 0.002073747, 0, 0),
               tolerance = 1e-6)
  expect_equal(ligex@H.norm[405, 1:5], c(2.756442e-03, 6.151422e-02, 9.366456e-03, 4.043478e-07, 1.632740e-03),
               tolerance = 1e-6)
})

ligex <- quantile_norm(ligex)

test_that("Dimensions and lengths are correct", {
  expect_equal(dim(ligex@H.norm), c(494, 15))
  expect_equal(length(ligex@clusters), 494)
  expect_equal(sum(as.character(ligex@clusters) == "5"), 54)
  expect_equal(levels(ligex@clusters), c("1" ,"10","11","12","13","14","15","2","3","4","5","6","7","8","9" ))
})

test_that("Alignment and clustering are correct", {
  expect_equal(ligex@H.norm[5, 1:5], c(0.004141647, 0, 0.002073747, 0, 0),
               tolerance = 1e-6)
  expect_equal(ligex@H.norm[405, 1:5], c(0.0022715258, 0.0194911522, 0.0077549767, 0, 0.0003304383),
               tolerance = 1e-6)
})

# Tests for Louvain Community Detection
####################################################################################
context("Louvain Community Detection")

ligex <- louvainCluster(ligex)

test_that("Dimensions and lengths are correct", {
  expect_equal(dim(ligex@H.norm), c(494, 15))
  expect_equal(length(ligex@clusters), 494)
  expect_equal(levels(ligex@clusters), c("0", "1", "2", "3", "4", "5", "6", "7", "8"))
  expect_equal(sum(as.character(ligex@clusters) == "5"), 50)
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
               c(NA, NA, 14.69475989, 14.46124833, NA, NA), tolerance = 1e-5)
  expect_equal(unname(plotgenes_plots[[2]]$data$gene[45:50]),
               c(NA, NA, NA, 13.61593681, NA, NA), tolerance = 1e-5)
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
# seqwell_c <- factor(sapply(colnames(ligex@raw.data[[2]]), function(x) {
#   strsplit(x, "_")[[1]][1]
# }))
# tenx_c <- factor(rep("tenx", ncol(ligex@raw.data[[1]])))
# names(tenx_c) <- colnames(ligex@raw.data[[1]])

# makeRiverplot(ligex, tenx_c, seqwell_c, min.frac = 0.05)

# Tests for subsetting, object conversion
# Again, included here because these functions are object dependent (see above)
####################################################################################
context("Object subsetting and conversion")

# Subsetting functionality
# First add a new cell.data column
ligex@cell.data[["clusters"]] <- ligex@clusters
ligex_subset <- subsetLiger(ligex, clusters.use = c(1, 2, 3, 4))

test_that("Returns correct subsetted object", {
  expect_equal(names(ligex_subset@raw.data), c('tenx', 'seqwell'))
  expect_equal(dim(ligex_subset@raw.data[[1]]), c(10090, 116))
  expect_equal(colnames(ligex_subset@raw.data[[1]])[1:3], c("ATGCCAGACAGTCA",
                                                            "CCAAAGTGTGAGAA", "TGATACCTCACTAG"))
  expect_equal(dim(ligex_subset@raw.data[[2]]), c(6689, 115))
  expect_equal(colnames(ligex_subset@raw.data[[2]])[1:3], c("Myeloid_457", "Myeloid_1683",
                                                            "Myeloid_345"))
  expect_equal(levels(ligex_subset@clusters), c("1", "2", "3", "4"))
  expect_equal(nrow(ligex_subset@cell.data), 231)
  expect_equal(rownames(ligex_subset@cell.data), rownames(ligex_subset@tsne.coords))
})

# create "pseudorandom" set of cells to downsample
cells.use = c('CACTGAGACAGTCA', 'CD8_124', 'Bcell_103', 'Bcell_17', 'Bcell_236', 'GGGCCAACCTTGGA',
              'ACCTCCGATATGCG', 'Bcell_242', 'CD4_407', 'CD8_265', 'GACAGTACGAGCTT', 'GACCCTACTAAAGG',
              'DC_37', 'CD4_35', 'ATACTCTGGTATGC', 'CAAAGCTGAAAGTG', 'AGCACTGATGCTTT', 'Bcell_280',
              'CD4_503', 'DC_97', 'NK_192', 'GGCACGTGGCTTAG', 'CGTTTAACTGGTCA', 'TATGCGGATAACCG',
              'TGTGATCTGACACT', 'CD4_500', 'GGCGGACTTGAACC', 'ATGTAAACACCTCC', 'CD4_539', 'DC_12')
ligex_subset <- subsetLiger(ligex, cells.use = cells.use)

test_that("Returns correct subsetted object", {
  expect_equal(names(ligex_subset@raw.data), c('tenx', 'seqwell'))
  expect_equal(dim(ligex_subset@raw.data[[1]]), c(5233, 14))
  expect_equal(colnames(ligex_subset@raw.data[[1]])[1:3], c("CACTGAGACAGTCA",
                                                            "GGGCCAACCTTGGA", "ACCTCCGATATGCG"))
  expect_equal(dim(ligex_subset@raw.data[[2]]), c(4670, 16))
  expect_equal(colnames(ligex_subset@raw.data[[2]])[1:3], c("CD8_124", "Bcell_103",
                                                            "Bcell_17"))
  expect_equal(levels(ligex_subset@clusters), c("0", "1", "2", "3", "4", "5", "6", "7"))
  expect_equal(nrow(ligex_subset@cell.data), 30)
  expect_equal(rownames(ligex_subset@cell.data), rownames(ligex_subset@tsne.coords))
})

# TODO: Add tests for ligerToSeurat and seuratToLiger functions
# after including functionality related to new cell.data slot

# Reorganization
# Make dummy grouping
ligex@cell.data[["broad_type"]] <- ligex@clusters
levels(ligex@cell.data[["broad_type"]]) <- c("non-myeloid", "non-myeloid", "myeloid", "non-myeloid",
                                             "myeloid", "non-myeloid", "non-myeloid", "non-myeloid", "myeloid")
ligex_reorg <- reorganizeLiger(ligex, by.feature = "broad_type", new.label = "protocol",
                               remove.missing = F)
gene_union <- length(union(rownames(ligex@raw.data[[1]]), rownames(ligex@raw.data[[2]])))

test_that("Returns correctly organized object", {
  expect_equal(names(ligex_reorg@raw.data), c('non-myeloid', 'myeloid'))
  expect_equal(dim(ligex_reorg@raw.data[[1]]), c(gene_union, 352))
  expect_equal(rownames(ligex@cell.data[ligex@cell.data$broad_type == "myeloid", ]),
               colnames(ligex_reorg@raw.data[["myeloid"]]))
  expect_equal(dim(ligex_reorg@raw.data[[2]]), c(gene_union, 142))
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

ligex.ds <- imputeKNN(ligex, reference = "seqwell", weight = FALSE, knn_k = 50)

test_that("List names and dimensions correct", {
  expect_equal(names(ligex.ds@raw.data), c("tenx", "seqwell"))
  expect_equal(dim(ligex.ds@raw.data[["tenx"]]), c(6712, 250))
})

test_that("Imputation results correct", {
  expect_equivalent(ligex.ds@raw.data[["tenx"]][1, 1:5], c(0.124664288, 0.167032324, 0.236469695,
                                                           0.124664288, 0.124664288),
    tolerance = 1e-8
  )
  expect_equivalent(ligex.ds@raw.data[["tenx"]][400, 1:5], c(0.222377206, 0.155470876, 0.239740853,
                                                             0.222377206, 0.222377206),
                    tolerance = 1e-8
  )
})

ligex.ds <- imputeKNN(ligex, reference = "seqwell", knn_k = 50)

test_that("List names and dimensions correct", {
  expect_equal(names(ligex.ds@raw.data), c("tenx", "seqwell"))
  expect_equal(dim(ligex.ds@raw.data[["tenx"]]), c(6712, 250))
})

test_that("Imputation results correct", {
  expect_equivalent(ligex.ds@raw.data[["tenx"]][1, 1:5], c(0.124636160, 0.167004123, 0.236442235,
                                                           0.124643728, 0.124640972),
    tolerance = 1e-8
  )
  expect_equivalent(ligex.ds@raw.data[["tenx"]][400, 1:5], c(0.222271152, 0.155450162, 0.239728702,
                                                             0.222305523, 0.222283342),
                    tolerance = 1e-8
  )
})

# Tests for running Wilcoxon tests
# Since this function depends on the cluster assignments, optimizeALS and quantileAlignSNF
# should be performed before this test
####################################################################################
context("Running Wilcoxon tests")

wilcox.results <- runWilcoxon(ligex, data.use = 2, compare.method = "clusters")

test_that("Wilcoxon test for 'clusters' results correct", {
  expect_equal(wilcox.results[1, 1], "AAED1")
  expect_equivalent(wilcox.results[1, 7], 0.526369034, tolerance = 1e-8)
})

wilcox.results <- runWilcoxon(ligex, compare.method = "datasets")

test_that("Wilcoxon test for 'datasets' results correct", {
  expect_equal(wilcox.results[1, 1], "ISG15")
  expect_equivalent(wilcox.results[1, 7], 0.390391202, tolerance = 1e-8)
})

# Tests for Creating gene-peak regulation network
# Since this function depends on the cluster assignments, optimizeALS and quantileAlignSNF
# should be performed before this test
####################################################################################
# context("Linking Genes and Peaks")

# temp_bed <- data.frame(
#   V1 = c("chr1", "chr1"),
#   V2 = c(1353799, 1337275),
#   V3 = c(1356824, 1342693),
#   V4 = c("ANKRD65", "MRPL20")
# )

# write.table(temp_bed,
#   file = "../testdata/temp_coords.bed", append = TRUE,
#   quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".",
#   row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
#   fileEncoding = ""
# )

# set.seed(17)
# psudo_data <- rnorm(200, mean = 0.5, sd = 0.1)
# gmat.small <- Matrix(
#   data = psudo_data, nrow = 2, ncol = 100,
#   dimnames = list(c("MRPL20", "ANKRD65"), paste0("cell_", seq(1:100))), sparse = T
# )
# pmat.small <- Matrix(
#   data = c(psudo_data[1:100], psudo_data[101:200] + 0.2), nrow = 2, ncol = 100,
#   dimnames = list(c("chr1:1821507-1822007", "chr1:1850611-1851111"), paste0("cell_", seq(1:100))), sparse = T
# )

# regnet <- linkGenesAndPeaks(gene_counts = gmat.small, peak_counts = pmat.small, dist = "spearman",
#                             alpha = 0.05, path_to_coords = "../testdata/temp_coords.bed") # about 40 mins

# test_that("Testing linkage between gene and peaks", {
#   expect_equivalent(regnet[1, 1], 0.6340474, tolerance = 1e-7)
#   expect_equivalent(regnet[2, 2], 0.6929733, tolerance = 1e-7)
# })
# unlink("../testdata/temp_coords.bed")

# Tests for runGSEA
# Since this function depends on the cluster assignments, optimizeALS and quantil_norm
# should be performed before this test
####################################################################################
# context("GSEA testing on metagenes")
# gsea <- runGSEA(ligex)

# test_that("Tests top pathways and NES values", {
#   expect_equal(gsea[[1]][1, 1][[1]], "Axon guidance")
#   expect_equivalent(gsea[[1]][1, 2][[1]], 9.999e-05, tolerance = 1e-7)
# })
