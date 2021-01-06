## Tests for object creation and preprocessing

set.seed(1)

# pbmc.file <- system.file('tests', 'testdata', 'small_pbmc_data.RDS', package = 'liger')
pbmc.file <- "../testdata/small_pbmc_data.RDS"
pbmc.small <- readRDS(pbmc.file)

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

# note that seqwell data is previously normalized, so nUMI is 10000 for all cells
test_that("cell.data created correctly", {
  expect_is(ligex@cell.data, "data.frame")
  expect_equal(rownames(ligex@cell.data)[1:10], colnames(ligex@raw.data[[1]])[1:10])
  expect_equal(unname(ligex@cell.data[["nUMI"]][3]), 2043)
  expect_equal(unname(ligex@cell.data[["nGene"]][253]), 1534)
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
  expect_equal(dim(ligex@raw.data[[1]]), c(12358, 250))
  expect_equal(dim(ligex@raw.data[[2]]), c(6712, 244))
  expect_equal(dim(ligex.nofil@raw.data[[1]]), c(32738, 250))
  expect_equal(dim(ligex.nofil@raw.data[[2]]), c(6713, 244))
  expect_equal(dim(ligex.union@raw.data[[1]]), c(12621, 250))
  expect_equal(dim(ligex.union@raw.data[[2]]), c(12621, 244))
  expect_equal(dim(ligex.union.nofil@raw.data[[1]]), c(32826, 250))
  expect_equal(dim(ligex.union.nofil@raw.data[[2]]), c(32826, 244))
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
  expect_equal(ligex@norm.data[[1]][6, 3], 0.001957905, tolerance = 1e-6)
  expect_equal(ligex@norm.data[[1]]['FCGR3A', 104], 0.002842063, tolerance = 1e-6)
  expect_equal(ligex@norm.data[[2]][8, 1], 0.000532198, tolerance = 1e-6)
  expect_equal(ligex@norm.data[[2]]['FCGR3A', 110], 0.0003521127, tolerance = 1e-6)
  expect_equal(sum(ligex@norm.data[[1]][, 1]), 1)
  expect_equal(sum(ligex@norm.data[[2]][, 1]), 1)
})

# Tests for gene selection
##########################################################################################
context('Gene selection')

ligex <- selectGenes(ligex, var.thresh = c(0.3, 0.9), do.plot = F)
ligex_higher <- selectGenes(ligex, var.thresh = c(0.5, 0.9), do.plot = F)

# Check for inclusion of significant genes
test_genes <- c('FCGR3A', 'GNLY', 'CD8A', 'CD3D', 'MS4A1')
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
  expect_equal(length(ligex@var.genes), 1907)
  expect_equal(length(ligex_higher@var.genes), 1794)
})

ligex_intersect <- selectGenes(ligex, var.thresh = c(0.3, 0.9), do.plot = F, combine = "intersection")

test_that("Number of genes is correct for intersection", {
  expect_equal(length(ligex_intersect@var.genes), 260)
})

test_that("Gives warning when no genes selected", {
  expect_warning(selectGenes(ligex, var.thresh = c(2.3, 2.3), do.plot = F, 
                             combine = "intersection"))
})

# # Keeping unique here would break iNMF later on but allows us to check number of genes
# ligex_higher <- selectGenes(ligex, num.genes = 950, keep.unique = T, datasets.use = 2)
# test_that("Returns same number of genes as requested", {
#   expect_equal(length(ligex_higher@var.genes), 950)
# })

rm(ligex_intersect)

# Tests for gene scaling
##########################################################################################
context('Gene scaling (no centering)')

ligex <- scaleNotCenter(ligex)

test_that("Dataset names passed correctly", {
  expect_identical(names(ligex@scale.data), c('tenx', 'seqwell'))
})

# Genes should now be columns
test_that("Dimensions are correct", {
  expect_equal(dim(ligex@scale.data[[1]]), c(250, 1907))
  expect_equal(dim(ligex@scale.data[[2]]), c(244, 1907))
})

test_that("Scaling is correct", {
  expect_equal(ligex@scale.data[[1]][3, 5], 0.1571564, tolerance = 1e-6)
  expect_equal(ligex@scale.data[[1]][113, 115], 0.6294048, tolerance = 1e-6)
  expect_equal(ligex@scale.data[[2]][3, 5], 1.253122, tolerance = 1e-6)
  expect_equal(ligex@scale.data[[2]][113, 115], 0.261834, tolerance = 1e-6)
})

# Tests for preliminary calculations
##########################################################################################
context('Calculating mitochondrial proportion')

mito_prop <- getProportionMito(ligex)

test_that("Values calculated correctly and names passed", {
  expect_equal(unname(mito_prop[1:10]), rep(0, 10))
  expect_equal(names(mito_prop), rownames(ligex@cell.data))
})
