## Tests for object creation and preprocessing

# pbmc.file <- system.file('tests', 'testdata', 'small_pbmc_data.RDS', package = 'liger')
#pbmc.file <- "../testdata/small_pbmc_data.RDS"
#pbmc.small <- readRDS(pbmc.file)
data("pbmc", package = "rliger")
rawDataList <- getMatrix(pbmc, "rawData")

withNewH5Copy <- function(fun) {
    ctrlpath.orig <- system.file("extdata/ctrl.h5", package = "rliger")
    stimpath.orig <- system.file("extdata/stim.h5", package = "rliger")
    if (file.exists("ctrltest.h5")) file.remove("ctrltest.h5")
    if (file.exists("stimtest.h5")) file.remove("stimtest.h5")
    file.copy(ctrlpath.orig, "ctrltest.h5")
    file.copy(stimpath.orig, "stimtest.h5")
    fun(list(ctrl = "ctrltest.h5", stim = "stimtest.h5"))
    if (file.exists("ctrltest.h5")) file.remove("ctrltest.h5")
    if (file.exists("stimtest.h5")) file.remove("stimtest.h5")
}

closeH5Liger <- function(object) {
    for (d in names(object)) {
        if (isH5Liger(object, d)) {
            h5file <- getH5File(object, d)
            h5file$close()
        }
    }
}

### IMPORTANT DEVELOPER NOTICE below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# When writing H5 related unit tests, please follow this template:
test_that("<topic> - on-disk", {
    withNewH5Copy(
        function(rawList, arg1, arg2) {
            pbmc <- createLiger(rawList)
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Then whatever test with pbmc. For example:
            expect_true(isH5Liger(pbmc))
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # And must close with:
            closeH5Liger(pbmc)
        }
    )
})
### IMPORTANT DEVELOPER NOTICE above %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Object creation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Object creation")
test_that("liger object creation - in memory", {
    pbmc2 <- createLiger(rawData = rawDataList)
    expect_is(pbmc2, "liger")
    expect_error(createLiger(rawData = "hi"),
                 "`rawData` has to be a named list.")
    expect_error(createLiger(rawData = rawDataList, modal = letters[1:3]),
                 "Wrong length of `modal`. ")
    ldList <- datasets(pbmc)
    cellmeta <- cellMeta(pbmc)
    pbmc2 <- createLiger(rawData = ldList, cellMeta = cellmeta,
                         addPrefix = FALSE)
    expect_identical(cellMeta(pbmc), cellMeta(pbmc2))

    pbmc <- removeMissing(pbmc)
    pbmc <- runGeneralQC(pbmc, pattern = "^S100",
                         features = rownames(ldList[[1]]))
    expect_true(all(c("featureSubset_pattern", "featureSubset_name") %in%
                      colnames(cellMeta(pbmc))))

    pbmc <- runGeneralQC(pbmc, pattern = list(p1 = "^S100", p2 = "^MT"),
                         features = list(f1 = letters,
                                         f2 = rownames(ldList[[2]])[6:10]))
    expect_true(all(c("p1", "p2", "f1", "f2") %in% colnames(cellMeta(pbmc))))
})

test_that("liger object creation - on disk", {
    withNewH5Copy(
        function(rawList) {
            expect_error(createLiger(rawList, formatType = "rliger"),
                         "Specified `formatType` '")

            # Customized paths
            barcodesName <- "matrix/barcodes"
            rawData <- "matrix/data"
            indicesName <- "matrix/indices"
            indptrName <- "matrix/indptr"
            genesName <- "matrix/features/name"
            pbmc <- createLiger(rawList,
                                formatType = NULL,
                                barcodesName = barcodesName,
                                dataName = rawData,
                                indicesName = indicesName,
                                indptrName = indptrName,
                                genesName = genesName)
            expect_is(pbmc, "liger")
            expect_true(isH5Liger(pbmc))
            closeH5Liger(pbmc)

            # Preset paths
            pbmc <- createLiger(rawList, formatType = "10X")
            expect_is(pbmc, "liger")
            expect_true(isH5Liger(pbmc))
            expect_is(rawData(dataset(pbmc, "ctrl")), "H5D")
            expect_is(rawData(dataset(pbmc, "stim")), "H5D")
            expect_is(getH5File(pbmc, "ctrl"), "H5File")
            expect_is(getH5File(pbmc, "stim"), "H5File")
            closeH5Liger(pbmc)
        }
    )
})

test_that("ligerDataset (in memory) object creation", {
    expect_error(createLigerDataset(),
                 "At least one type of")

    ld <- createLigerDataset(rawData = rawDataList[[1]], modal = "atac")
    expect_is(ld, "ligerATACDataset")
    data("pbmc")
    pbmc <- normalize(pbmc)
    normDataList <- getMatrix(pbmc, "normData")
    ld <- createLigerDataset(normData = normDataList[[1]], modal = "rna")
    expect_is(ld, "ligerRNADataset")
    expect_null(rawData(ld))

    pbmc <- selectGenes(pbmc)
    pbmc <- scaleNotCenter(pbmc)
    scaledMat <- scaleData(pbmc, dataset = "ctrl")
    featuremeta <- featureMeta(dataset(pbmc, "ctrl"))
    ld <- createLigerDataset(scaleData = scaledMat, featureMeta = featuremeta)
    expect_equal(length(varFeatures(pbmc)), nrow(ld))
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Tests for data merging
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Merged sparse matrix", {
    # create fake datasets
    dataset1 <- as(matrix(0, nrow = 6, ncol = 5), 'CsparseMatrix')
    dataset1[c(1, 5, 14, 18, 21, 28)] <- 1:6
    rownames(dataset1) <- paste0('gene', 11:16)
    colnames(dataset1) <- paste0('cell', 1:5)

    dataset2 <- as(matrix(0, nrow = 6, ncol = 6), 'CsparseMatrix')
    dataset2[c(3, 8, 12, 14, 20, 21, 35)] <- 1:7
    rownames(dataset2) <- c(paste0('gene', 11:13), paste0('gene', 7:9))
    colnames(dataset2) <- paste0('cell', 6:11)

    merged <- mergeSparseAll(list(dataset1, dataset2))

    expect_equal(unname(merged[, 'cell2']), rep(0, 9))
    expect_equal(unname(merged[, 'cell7']), c(0, 2, 0, 0, 0, 0, 0, 0, 3))
    expect_equal(unname(merged['gene12', ]), c(0, 0, 3, 0, 0, 0, 2, 4, 5, 0, 0))
    expect_equal(unname(merged['gene7', ]), rep(0, 11))
    expect_equal(merged['gene13', 'cell9'], 6)
    expect_equal(merged['gene14', 'cell5'], 6)
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Normalization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("normalization")
test_that("Normalization - in-memory", {
    pbmc2 <- normalize(pbmc, log = TRUE, scaleFactor = 1e4)
    ctrl.norm <- normData(dataset(pbmc2, "ctrl"))
    expect_gt(sum(ctrl.norm[,1]), 345)

    pbmc <- normalize(pbmc, scaleFactor = 1)
    expect_identical(dim(rawData(dataset(pbmc, "ctrl"))),
                     dim(normData(dataset(pbmc, "ctrl"))))
    expect_identical(dim(rawData(dataset(pbmc, "stim"))),
                     dim(normData(dataset(pbmc, "stim"))))
    ld <- dataset(pbmc, "ctrl")
    for (i in seq_len(ncol(ld))) {
        expect_equal(sum(normData(ld)[, i]), 1, tolerance = 1e-6)
    }
    ld <- dataset(pbmc, "stim")
    for (i in seq_len(ncol(ld))) {
      expect_equal(sum(normData(ld)[, i]), 1, tolerance = 1e-6)
    }

    # For atac peak normalization
    fakePeak <- rawData(ld)
    ld <- as.ligerDataset(ld, modal = "atac")
    rawPeak(ld) <- fakePeak
    datasets(pbmc)[["stim"]] <- ld
    pbmc <- normalizePeak(pbmc, useDatasets = "stim")
    expect_identical(dim(rawPeak(dataset(pbmc, "stim"))),
                     dim(normPeak(dataset(pbmc, "stim"))))
    ld <- dataset(pbmc, "stim")
    for (i in seq_len(ncol(ld))) {
      expect_equal(sum(normPeak(ld)[, i]), 1, tolerance = 1e-6)
    }
})

test_that("Normalize - on disk", {
    withNewH5Copy(
        function(rawList) {
            pbmc <- createLiger(rawList, formatType = "10X")
            pbmc1 <- normalize(pbmc, chunk = 100, log = TRUE, scaleFactor = 1e4)
            expect_is(normData(dataset(pbmc1, "ctrl")), "H5D")
            expect_is(normData(dataset(pbmc1, "stim")), "H5D")
            expect_equal(normData(dataset(pbmc1, "ctrl"))$dims,
                         rawData(dataset(pbmc1, "ctrl"))$dims)
            expect_equal(normData(dataset(pbmc1, "stim"))$dims,
                         rawData(dataset(pbmc1, "stim"))$dims)

            closeH5Liger(pbmc)
        }
    )
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Variable gene selection
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Select variable genes")
test_that("selectGenes", {
    pbmc <- normalize(pbmc, useDatasets = 1)
    expect_error(selectGenes(pbmc, var.thresh = 1:3),
                 "Wrong length of `var.thresh`.")
    expect_warning(selectGenes(pbmc, var.thresh = 0.1),
                   "Dataset \"stim\" is not normalized, skipped")
    pbmc <- normalize(pbmc, useDatasets = 2)
    expect_error(selectGenes(pbmc, unshared = TRUE, unshared.thresh = 1:3),
                 "Wrong length of `unshared.thresh`. ")

    pbmc <- selectGenes(pbmc, unshared = TRUE, unshared.thresh = 0)
    expect_identical(dataset(pbmc, "ctrl")@varUnsharedFeatures, character())
    expect_identical(dataset(pbmc, "ctrl")@varUnsharedFeatures, character())

    pbmc <- selectGenes(pbmc, combine = "inters")
    expect_equal(length(varFeatures(pbmc)), 161)

    expect_warning(selectGenes(pbmc, var.thresh = 3),
                   "No genes were selected.")
    expect_warning(selectGenes(pbmc, num.genes = 100),
                   "Returned number of genes for dataset ctrl differs from ")
    expect_warning(selectGenes(pbmc, num.genes = 200),
                   "Cannot optimize the number of selected genes for ")

    pbmc <- selectGenes(pbmc)
    expect_equal(length(varFeatures(pbmc)), 173)

    glist <- plotVarFeatures(pbmc, combinePlot = FALSE)
    expect_is(glist, "list")
    expect_is(glist[[1]], "ggplot")
    g <- plotVarFeatures(pbmc, combinePlot = TRUE)
    expect_is(g, "ggplot")
})

test_that("selectGenes - on disk", {
    withNewH5Copy(
        function(rawList) {
            pbmc <- createLiger(rawList)
            pbmc <- normalize(pbmc)
            pbmc <- selectGenes(pbmc)
            expect_equal(length(varFeatures(pbmc)), 173)
            closeH5Liger(pbmc)
        }
    )
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Scaling
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
context('Gene scaling (no centering)')
test_that("scaleNotCenter - in-memory", {
    expect_error(scaleNotCenter(pbmc),
                 "No variable feature found. ")
    pbmc <- normalize(pbmc)
    pbmc <- selectGenes(pbmc)
    pbmc <- scaleNotCenter(pbmc)
    expect_equal(length(varFeatures(pbmc)),
                 nrow(scaleData(pbmc, 1)))
    expect_equal(length(varFeatures(pbmc)),
                 nrow(scaleData(pbmc, 2)))
    # Add false unshared features in order to cover the code that scales the
    # unshared features
    pbmc@datasets$ctrl@varUnsharedFeatures <- varFeatures(pbmc)[1:5]
    pbmc <- scaleNotCenter(pbmc)
    expect_equal(nrow(scaleUnsharedData(pbmc, "ctrl")), 5)
    expect_null(scaleUnsharedData(pbmc, "stim"))

    expect_equal(scaleData(pbmc, "ctrl")[3, 5], 0.4693316, tolerance = 1e-6)
    expect_equal(scaleData(pbmc, "stim")[7, 9], 2.360295, tolerance = 1e-6)
})

test_that("scaleNotCenter - on-disk", {
    withNewH5Copy(
        function(rawList) {
            pbmc <- createLiger(rawList)
            pbmc <- normalize(pbmc)
            pbmc <- selectGenes(pbmc)
            # Add false unshared features in order to cover the code that scales the
            # unshared features
            pbmc@datasets$ctrl@varUnsharedFeatures <- varFeatures(pbmc)[1:5]
            pbmc <- scaleNotCenter(pbmc)
            expect_is(scaleData(pbmc, "ctrl"), "H5D")
            expect_equal(scaleData(pbmc, "ctrl")$dims, c(173, 300))
            expect_is(scaleData(pbmc, "stim"), "H5D")
            expect_equal(scaleData(pbmc, "stim")$dims, c(173, 300))
        }
    )
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# liger object methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
context("liger object S3/S4 methods")

test_that("liger S4 methods - getters", {
    expect_output(show(pbmc), "An object of class liger with 600 cells")
    expect_equal(dim(pbmc), c(NA, 600))
    expect_null(rownames(pbmc))
    expect_identical(colnames(pbmc), rownames(cellMeta(pbmc)))
    expect_equal(ncol(pbmc[,pbmc$dataset == "ctrl"]), 300)
    expect_is(datasets(pbmc), "list")
    expect_is(dataset(pbmc), "ligerDataset")
    expect_is(dataset(pbmc, "ctrl"), "ligerDataset")
    expect_is(dataset(pbmc, 2), "ligerDataset")
    expect_equal(names(pbmc), c("ctrl", "stim"))
    expect_equal(length(pbmc), 2)
    expect_is(cellMeta(pbmc), "DFrame")
    expect_null(cellMeta(pbmc, NULL))
    expect_is(cellMeta(pbmc, "dataset"), "factor")
    expect_is(pbmc[["nUMI"]], "numeric")
    expect_is(pbmc$mito, "numeric")
    expect_is(varFeatures(pbmc), "character")
    expect_is(c(pbmc, pbmc), "liger")
    expect_is(fortify(pbmc), "data.frame")
})


