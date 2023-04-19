## Tests for object creation and preprocessing

# pbmc.file <- system.file('tests', 'testdata', 'small_pbmc_data.RDS', package = 'liger')
#pbmc.file <- "../testdata/small_pbmc_data.RDS"
#pbmc.small <- readRDS(pbmc.file)
data("pbmc", package = "rliger2")
rawDataList <- getMatrix(pbmc, "rawData")

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

  if (file.exists(ctrlpath)) file.remove(ctrlpath)
  if (file.exists(stimpath)) file.remove(stimpath)
}

closeH5Liger <- function(object) {
  for (d in names(object)) {
    if (isH5Liger(object, d)) {
      h5file <- getH5File(object, d)
      h5file$close()
    }
  }
}

process <- function(object) {
    object <- normalize(object)
    object <- selectGenes(object)
    object <- scaleNotCenter(object)
    object <- online_iNMF(object, k = 20, miniBatch_size = 100)
    object <- quantileNorm(object)
    object <- runUMAP(object)
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
                 "`var.thresh` has to be a vector of length 2")
    expect_warning(selectGenes(pbmc, var.thresh = 0.1),
                   "Dataset \"stim\" is not normalized, skipped")
    pbmc <- normalize(pbmc, useDatasets = 2)
    expect_error(selectGenes(pbmc, unshared = TRUE, unshared.thresh = 1:3),
                 "`unshared.thresh` has to be a vector of length 2")

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



