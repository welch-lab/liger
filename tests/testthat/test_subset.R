data("pbmc", package = "rliger")

withNewH5Copy <- function(fun) {
    ctrlpath.orig <- system.file("extdata/ctrl.h5", package = "rliger")
    stimpath.orig <- system.file("extdata/stim.h5", package = "rliger")
    if (!file.exists(ctrlpath.orig))
        stop("Cannot find original h5 file at: ", ctrlpath.orig)
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
    on.exit({
        if (file.exists(ctrlpath)) unlink(ctrlpath)
        if (file.exists(stimpath)) unlink(stimpath)
    })
    fun(list(ctrl = ctrlpath, stim = stimpath))
}

process <- function(object) {
    object <- normalize(object)
    object <- selectGenes(object)
    object <- scaleNotCenter(object)
    object <- runOnlineINMF(object, k = 20, minibatchSize = 100)
    object <- quantileNorm(object)
    object <- runUMAP(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# liger object creation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("subset liger object")
test_that("subsetLiger", {
    expect_message(a <- subsetLiger("a"), "`object` is not a ")
    expect_identical(a, "a")
    skip_if_not_installed("RcppPlanc")
    pbmc <- process(pbmc)
    expect_error(subsetLiger(pbmc, featureIdx = 1:3),
                 "Feature subscription from a")
    expect_error(
        expect_warning(subsetLiger(pbmc,
                                   featureIdx = c("fakeGene1", "fakeGene2")),
                       "2 out of 2 given features were not found"),
        "No feature can be retrieved"
    )
    expect_is(retrieveCellFeature(pbmcPlot, 1, slot = "H"), "data.frame")
})

context("subset ligerDataset object")
test_that("subsetH5LigerDataset", {
    skip_if_not_installed("RcppPlanc")
    withNewH5Copy(
        function(rawList, arg1, arg2) {
            ctrlfile <- rawList$ctrl
            stimfile <- rawList$stim
            expect_warning(pbmcH5 <- createLiger(rawList), 'deprecated')
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Then whatever test with pbmc. For example:
            pbmcH5 <- process(pbmcH5)
            ctrl <- dataset(pbmcH5, "ctrl")
            ctrlSmall <- subsetLiger(ctrl, featureIdx = 1:10, cellIdx = 1:10, newH5 = FALSE)
            expect_false(isH5Liger(ctrlSmall))
            path <- dirname(h5fileInfo(ctrl, "filename"))
            newName <- file.path(path, "ctrltest.h5.small.h5")
            expect_message(
                subsetLigerDataset(ctrl, featureIdx = 1:10, cellIdx = 1:10,
                                   newH5 = TRUE,
                                   filename = newName,
                                   returnObject = FALSE),
                "Cannot set `returnObject = FALSE`"
            )
            expect_true(file.exists(newName))
            unlink(newName)
            expect_message(
                rliger:::subsetH5LigerDatasetToMem(letters),
                "`object` is not a "
            )
            expect_message(
                rliger:::subsetH5LigerDatasetToMem(dataset(pbmc, "ctrl")),
                "`object` is not HDF5 based."
            )
            valueList <- rliger:::subsetH5LigerDatasetToMem(
                ctrl, 1:20, 1:20, returnObject = FALSE
            )
            expect_is(valueList, "list")

            expect_message(
                rliger:::subsetH5LigerDatasetToH5(letters),
                "`object` is not a"
            )
            expect_message(
                rliger:::subsetH5LigerDatasetToH5(dataset(pbmc, "ctrl")),
                "`object` is not HDF5 based."
            )
            expect_no_error(
                subsetH5LigerDataset(ctrl, 1:20, 1:20)
            )

            ctrlSmallH5 <- rliger:::subsetH5LigerDataset(
                ctrl, 1:20, 1:20, filenameSuffix = "small2"
            )
            newPath <- paste0(ctrlfile, ".small2.h5")
            expect_true(file.exists(newPath))
            unlink(newPath)
            expect_no_error(
                rliger:::subsetH5LigerDataset(ctrl, 1:20, 1:20, newH5 = TRUE,
                                               useSlot = "normData",
                                               filenameSuffix = "small3")
            )
            newPath <- paste0(ctrlfile, ".small3.h5")
            expect_true(file.exists(newPath))
            unlink(newPath)
            expect_no_error(
                subsetH5LigerDataset(ctrl, 1:20, 1:20, useSlot = "scaleData")
            )
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # And must close with:
            closeAllH5(pbmcH5)
        }
    )
})

unlink(grep("ctrltest", list.files(), value = TRUE))
unlink(grep("stimtest", list.files(), value = TRUE))
