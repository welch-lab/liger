data("pbmc", package = "rliger2")

withNewH5Copy <- function(fun) {
    ctrlpath.orig <- system.file("extdata/ctrl.h5", package = "rliger2")
    stimpath.orig <- system.file("extdata/stim.h5", package = "rliger2")
    if (!file.exists(ctrlpath.orig))
        stop("Cannot find original h5 file at: ", ctrlpath.orig)
    if (file.exists("ctrltest.h5")) file.remove("ctrltest.h5")
    if (file.exists("stimtest.h5")) file.remove("stimtest.h5")
    pwd <- getwd()
    ctrlpath <- file.path(pwd, "ctrltest.h5")
    stimpath <- file.path(pwd, "stimtest.h5")
    file.copy(ctrlpath.orig, ctrlpath)
    file.copy(stimpath.orig, stimpath)
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# liger object creation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("subset liger object")
test_that("subsetLiger", {
    expect_warning(a <- subsetLiger("a"), "`object` is not a liger obejct")
    expect_identical(a, "a")
    pbmc <- process(pbmc)
    expect_error(subsetLiger(pbmc, featureIdx = 1:3),
                 "Feature subscription from liger object")
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
    withNewH5Copy(
        function(rawList, arg1, arg2) {
            pbmcH5 <- createLiger(rawList)
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Then whatever test with pbmc. For example:
            pbmcH5 <- process(pbmcH5)
            ctrl <- dataset(pbmcH5, "ctrl")
            ctrlSmall <- subsetLiger(ctrl, featureIdx = 1:10, cellIdx = 1:10)
            expect_false(isH5Liger(ctrlSmall))
            expect_warning(
                subsetLigerDataset(ctrl, newH5 = TRUE,
                                   filename = "ctrltest.h5.small.h5",
                                   returnObject = FALSE),
                "Cannot set `returnObject = FALSE`"
            )
            expect_true(hdf5r::is.h5file("ctrltest.h5.small.h5"))
            expect_warning(
                subsetH5LigerDatasetToMem(letters),
                "`object` is not a ligerDataset obejct."
            )
            expect_warning(
                subsetH5LigerDatasetToMem(dataset(pbmc, "ctrl")),
                "`object` is not HDF5 based."
            )
            valueList <- subsetH5LigerDatasetToMem(ctrl, 1:20, 1:20,
                                                   returnObject = FALSE)
            expect_is(valueList, "list")

            expect_warning(
                subsetH5LigerDatasetToH5(letters),
                "`object` is not a ligerDataset obejct."
            )
            expect_warning(
                subsetH5LigerDatasetToH5(dataset(pbmc, "ctrl")),
                "`object` is not HDF5 based."
            )
            expect_no_error(
                subsetH5LigerDatasetToH5(ctrl, 1:20, 1:20)
            )

            ctrlSmallH5 <- subsetH5LigerDatasetToH5(ctrl, 1:20, 1:20,
                                                    filenameSuffix = "small")
            expect_true(hdf5r::is.h5file("ctrltest.h5.small.h5"))
            expect_no_error(
                subsetH5LigerDatasetToH5(ctrl, 1:20, 1:20, useSlot = "normData",
                                         filenameSuffix = "small")
            )
            file.remove("ctrltest.h5.small.h5")
            expect_no_error(
                subsetH5LigerDatasetToH5(ctrl, 1:20, 1:20,
                                         useSlot = "scaleData")
            )
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # And must close with:
            closeH5Liger(pbmc)
        }
    )
})

