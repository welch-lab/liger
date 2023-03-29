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
    expect_error(subsetLiger(pbmc, featureIdx = 1:3),
                 "Feature subscription from liger object")
    expect_error(
        expect_warning(subsetLiger(pbmc,
                                   featureIdx = c("fakeGene1", "fakeGene2")),
                       "2 out of 2 given features were not found"),
        "No feature can be retrieved"
    )
    # TODO TwT
})

context("subset ligerDataset object")
test_that("subsetH5LigerDataset", {
    withNewH5Copy(
        function(rawList, arg1, arg2) {
            pbmc <- createLiger(rawList)
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # Then whatever test with pbmc. For example:
            pbmc <- process(pbmc)
            ctrl <- dataset(pbmc, "ctrl")
            ctrlSmall <- subsetLiger(ctrl, featureIdx = 1:10, cellIdx = 1:10)
            expect_false(isH5Liger(ctrlSmall))
            #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # And must close with:
            closeH5Liger(pbmc)
        }
    )
})
