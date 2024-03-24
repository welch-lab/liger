data("pbmc", package = "rliger")
rawDataList <- getMatrix(pbmc, "rawData")

withNewH5Copy <- function(fun) {
    ctrlpath.orig <- system.file("extdata/ctrl.h5", package = "rliger")
    stimpath.orig <- system.file("extdata/stim.h5", package = "rliger")
    if (!file.exists(ctrlpath.orig))
        stop("Cannot find original h5 file at: ", ctrlpath.orig)
    # if (file.exists("ctrltest.h5")) file.remove("ctrltest.h5")
    # if (file.exists("stimtest.h5")) file.remove("stimtest.h5")
    # pwd <- getwd()
    # # Temp setting for GitHub Actions
    # fsep <- ifelse(Sys.info()["sysname"] == "Windows", "\\", "/")
    # if (Sys.info()["sysname"] == "Windows") {
    #     pwd <- file.path("C:\\Users", Sys.info()["user"], "Documents", fsep = fsep)
    # }

    # ctrlpath <- file.path(pwd, "ctrltest.h5", fsep = fsep)
    # stimpath <- file.path(pwd, "stimtest.h5", fsep = fsep)
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

    fun(list(ctrl = ctrlpath, stim = stimpath))

    if (file.exists(ctrlpath)) unlink(ctrlpath)
    if (file.exists(stimpath)) unlink(stimpath)
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
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# optimizeALS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("iNMF")
test_that("iNMF - in-memory", {
    skip_if_not_installed("RcppPlanc")
    pbmc <- normalize(pbmc)
    pbmc <- selectGenes(pbmc)

    expect_error(pbmc <- runINMF(pbmc),
                 "Scaled data not available. ")
    pbmc <- scaleNotCenter(pbmc)
    expect_error(pbmc <- runINMF(pbmc, k = 5000),
                 "Number of factors")
    pbmc <- runIntegration(pbmc, k = 10, nIteration = 2)
    expect_no_error(.checkValidFactorResult(pbmc))

    pbmc@datasets$ctrl@scaleData <- as.matrix(pbmc@datasets$ctrl@scaleData)
    expect_error(runINMF(pbmc),
                 "Currently the scaledData of all datasets have to be of the same")
})

# test_that("optimizeALS - on-disk", {
#     withNewH5Copy(
#         function(rawList) {
#             pbmc <- createLiger(rawList)
#             pbmc <- process(pbmc)
#             expect_warning(optimizeALS(pbmc, k = 10, readH5 = "auto",
#                                        maxIter = 1),
#                            "Automatically reading H5 based ")
#             expect_error(optimizeALS(pbmc, k = 10, readH5 = "hello"),
#                          "Can only set `readH5` to TRUE, FALSE")
#             expect_error(optimizeALS(pbmc, k = 10, readH5 = FALSE),
#                          "H5 based dataset detected while")
#             pbmc <- optimizeALS(pbmc, k = 10, readH5 = TRUE, maxIter = 2)
#             expect_no_error(.checkValidFactorResult(pbmc))
#         }
#     )
# })

test_that("UINMF", {
    skip_if_not_installed("RcppPlanc")
    # Need to fake the situation because test dataset doesn't have real
    # unshared var feature
    pbmc <- normalize(pbmc)
    pbmc <- selectGenes(pbmc, 1.2) # Should found 22 shared var features
    allFeature <- union(rownames(pbmc@datasets$ctrl), rownames(pbmc@datasets$stim))
    unshare <- setdiff(allFeature, varFeatures(pbmc))
    expect_error(runUINMF(pbmc),
                 "Scaled data not available. ")
    pbmc <- scaleNotCenter(pbmc)
    expect_error(runUINMF(pbmc),
                 "No scaled data for unshared feature found. ")
    varUnsharedFeatures(pbmc, "ctrl") <- unshare[1:10]


    pbmc <- scaleNotCenter(pbmc)
    # This aims at testing nVarFeature < k but we have enough unshared var feature
    pbmc <- runUINMF(pbmc, k = 30, nIteration = 2, nRandomStarts = 2)
    expect_no_error(rliger:::.checkValidFactorResult(pbmc))

    varUnsharedFeatures(pbmc, "stim") <- unshare[11:20]
    pbmc <- scaleNotCenter(pbmc, useDatasets = "stim")
    # Normal case
    expect_no_error(pbmc <- runUINMF(pbmc, k = 10, nIteration = 2))
})

test_that("Optimize new parameters", {
    skip_if_not_installed("RcppPlanc")
    pbmc <- process(pbmc)
    pbmc <- runOnlineINMF(pbmc, k = 10, minibatchSize = 100)
    pbmc0 <- optimizeNewK(pbmc, kNew = 10, nIteration = 2)
    expect_equal(pbmc@W, pbmc0@W)
    pbmc1 <- optimizeNewK(pbmc, kNew = 15, nIteration = 2)
    expect_no_error(.checkValidFactorResult(pbmc1))
    pbmc2 <- optimizeNewK(pbmc, kNew = 8, nIteration = 2)
    expect_no_error(.checkValidFactorResult(pbmc2))
    pbmc3 <- optimizeNewLambda(pbmc, lambdaNew = 5.5, nIteration = 2)
    expect_no_error(.checkValidFactorResult(pbmc3))
    expect_message(optimizeNewLambda(pbmc, lambdaNew = 4, nIteration = 2),
                  "New lambda less than current lambda")
    pbmc4 <- optimizeSubset(pbmc, cellIdx = sort(sample(ncol(pbmc), 300)),
                            nIteration = 2, scaleDatasets = names(pbmc))
    expect_no_error(.checkValidFactorResult(pbmc4))
    expect_equal(ncol(pbmc4), 300)

    ctrl2 <- rawData(dataset(pbmc, "ctrl"))
    ctrl2@x <- ctrl2@x + 1
    colnames(ctrl2) <- paste0(colnames(ctrl2), 2)
    pbmc5 <- optimizeNewData(pbmc, dataNew = list(ctrl2 = ctrl2),
                             useDatasets = "ctrl", merge = TRUE, nIteration = 2)
    expect_no_error(.checkValidFactorResult(pbmc5))
    expect_equal(ncol(pbmc5), 900)

    pbmc6 <- optimizeNewData(pbmc, dataNew = list(ctrl2 = ctrl2),
                             useDatasets = "ctrl", merge = FALSE, nIteration = 2)
    expect_no_error(.checkValidFactorResult(pbmc5))
    expect_equal(ncol(pbmc5), 900)
})

test_that("Online iNMF - in-memory", {
    skip_if_not_installed("RcppPlanc")
    expect_error(runOnlineINMF(pbmc, k = 20, minibatchSize = 100),
                 "Scaled data not available. ")
    pbmc <- process(pbmc)
    # Scenario 1
    pbmc <- runOnlineINMF(pbmc, k = 20, minibatchSize = 100)
    expect_no_error(.checkValidFactorResult(pbmc))
    # Scenario 2
    ctrl2 <- rawData(dataset(pbmc, "ctrl"))
    ctrl2@x <- ctrl2@x + 1
    colnames(ctrl2) <- paste0(colnames(ctrl2), 2)
    ctrlA <- getMatrix(pbmc, "A", "ctrl")
    pbmc@datasets$ctrl@A <- NULL
    expect_error(runOnlineINMF(pbmc, k = 20, newDatasets = list(ctrl2 = ctrl2),
                               minibatchSize = 100),
                 "Cannot find complete online iNMF result")
    pbmc@datasets$ctrl@A <- ctrlA
    expect_error(runOnlineINMF(pbmc, k = 20, newDatasets = list(ctrl = ctrl2),
                               minibatchSize = 100),
                 "Names of `newDatasets` overlap")
    expect_error(runOnlineINMF(pbmc, k = 20, newDatasets = list(ctrl2),
                               minibatchSize = 100),
                 "The list of new datasets must be named")
    expect_error(runOnlineINMF(pbmc, k = 20, newDatasets = list(ctrl2 = 1),
                               minibatchSize = 100),
                 "Cannot interpret")
    pbmc2 <- runOnlineINMF(pbmc, k = 20, newDatasets = list(ctrl2 = ctrl2),
                           minibatchSize = 100)
    expect_no_error(.checkValidFactorResult(pbmc2))
    W <- getMatrix(pbmc, "W")
    V <- getMatrix(pbmc, "V")
    A <- getMatrix(pbmc, "A")
    B <- getMatrix(pbmc, "B")
    pbmc2 <- runOnlineINMF(pbmc, k = 20, newDatasets = list(ctrl2 = ctrl2),
                           minibatchSize = 100, WInit = W, VInit = V,
                           AInit = A, BInit = B)
    expect_no_error(.checkValidFactorResult(pbmc2))
    # Scenario 3
    pbmc3 <- runOnlineINMF(pbmc, k = 20, newDatasets = list(ctrl2 = ctrl2),
                           minibatchSize = 100, projection = TRUE, WInit = W)
})

test_that("quantileNorm", {
    skip_if_not_installed("RcppPlanc")
    pbmc <- process(pbmc)
    pbmc <- runOnlineINMF(pbmc, k = 20, minibatchSize = 100)

    # For quantileNorm,liger method
    expect_error(quantileNorm(pbmc, reference = names(pbmc)),
                 "Should specify only one reference dataset.")

    pbmc2 <- quantileNorm(pbmc)
    expect_equal(dim(getMatrix(pbmc2, "H.norm")), c(ncol(pbmc), 20))

    pbmc2 <- quantileNorm(pbmc, reference = "ctrl")
    expect_equal(dim(getMatrix(pbmc2, "H.norm")), c(ncol(pbmc), 20))

    # For quantileNorm,list method
    # expect_error(quantileNorm(list("hello")),
    #              "`object` should be a named list of matrices")
    # expect_error(quantileNorm(list(foo = "bar")),
    #              "All values in `object` must be a matrix")
    # mat1 <- matrix(runif(20), 4, 5)
    # mat2 <- matrix(runif(20), 5, 4)
    # expect_error(quantileNorm(list(a = mat1, b = mat2)),
    #              "All matrices must have the same number of rows")
    # Hs <- getMatrix(pbmc, "H")
    # expect_error(quantileNorm(Hs, reference = "Hey"),
    #              "Should specify one existing dataset")
    # expect_error(quantileNorm(Hs, reference = 114514),
    #              "Should specify one dataset within the range. ")
    # expect_error(quantileNorm(Hs, reference = c(TRUE, FALSE, TRUE)),
    #              "logical `reference` wrong length or")
    # expect_error(quantileNorm(Hs, reference = dataset(pbmc, "ctrl")),
    #              "Unable to understand `reference`.")
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Seurat wrapper for everything
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("Seurat iNMF wrapper")
test_that("Seurat wrapper", {
    skip_if_not_installed("RcppPlanc")
    skip_if_not_installed("Seurat")
    skip_if_not_installed("SeuratObject")
    seu <- ligerToSeurat(pbmc)
    seu <- normalize(seu)
    seu <- selectGenes(seu)
    seu <- scaleNotCenter(seu)
    seu1 <- runIntegration(seu)
    expect_in("inmf", SeuratObject::Reductions(seu1))

    x1 <- seu@assays$RNA@layers$ligerScaleData.ctrl@x[1]
    seu@assays$RNA@layers$ligerScaleData.ctrl@x[1] <- -1
    expect_error(runINMF(seu),
                 "Negative data encountered")
    expect_error(runOnlineINMF(seu, k = 10, minibatchSize = 100),
                 "Negative data encountered")
    seu@assays$RNA@layers$ligerScaleData.ctrl@x[1] <- x1
    seu$orig.ident <- as.character(seu$orig.ident)
    seu <- runINMF(seu)
    seu <- runOnlineINMF(seu, k = 10, minibatchSize = 100)
    expect_in("inmf", SeuratObject::Reductions(seu))
    expect_in("onlineINMF", SeuratObject::Reductions(seu))

    expect_error(quantileNorm(seu, reduction = "orig.ident"),
                 "Specified `reduction` does not points to a")
    seu <- quantileNorm(seu, reduction = "inmf")
    expect_in("inmfNorm", SeuratObject::Reductions(seu))

    expect_error(quantileNorm(seu, reference = "hello"),
                 "Should specify one existing dataset")
    expect_error(quantileNorm(seu, reference = 114514),
                 "Should specify one existing dataset as reference")
    expect_error(quantileNorm(seu, reference = c(TRUE, FALSE, TRUE)),
                 "Should specify one existing dataset as reference")
})
