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
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# optimizeALS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

context("optimizeALS")
test_that("optimizeALS - in-memory", {
    pbmc <- normalize(pbmc)
    pbmc <- selectGenes(pbmc)
    expect_error(pbmc <- optimizeALS(pbmc),
                 "Scaled data not available. ")
    pbmc <- scaleNotCenter(pbmc)
    expect_error(pbmc <- optimizeALS(pbmc, k = 5000),
                 "Select k lower than the number of variable genes: ")
    pbmc <- optimizeALS(pbmc, k = 10, max.iters = 2)
    expect_no_error(.checkValidFactorResult(pbmc))
})

test_that("optimizeALS - on-disk", {
    withNewH5Copy(
        function(rawList) {
            pbmc <- createLiger(rawList)
            pbmc <- process(pbmc)
            expect_warning(optimizeALS(pbmc, k = 10, readH5 = "auto",
                                       maxIter = 1),
                           "Automatically reading H5 based ")
            expect_error(optimizeALS(pbmc, k = 10, readH5 = "hello"),
                         "Can only set `readH5` to TRUE, FALSE")
            expect_error(optimizeALS(pbmc, k = 10, readH5 = FALSE),
                         "H5 based dataset detected while")
            pbmc <- optimizeALS(pbmc, k = 10, readH5 = TRUE, maxIter = 2)
            expect_no_error(.checkValidFactorResult(pbmc))
        }
    )
})

test_that("optimizeUANLS", {
    # Need to fake the situation because test dataset doesn't have real
    # unshared var feature
    pbmc <- normalize(pbmc)
    pbmc <- selectGenes(pbmc, 0.8)
    allFeature <- union(rownames(pbmc@datasets$ctrl), rownames(pbmc@datasets$stim))
    unshare <- setdiff(allFeature, varFeatures(pbmc))
    pbmc@datasets$ctrl@varUnsharedFeatures <- unshare[1:10]
    pbmc@datasets$stim@varUnsharedFeatures <- unshare[11:20]
    pbmc <- scaleNotCenter(pbmc)
    pbmc <- optimizeALS(pbmc, k = 10, maxIter = 2, useUnshared = TRUE)
    expect_no_error(.checkValidFactorResult(pbmc))
})

test_that("Optimize new parameters", {
    pbmc <- process(pbmc)
    pbmc <- online_iNMF(pbmc, k = 10, miniBatch_size = 100)
    pbmc0 <- optimizeNewK(pbmc, k.new = 10, max.iters = 2)
    expect_equal(pbmc@W, pbmc0@W)
    pbmc1 <- optimizeNewK(pbmc, k.new = 15, max.iters = 2)
    expect_no_error(.checkValidFactorResult(pbmc1))
    pbmc2 <- optimizeNewK(pbmc, k.new = 8, max.iters = 2)
    expect_no_error(.checkValidFactorResult(pbmc2))
    pbmc3 <- optimizeNewLambda(pbmc, new.lambda = 5.5, max.iters = 2)
    expect_no_error(.checkValidFactorResult(pbmc3))
    expect_message(optimizeNewLambda(pbmc, new.lambda = 4, max.iters = 2),
                  "New lambda less than current lambda")
    pbmc4 <- optimizeSubset(pbmc, sort(sample(ncol(pbmc), 300)), max.iters = 2,
                            datasets.scale = names(pbmc))
    expect_no_error(.checkValidFactorResult(pbmc4))
    expect_equal(ncol(pbmc4), 300)

    ctrl2 <- rawData(dataset(pbmc, "ctrl"))
    ctrl2@x <- ctrl2@x + 1
    colnames(ctrl2) <- paste0(colnames(ctrl2), 2)
    pbmc5 <- optimizeNewData(pbmc, new.data = list(ctrl2 = ctrl2),
                             which.datasets = "ctrl",
                             add.to.existing = TRUE, max.iters = 2)
    expect_no_error(.checkValidFactorResult(pbmc5))
    expect_equal(ncol(pbmc5), 900)

    pbmc6 <- optimizeNewData(pbmc, new.data = list(ctrl2 = ctrl2),
                             which.datasets = "ctrl",
                             add.to.existing = FALSE, max.iters = 2)
    expect_no_error(.checkValidFactorResult(pbmc5))
    expect_equal(ncol(pbmc5), 900)
})

test_that("Online iNMF - in-memory", {
    pbmc <- process(pbmc)
    # Scenario 1
    pbmc <- online_iNMF(pbmc, k = 20, miniBatch_size = 100)
    expect_no_error(.checkValidFactorResult(pbmc))
    # Scenario 2
    ctrl2 <- rawData(dataset(pbmc, "ctrl"))
    ctrl2@x <- ctrl2@x + 1
    colnames(ctrl2) <- paste0(colnames(ctrl2), 2)
    pbmc2 <- online_iNMF(pbmc, k = 20, X_new = list(ctrl2 = ctrl2),
                         miniBatch_size = 100)
    expect_no_error(.checkValidFactorResult(pbmc2))
    W <- getMatrix(pbmc, "W")
    V <- getMatrix(pbmc, "V")
    A <- getMatrix(pbmc, "A")
    B <- getMatrix(pbmc, "B")
    pbmc2 <- online_iNMF(pbmc, k = 20, X_new = list(ctrl2 = ctrl2),
                         miniBatch_size = 100, W.init = W, V.init = V,
                         A.init = A, B.init = B)
    expect_no_error(.checkValidFactorResult(pbmc2))
    # Scenario 3
    pbmc3 <- online_iNMF(pbmc, k = 20, X_new = list(ctrl2 = ctrl2),
                         miniBatch_size = 100, projection = TRUE,
                         W.init = W)
    expect_no_error(.checkValidFactorResult(pbmc3))
})

test_that("quantileNorm", {
    pbmc <- process(pbmc)
    pbmc <- online_iNMF(pbmc, k = 20, miniBatch_size = 100)

    # For quantileNorm,liger method
    expect_error(quantileNorm(pbmc, reference = names(pbmc)),
                 "Should specify only one reference dataset.")

    pbmc2 <- quantileNorm(pbmc)
    expect_equal(dim(getMatrix(pbmc2, "H.norm")), c(ncol(pbmc), 20))

    pbmc2 <- quantileNorm(pbmc, reference = "ctrl")
    expect_equal(dim(getMatrix(pbmc2, "H.norm")), c(ncol(pbmc), 20))

    # For quantileNorm,list method
    expect_error(quantileNorm(list("hello")),
                 "`object` should be a named list of matrices")
    expect_error(quantileNorm(list(foo = "bar")),
                 "All values in `object` must be a matrix")
    mat1 <- matrix(runif(20), 4, 5)
    mat2 <- matrix(runif(20), 5, 4)
    expect_error(quantileNorm(list(a = mat1, b = mat2)),
                 "All matrices must have the same number of rows")
    Hs <- getMatrix(pbmc, "H")
    expect_error(quantileNorm(Hs, reference = "Hey"),
                 "Should specify one existing dataset")
    expect_error(quantileNorm(Hs, reference = 114514),
                 "Should specify one dataset within the range. ")
    expect_error(quantileNorm(Hs, reference = c(TRUE, FALSE, TRUE)),
                 "logical `reference` wrong length or")
    expect_error(quantileNorm(Hs, reference = dataset(pbmc, "ctrl")),
                 "Unable to understand `reference`.")
    expect_warning(quantile_norm(pbmc), "was deprecated in")
})
