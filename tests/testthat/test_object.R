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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# liger object creation
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
            expect_error(createLiger(rawList, formatType = "Hello"),
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# liger object methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
context("liger object S3/S4 methods")

test_that("liger S3/S4 methods", {
    pbmc <- process(pbmc)
    expect_output(show(pbmc), "An object of class liger with 600 cells")
    expect_equal(dim(pbmc), c(NA, 600))
    expect_null(rownames(pbmc))

    bc <- colnames(pbmc)
    expect_identical(bc, rownames(cellMeta(pbmc)))
    expect_no_error(colnames(pbmc) <- bc)

    expect_equal(dim(pbmc[varFeatures(pbmc)[1:5],]), c(NA, 600))
    expect_equal(dim(pbmc[varFeatures(pbmc)[1:5], pbmc$dataset == "ctrl"]),
                 c(NA, 300))
    expect_equal(dim(pbmc[,pbmc$dataset == "ctrl"]), c(NA, 300))

    ldList <- datasets(pbmc)
    expect_is(ldList, "list")
    expect_no_error(datasets(pbmc)[["ctrl"]] <- ldList[["ctrl"]])
    expect_no_error(datasets(pbmc, check = FALSE)[["ctrl"]] <- ldList[["ctrl"]])

    expect_is(dataset(pbmc), "ligerDataset")
    expect_is(dataset(pbmc, "ctrl"), "ligerDataset")
    expect_is(dataset(pbmc, NULL), "ligerDataset")
    expect_is(dataset(pbmc, 2), "ligerDataset")

    dataset(pbmc, "ctrl") <- ldList$ctrl
    expect_equal(names(pbmc), c("stim", "ctrl"))
    dataset(pbmc, "stim2") <- rawData(dataset(pbmc, "stim"))
    expect_equal(names(pbmc), c("stim", "ctrl", "stim2"))
    dataset(pbmc, "stim2") <- NULL
    expect_equal(length(pbmc), 2)

    names(pbmc) <- c("STIM", "CTRL")
    expect_identical(levels(pbmc$dataset), c("STIM", "CTRL"))

    meta <- cellMeta(pbmc)
    expect_is(meta, "DFrame")
    expect_null(cellMeta(pbmc, NULL))
    expect_is(cellMeta(pbmc, "dataset"), "factor")
    expect_warning(cellMeta(pbmc, "UMAP.1"),
                   "Specified variables from cellMeta not found: UMAP.1")
    expect_is(cellMeta(pbmc, "UMAP.1", cellIdx = 1:500, as.data.frame = TRUE),
              "numeric")
    expect_is(pbmc[["nUMI"]], "numeric")
    expect_is(pbmc$mito, "numeric")

    expect_no_error(cellMeta(pbmc) <- meta)
    cellMeta(pbmc, "newZeros") <- 0
    expect_true("newZeros" %in% colnames(cellMeta(pbmc)))
    pbmc[["newOnes"]] <- 1
    expect_true("newOnes" %in% colnames(cellMeta(pbmc)))
    pbmc$newTwos <- 2
    expect_true("newTwos" %in% colnames(cellMeta(pbmc)))

    expect_is(varFeatures(pbmc), "character")
    expect_no_error(varFeatures(pbmc) <- varFeatures(pbmc))

    expect_is(c(pbmc, pbmc), "liger")
    expect_is(fortify(pbmc), "data.frame")
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ligerDataset object creation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
# ligerDataset object methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("ligerDataset methods", {
    pbmc <- process(pbmc)
    expect_false(isH5Liger(pbmc))
    ctrl <- dataset(pbmc, "ctrl")
    expect_false(isH5Liger(ctrl))
    expect_warning(isH5Liger("hi"), "Given object is not of liger")

    expect_identical(modalOf(ctrl), "default")
    expect_identical(modalOf(pbmc), c(ctrl = "default", stim = "default"))
    expect_output(show(ctrl), "An object of class ligerDataset with 300 cells")

    expect_equal(dim(ctrl), c(266, 300))
    expect_is(dimnames(ctrl), "list")
    expect_identical(rownames(ctrl)[1:3], c("ISG15", "ID3", "RPL11"))
    expect_identical(colnames(ctrl)[1:2],
                     c("ctrl_AAACATACCTCGCT.1", "ctrl_AAACGGCTCTTCGC.1"))
    expect_no_error(colnames(ctrl) <- colnames(ctrl))

    expect_equal(dim(ctrl[1:10,]), c(10, 300))
    expect_equal(dim(ctrl[,1:10]), c(266, 10))
    expect_equal(dim(ctrl[1:10, 1:10]), c(10, 10))

    expect_is(rawData(ctrl), "dgCMatrix")
    expect_no_error(rawData(ctrl) <- rawData(ctrl))
    expect_error(rawData(ctrl) <- t(rawData(ctrl)), "invalid class")

    expect_is(normData(ctrl), "dgCMatrix")
    expect_is(scaleData(ctrl), "matrix")
    expect_is(scaleData(pbmc, 2), "matrix")
    expect_null(scaleUnsharedData(ctrl))
    expect_null(scaleUnsharedData(pbmc, 2))
    expect_no_error(scaleUnsharedData(ctrl) <- scaleData(ctrl))

    expect_is(getMatrix(ctrl, "rawData"), "dgCMatrix")
    expect_is(getMatrix(pbmc, "W"), "matrix")
    expect_is(getMatrix(pbmc, "H.norm"), "matrix")
    expect_is(getMatrix(pbmc, "V"), "list")
    expect_is(getMatrix(pbmc, "H", dataset = 1, returnList = TRUE), "list")
    expect_is(getMatrix(pbmc, "H", dataset = 1, returnList = FALSE), "matrix")
    expect_is(getMatrix(pbmc, "H", dataset = 1:2), "list")

    expect_is(featureMeta(ctrl), "DFrame")
    expect_no_error(featureMeta(ctrl) <- featureMeta(ctrl))
    expect_no_error(featureMeta(ctrl) <- as.data.frame(featureMeta(ctrl)))

    stim <- dataset(pbmc, "stim")
    merged <- cbind(ctrl, stim)
    expect_equal(dim(merged), c(279, 600))
})

test_that("H5 ligerDataset methods", {
    withNewH5Copy(
        function(rawList) {
            pbmc <- createLiger(rawList)
            pbmc <- process(pbmc)
            expect_true(isH5Liger(pbmc))
            ctrl <- dataset(pbmc, "ctrl")
            expect_output(show(ctrl),
                          "An object of class ligerDataset with 300 cells")
            expect_true(isH5Liger(ctrl))
            h5file <- getH5File(pbmc, "ctrl")
            ctrl.h5 <- getH5File(ctrl)
            expect_is(ctrl.h5, "H5File")
            expect_is(getH5File(pbmc), "list")
            expect_identical(ctrl.h5, h5file)

            expect_no_error(rawData(ctrl) <- ctrl.h5[["matrix/data"]])
            expect_no_error(normData(ctrl) <- ctrl.h5[["normData"]])
            expect_no_error(scaleData(ctrl) <- ctrl.h5[["scaleData"]])
            expect_no_error(scaleUnsharedData(ctrl) <- ctrl.h5[["scaleData"]])
            expect_error(rawData(ctrl) <- matrix(1),
                         "Cannot replace slot with in-memory")
            expect_error(normData(ctrl) <- matrix(1),
                         "Cannot replace slot with in-memory")
            expect_error(scaleData(ctrl) <- matrix(1),
                         "Cannot replace slot with in-memory")
            expect_error(scaleUnsharedData(ctrl) <- matrix(1),
                         "Cannot replace slot with in-memory")

            expect_is(h5fileInfo(ctrl), "list")
            expect_identical(h5fileInfo(ctrl, "formatType"), "10X")
            expect_identical(h5fileInfo(ctrl, c("indicesName", "indptrName")),
                             list(indicesName = "matrix/indices",
                                  indptrName = "matrix/indptr"))
            expect_error(h5fileInfo(ctrl, c("indicesName", "hello")),
                         "Specified h5file info not found: hello")

            expect_error(h5fileInfo(ctrl, info = 1:2) <- "hey",
                         "`info` has to be a single character.")
            expect_error(h5fileInfo(ctrl, "indicesName") <- "hey",
                         "Specified info is invalid,")
            expect_no_error(h5fileInfo(ctrl, "barcodesName") <-
                                "matrix/barcodes")

            ctrl.h5$close()
            expect_warning(show(ctrl), "Link to HDF5 file fails.")
        }
    )
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# class conversion
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("as.liger methods", {
    # dgCMatrix
    ctrlRaw <- rawDataList$ctrl
    lig <- as.liger(ctrlRaw, sampleName = "ctrl")
    expect_is(lig, "liger")

    lig <- as.liger(ctrlRaw, sampleName = rep("ctrl", ncol(ctrlRaw)))
    expect_is(lig, "liger")
})

test_that("as.ligerDataset methods", {
    # ligerDataset
    ctrlLD <- dataset(pbmc, "ctrl")
    ld <- as.ligerDataset(ctrlLD)
    expect_is(ld, "ligerDataset")
    ld <- as.ligerDataset(ctrlLD, modal = "atac")
    expect_is(ld, "ligerATACDataset")
    expect_warning(ld <- as.ligerDataset(ld, modal = "rna"),
                   "Will remove information in the following slots when ")
    expect_is(ld, "ligerDataset")
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Importing data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("Importing data", {
    expect_error(importDataset("hello"),
                 "Requested dataset")
    obj <- importDataset("bmmc")
    expect_is(obj, "liger")
    expect_is(obj@datasets[[1]], "ligerDataset")
    expect_is(obj@datasets[[3]], "ligerATACDataset")
    expect_equal(ncol(obj), 16710)
    expect_warning(
        obj <- importDataset("bmmc"),
        "File already exists"
    )
})
