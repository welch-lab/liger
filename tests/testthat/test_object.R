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
    object <- runOnlineINMF(object, k = 10, minibatchSize = 100)
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
                 "`modal` has to be a length 1 or 2 object of class")
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
            expect_message(createLiger(rawList, formatType = "Hello"),
                         "Specified `formatType`")

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
    skip_if_not_installed("RcppPlanc")
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
    expect_message(cellMeta(pbmc, "UMAP.1"),
                   "Specified variables from cellMeta not found")
    expect_is(cellMeta(pbmc, "nUMI", cellIdx = 1:500, as.data.frame = TRUE),
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

    expect_error(expect_is(c(pbmc, pbmc), "liger"))
    expect_is(ggplot2::fortify(pbmc), "data.frame")

    expect_no_error(print(commands(pbmc, "normalize")))
    pbmc <- normalize(pbmc, scaleFactor = 10, log = TRUE)
    pbmc <- normalize(pbmc, scaleFactor = 100, log = TRUE)
    normCmds <- commands(pbmc, "normalize")
    expect_equal(commandDiff(pbmc, names(normCmds)[2], names(normCmds)[3]),
                 "Argument not identical: scaleFactor")
    expect_identical(commands(pbmc, names(normCmds)[3], "scaleFactor"),
                     c(scaleFactor = 100))
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ligerDataset object creation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("ligerDataset (in memory) object creation", {
    expect_error(createLigerDataset(),
                 "At least one of")

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
    expect_error(
        ld <- createLigerDataset(scaleData = scaledMat, featureMeta = featuremeta),
        "At least one of "
    )
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ligerDataset object methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("ligerDataset methods", {
    skip_if_not_installed("RcppPlanc")
    pbmc <- process(pbmc)
    expect_false(isH5Liger(pbmc))
    ctrl <- dataset(pbmc, "ctrl")
    expect_false(isH5Liger(ctrl))
    expect_message(isH5Liger("hi"), "Given object is not ")

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
    expect_is(scaleData(ctrl), "dgCMatrix")
    expect_is(scaleData(pbmc, 2), "dgCMatrix")
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
    expect_no_error(featureMeta(ctrl) <- rliger:::.DataFrame.as.data.frame(featureMeta(ctrl)))

    stim <- dataset(pbmc, "stim")
    merged <- cbind(ctrl, stim)
    expect_equal(dim(merged), c(279, 600))


    # ligerATACDataset related
    expect_error(rawPeak(pbmc, "stim"),
                 "unable to find an inherited")
    expect_error(rawPeak(pbmc, "stim") <- rawData(ctrl),
                 "unable to find an inherited")
    ctrl <- as.ligerDataset(ctrl, modal = "atac")
    pbmc@datasets$ctrl <- ctrl
    rawPeak(pbmc, "ctrl") <- rawData(ctrl)

    expect_error(normPeak(pbmc, "stim"),
                 "unable to find an inherited")
    expect_error(normPeak(pbmc, "stim") <- normData(stim),
                 "unable to find an inherited")
    normPeak(pbmc, "ctrl") <- normData(ctrl)
    expect_true(identical(normPeak(pbmc, "ctrl"), normData(ctrl, "ctrl")))
    expect_true(validObject(ctrl))
    # ligerSpatialDataset related

    expect_message(ctrl <- as.ligerDataset(ctrl, modal = "spatial"),
                   "Will remove information in the following slots when ")
    pbmc@datasets$ctrl <- ctrl
    coords <- matrix(rnorm(300*2), 300, 2)
    rownames(coords) <- colnames(ctrl)
    colnames(coords) <- c("x", "y")
    expect_error(coordinate(pbmc, "stim"),
                 "unable to find an inherited")
    expect_error(coordinate(pbmc, "stim") <- coords,
                 "unable to find an inherited")
    coordinate(pbmc, "ctrl") <- coords
    expect_true(identical(coordinate(pbmc, "ctrl"), coords))
    expect_true(validObject(ctrl))

    coords <- matrix(rnorm(300*3), 300, 3)
    expect_message(coordinate(ctrl) <- coords,
                   "No rownames with given spatial coordinate")
    coords <- matrix(rnorm(300*4), 300, 4)
    rownames(coords) <- colnames(ctrl)
    expect_error(coordinate(ctrl) <- coords,
                 "More than 3 dimensions for the coordinates")

    coords <- matrix(rnorm(300*2), 300, 2)
    rownames(coords) <- c(colnames(ctrl)[1:299], "hello")
    colnames(coords) <- c("x", "y")
    ctrl@coordinate <- coords
    expect_error(validObject(ctrl), "Inconsistant cell identifiers")
    expect_message(coordinate(ctrl) <- coords,
                   "NA generated for missing cells")
    # ligerMethDataset related
    expect_message(ctrl <- as.ligerDataset(ctrl, modal = "meth"),
                   "Will remove information in the following slots when ")
    expect_no_error(validObject(ctrl))
})

test_that("H5 ligerDataset methods", {
    skip_if_not_installed("RcppPlanc")
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
            expect_no_error(scaleData(ctrl) <- ctrl.h5[["scaleDataSparse"]])
            expect_no_error(scaleUnsharedData(ctrl) <- ctrl.h5[["scaleDataSparse"]])
            expect_error(rawData(ctrl) <- matrix(1),
                         "Cannot replace slot with in-memory")
            expect_error(normData(ctrl) <- matrix(1),
                         "Cannot replace slot with in-memory")
            expect_error(scaleData(ctrl) <- matrix(1),
                         "Cannot replace slot with in-memory")
            expect_error(scaleUnsharedData(ctrl) <- matrix(1),
                         "Cannot replace slot with in-memory")

            expect_is(h5fileInfo(ctrl), "list")
            expect_identical(h5fileInfo(ctrl, "formatType"), "10x")
            expect_identical(h5fileInfo(ctrl, c("indicesName", "indptrName")),
                             list(indicesName = "matrix/indices",
                                  indptrName = "matrix/indptr"))
            expect_error(h5fileInfo(ctrl, c("indicesName", "hello")),
                         "Specified `info` not found")

            expect_error(h5fileInfo(ctrl, info = 1:2) <- "hey",
                         "`info` has to be a single character.")
            expect_error(h5fileInfo(ctrl, "indicesName") <- "hey",
                         "Specified `info`")
            expect_no_error(h5fileInfo(ctrl, "barcodesName") <-
                                "matrix/barcodes")

            # ctrl.h5$close()
            closeAllH5(ctrl)
            expect_message(show(ctrl), "Link to HDF5 file fails.")
        }
    )
})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# class conversion
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("as.liger methods", {
    # dgCMatrix
    ctrlRaw <- rawDataList$ctrl
    lig <- as.liger(ctrlRaw)
    expect_equal(names(lig), "sample")

    lig <- as.liger(ctrlRaw, datasetVar = "ctrl")
    expect_equal(names(lig), "ctrl")

    lig <- as.liger(ctrlRaw, datasetVar = c(rep("ctrl", 150), rep("stim", 150)))
    expect_true(identical(names(lig), c("ctrl", "stim")))

    # SCE
    if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = ctrlRaw),
            colData = data.frame(dataset = factor(rep(c("a", "b"), each = 150)))
        )
        sce$useless <- 1
        expect_message(lig <- as.liger(sce))
        expect_equal(names(lig), "SCE")

        expect_message(lig <- as.liger(sce, datasetVar = "dataset"))
        expect_true(all.equal(sapply(datasets(lig), ncol), c(a = 150, b = 150)))
    }

    skip_if_not_installed("Seurat")
    skip_if_not_installed("SeuratObject")
    # Seurat
    seu <- SeuratObject::CreateSeuratObject(
        ctrlRaw,
        meta.data = data.frame(orig.ident = factor(rep(c("a", "b"), each = 150)),
                               nUMI = 0,
                               row.names = colnames(ctrlRaw))
    )

    seu <- Seurat::NormalizeData(seu) %>%
        Seurat::FindVariableFeatures() %>%
        Seurat::ScaleData() %>%
        Seurat::RunPCA()
    expect_message(lig <- as.liger(seu))
    expect_true(all.equal(sapply(datasets(lig), ncol), c(a = 150, b = 150)))

    expect_in("pca", names(dimReds(lig)))
})

test_that("as.ligerDataset methods", {
    # ligerDataset
    ctrlLD <- dataset(pbmc, "ctrl")
    ld <- as.ligerDataset(ctrlLD)
    expect_is(ld, "ligerDataset")
    ld <- as.ligerDataset(ctrlLD, modal = "atac")
    expect_is(ld, "ligerATACDataset")
    expect_message(ld <- as.ligerDataset(ld, modal = "rna"),
                   "Will remove information in the following slots when ")
    expect_is(ld, "ligerDataset")

    # matrix
    mat <- matrix(rnorm(26*26), 26, 26, dimnames = list(letters, letters))
    ld <- as.ligerDataset(mat, normData = mat, scaleData = mat,
                          featureMeta = data.frame(id = 1:26, row.names = letters))
    expect_true(all.equal(rownames(ld), letters))
    expect_true(all.equal(colnames(ld), letters))

    if (requireNamespace("Seurat", quietly = TRUE)) {
        # Seurat
        seu <- SeuratObject::CreateSeuratObject(rawData(ctrlLD))
        ld <- as.ligerDataset(seu)
        expect_is(ld, "ligerDataset")
    }

    # SCE
    if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = rawData(ctrlLD))
        )
        ld <- as.ligerDataset(sce)
        expect_is(ld, "ligerDataset")
    }
})

test_that("ligerToSeurat", {
    skip_if_not_installed("Seurat")
    skip_if_not_installed("SeuratObject")
    seu <- ligerToSeurat(pbmc)
    expect_equal(SeuratObject::Assays(seu), "RNA")

    pbmc@datasets$stim <- as.ligerDataset(pbmc@datasets$stim, modal = "atac")
    pbmc <- normalize(pbmc, useDatasets = "ctrl")
    seu <- ligerToSeurat(pbmc)
    expect_equal(SeuratObject::Assays(seu), "LIGER")
    expect_true(all.equal(SeuratObject::Layers(seu),
                          c("counts.ctrl", "counts.stim", "ligerNormData.ctrl")))

    expect_error(seu <- ligerToSeurat(pbmcPlot), "rawData not found")

    rawData(pbmcPlot, "ctrl") <- rawData(pbmc, "ctrl")
    rawData(pbmcPlot, "stim") <- rawData(pbmc, "stim")
    seu <- ligerToSeurat(pbmcPlot, identByDataset = TRUE)
})
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Importing data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# test_that("Importing data", {
#     obj <- importBMMC()
#     expect_is(obj, "liger")
#     expect_is(obj@datasets[[1]], "ligerDataset")
#     expect_is(obj@datasets[[3]], "ligerATACDataset")
#     expect_equal(ncol(obj), 16710)
#     expect_warning(
#         obj <- importBMMC(),
#         "File already exists"
#     )
#     unlink("liger_BMMC_rna_D1T1.rds")
#     unlink("liger_BMMC_rna_D1T2.rds")
#     unlink("liger_BMMC_atac_D5T1.rds")
#     unlink("liger_BMMC_atac_D5T1_peak.rds")
# })
