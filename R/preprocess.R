#################################### qc ########################################

#' General QC for liger object
#' @description Calculate number of UMIs, number of detected features and
#' percentage of feature subset (e.g. mito) expression per cell.
#' @param object \linkS4class{liger} object with \code{rawData} available in
#' each \linkS4class{ligerDataset} embedded
#' @param mito,ribo,hemo Whether to calculate the expression percentage of
#' mitochondrial, ribosomal or hemoglobin genes, respectively. Default
#' \code{TRUE}.
#' @param features Feature names matching the feature subsets that users want to
#' calculate the expression percentage with. A vector for a single subset, or a
#' named list for multiple subset. Default \code{NULL}.
#' @param pattern Regex patterns for matching the feature subsets that users
#' want to calculate the expression percentage with. A vector for a single
#' subset, or a named list for multiple subset. Default \code{NULL}.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to be included for QC. Default
#' \code{NULL} performs QC on all datasets.
#' @param chunkSize Integer number of cells to include in a chunk when working
#' on HDF5 based dataset. Default \code{1000}
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @return Updated \code{object} with \code{nUMI}, \code{nGene} updated
#' in \code{cellMeta(object)}, as well as expression percentage value for each
#' feature subset.
#' @export
#' @examples
#' pbmc <- runGeneralQC(pbmc)
runGeneralQC <- function(
        object,
        mito = TRUE,
        ribo = TRUE,
        hemo = TRUE,
        features = NULL,
        pattern = NULL,
        useDatasets = NULL,
        chunkSize = 1000,
        verbose = getOption("ligerVerbose", TRUE)
) {
    .checkObjVersion(object)
    useDatasets <- .checkUseDatasets(object, useDatasets)
    # Process the the two arguments all into one named list of feature names
    # before exactly calculate the percentage
    featureSubsets <- list()
    allFeatures <- unique(unlist(lapply(datasets(object), rownames),
                                 use.names = FALSE))

    # Work on the presets
    if (isTRUE(mito))
        featureSubsets$mito <- grep("^MT-", allFeatures, value = TRUE)
    if (isTRUE(ribo))
        featureSubsets$ribo <- grep("^RP[SL]", allFeatures, value = TRUE)
    if (isTRUE(hemo))
        featureSubsets$hemo <- grep("^HB[^(P)]", allFeatures, value = TRUE)

    # Then process the user specified gene sets
    if (!is.null(features)) {
        if (is.list(features)) {
            featureSubsets <- c(featureSubsets, features)
        } else if (is.vector(features)) {
            featureSubsets[["featureSubset_name"]] <- features
        }
    }
    if (!is.null(pattern)) {
        if (is.list(pattern)) {
            pattern <- lapply(pattern, function(x) {
                grep(x, allFeatures, value = TRUE)
            })
            featureSubsets <- c(featureSubsets, pattern)
        } else if (is.vector(pattern)) {
            pattern <- grep(pattern, allFeatures, value = TRUE)
            featureSubsets[["featureSubset_pattern"]] <- pattern
        }
    }

    # Start calculation on each dataset
    newResultNames <- c("nUMI", "nGene", names(featureSubsets))

    for (d in useDatasets) {
        ld <- dataset(object, d)
        if (isTRUE(verbose))
            cliID <- cli::cli_process_start("calculating QC for dataset {.val {d}}")
        if (isH5Liger(ld))
            results <- runGeneralQC.h5(
                ld,
                featureSubsets = featureSubsets,
                chunkSize = chunkSize,
                verbose = verbose
            )
        else
            results <- runGeneralQC.Matrix(
                ld,
                featureSubsets = featureSubsets,
                verbose = verbose
            )
        object@cellMeta[object$dataset == d, newResultNames] <- results$cell
        featureMeta(ld, check = FALSE)$nCell <- results$feature
        datasets(object, check = FALSE)[[d]] <- ld
        if (isTRUE(verbose)) cli::cli_process_done(id = cliID)
    }

    return(object)
}

#' Calculate general QC on H5 based ligerDataset object
#' @param object ligerDataset object
#' @param featureSubsets Named list passed from \code{runGeneralQC}
#' @param chunkSize Integer
#' @return data.frame
#' @noRd
runGeneralQC.h5 <- function(
        object,
        featureSubsets = NULL,
        chunkSize = 1000,
        verbose = getOption("ligerVerbose", TRUE)) {
    allFeatures <- rownames(object)
    # Initialize results
    cell <- data.frame(row.names = colnames(object))
    cell$nUMI <- 0
    cell$nGene <- 0
    for (i in names(featureSubsets)) {
        cell[[i]] <- 0
    }
    nCell <- rep(0, nrow(object))
    rowIndices <- lapply(featureSubsets, function(x) allFeatures %in% x)

    # Calculate in only one iteration
    H5Apply(
        object,
        init = list(cell = cell, feature = nCell),
        useData = "rawData",
        chunkSize = chunkSize,
        verbose = verbose,
        FUN = function(chunk, sparseXIdx, cellIdx, values) {
            nUMI <- colSums(chunk)
            values$cell$nUMI[cellIdx] <- nUMI
            nonzero <- methods::as(chunk, "lMatrix")
            values$cell$nGene[cellIdx] <- colSums(nonzero)
            for (fs in names(rowIndices)) {
                values$cell[[fs]][cellIdx] <-
                    colSums(chunk[rowIndices[[fs]], , drop = FALSE]) / nUMI *
                    100
            }
            values$feature <- values$feature + Matrix::rowSums(nonzero)
            return(values)
        }
    )
}

#' Calculate general QC on in memory matrix based ligerDataset object
#' @param object ligerDataset object
#' @param featureSubsets  Named list passed from \code{runGeneralQC}
#' @return data.frame
#' @noRd
runGeneralQC.Matrix <- function(
        object,
        featureSubsets = NULL,
        verbose = getOption("ligerVerbose", TRUE)) {
    nUMI <- Matrix::colSums(rawData(object))
    # Instead of using `nonzero <- rawData > 0` which generates dense logical
    # matrix, keep it sparse with 1 for TRUE
    # nonzero <- rawData(object)
    # nonzero@x <- rep(1, length(nonzero@x))
    nonzero <- methods::as(rawData(object), "lMatrix")
    nGene <- Matrix::colSums(nonzero)
    nCell <- Matrix::rowSums(nonzero)
    results <- data.frame(nUMI = nUMI, nGene = nGene,
                          row.names = colnames(object))
    if (length(featureSubsets) > 0) {
        percentages <- lapply(featureSubsets, function(x) {
            rowIdx <- rownames(object) %in% x
            if (sum(rowIdx) == 0) {
                return(rep(0, ncol(object)))
            } else {
                return(colSums(rawData(object)[rowIdx, , drop = FALSE]) /
                           colSums(rawData(object)) * 100)
            }
        })
        results <- cbind(results, as.data.frame(percentages))
    }
    list(cell = results, feature = nCell)
}

#' Calculate proportion mitochondrial contribution
#' @description
#' Calculates proportion of mitochondrial contribution based on raw or
#' normalized data.
#' @param object \code{liger} object.
#' @param use.norm \bold{Deprecated} Whether to use cell normalized data in
#' calculating contribution. Default \code{FALSE}.
#' @param pattern Regex pattern for identifying mitochondrial genes. Default
#' \code{"^mt-"} for mouse.
#' @return Named vector containing proportion of mitochondrial contribution for
#' each cell.
#' @export
#' @note
#' \code{getProportionMito} will be deprecated because
#' \code{\link{runGeneralQC}} generally covers and expands its use case.
#' @examples
#' # Example dataset does not contain MT genes, expected to see a message
#' pbmc$mito <- getProportionMito(pbmc)
getProportionMito <- function(object, use.norm = FALSE, pattern = "^mt-") {
    lifecycle::deprecate_warn("1.99.0", "getProportionMito()",
                              "runGeneralQC()")
    result <- numeric()
    for (d in names(object)) {
        ld <- dataset(object, d)
        mitoGeneIdx <- grep(pattern, rownames(ld))
        if (isTRUE(use.norm)) {
            pctMT <- colSums(normData(ld)[mitoGeneIdx, , drop = FALSE]) /
                colSums(normData(ld))
        } else {
            pctMT <- colSums(rawData(ld)[mitoGeneIdx, , drop = FALSE]) /
                colSums(rawData(ld))
        }
        names(pctMT) <- colnames(ld)
        result <- c(result, pctMT)
    }
    if (all(result == 0)) {
        cli::cli_alert_warning("Zero proportion detected in all cells")
    }
    return(result)
}

#' Remove missing cells or features from liger object
#' @param object \linkS4class{liger} object
#' @param orient Choose to remove non-expressing features (\code{"feature"}),
#' empty barcodes (\code{"cell"}), or both of them (\code{"both"}). Default
#' \code{"both"}.
#' @param minCells Keep features that are expressed in at least this number of
#' cells, calculated on a per-dataset base. A single value for all datasets or
#' a vector for each dataset. Default \code{NULL} only removes none expressing
#' features.
#' @param minFeatures Keep cells that express at least this number of features,
#' calculated on a per-dataset base. A single value for all datasets or a vector
#' for each dataset. Default \code{NULL} only removes none expressing cells.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to be processed. Default
#' \code{NULL} removes empty entries from all datasets.
#' @param newH5 Logical, whether to create a new H5 file on disk for each
#' H5-based dataset on subset. Default \code{TRUE}
#' @param filenameSuffix When subsetting H5-based datasets to new H5 files, this
#' suffix will be added to all the filenames. Default \code{"removeMissing"}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param ... Arguments passed to \code{\link{subsetLigerDataset}}
#' @return Updated (subset) \code{object}.
#' @export
#' @rdname removeMissing
#' @examples
#' # The example dataset does not contain non-expressing genes or empty barcodes
#' pbmc <- removeMissing(pbmc)
removeMissing <- function(
        object,
        orient = c("both", "feature", "cell"),
        minCells = NULL,
        minFeatures = NULL,
        useDatasets = NULL,
        newH5 = TRUE,
        filenameSuffix = "removeMissing",
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    orient <- match.arg(orient)
    useDatasets <- .checkUseDatasets(object, useDatasets)
    minCells <- minCells %||% rep(0, length(useDatasets))
    minCells <- .checkArgLen(minCells, length(useDatasets), class = "numeric")
    names(minCells) <- useDatasets
    minFeatures <- minFeatures %||% rep(0, length(useDatasets))
    minFeatures <- .checkArgLen(minFeatures, length(useDatasets), class = "numeric")
    names(minFeatures) <- useDatasets
    rmFeature <- ifelse(orient %in% c("both", "feature"), TRUE, FALSE)
    rmCell <- ifelse(orient %in% c("both", "cell"), TRUE, FALSE)
    datasets.new <- list()
    subsetted <- c()
    for (d in useDatasets) {
        ld <- dataset(object, d)
        if (rmFeature) {
            featureIdx <- which(featureMeta(ld)$nCell > minCells[d])
        } else {
            featureIdx <- seq_len(nrow(ld))
        }
        rmFeatureDataset <- length(featureIdx) != nrow(ld)
        if (rmCell) {
            cellIdx <- object$dataset == d & object$nGene > minFeatures[d]
            cellIdx <- colnames(object)[cellIdx]
            cellIdx <- which(colnames(ld) %in% cellIdx)
        } else {
            cellIdx <- seq_len(ncol(ld))
        }
        rmCellDataset <- length(cellIdx) != ncol(ld)
        subsetted <- c(subsetted, any(c(rmFeatureDataset, rmCellDataset)))
        if (any(c(rmFeatureDataset, rmCellDataset))) {
            if (isTRUE(verbose)) {
                cli::cli_alert_info("Removing missing in dataset: {.val {d}}")
            }

            datasets.new[[d]] <- subsetLigerDataset(
                ld,
                featureIdx = featureIdx,
                cellIdx = cellIdx,
                filenameSuffix = filenameSuffix,
                verbose = verbose,
                newH5 = newH5,
                ...
            )
        } else {
            datasets.new[[d]] <- ld
        }
    }
    if (any(subsetted)) {
        allCells <- unlist(lapply(datasets.new, colnames), use.names = FALSE)
        object <- methods::new(
            "liger",
            datasets = datasets.new,
            cellMeta = cellMeta(object, cellIdx = allCells,
                                drop = FALSE),
            varFeatures = character(),
            H.norm = object@H.norm[allCells, , drop = FALSE]
        )
    }
    return(object)
}

#' @rdname removeMissing
#' @export
#' @param slot.use \bold{Deprecated}. Always look at \code{rawData} slot of
#' inner \linkS4class{ligerDataset} objects.
#' @param use.cols \bold{Deprecated}. Previously means "treating each column as
#' a cell" when \code{TRUE}, now means \code{orient="cell"}.
#' @note
#' \code{removeMissingObs} will be deprecated. \code{removeMissing} covers and
#' expands the use case and should be easier to understand.
removeMissingObs <- function(
        object,
        slot.use = NULL,
        use.cols = TRUE,
        verbose = getOption("ligerVerbose", TRUE)) {
    lifecycle::deprecate_warn("1.99.0", "removeMissingObs()",
                              "removeMissing()")
    if (!missing(slot.use)) {
        cli::cli_alert_warning("Argument {.code slot.use} is deprecated and ignored.")
    }
    orient <- ifelse(isTRUE(use.cols), "cell", "gene")
    object <- removeMissing(object, orient, verbose = verbose)
    return(object)
}



################################ Normalize #####################################

#' Normalize raw counts data
#' @description Perform library size normalization on raw counts input. As for
#' the preprocessing step of iNMF integration, by default we don't multiply the
#' normalized values with a scale factor, nor do we take the log transformation.
#' Applicable S3 methods can be found in Usage section.
#'
#' \code{normalizePeak} is designed for datasets of "atac" modality, i.e. stored
#' in \linkS4class{ligerATACDataset}. S3 method for various container object is
#' not supported yet due to difference in architecture design.
#' @param object \linkS4class{liger} object
#' @param ... Arguments to be passed to S3 methods. The "liger" method calls
#' the "ligerDataset" method, which then calls "dgCMatrix" method.
#' \code{normalizePeak} directly calls \code{normalize.dgCMatrix}.
#' @return Updated \code{object}.
#' \itemize{
#'  \item{dgCMatrix method - Returns processed dgCMatrix object}
#'  \item{ligerDataset method - Updates the \code{normData} slot of the object}
#'  \item{liger method - Updates the \code{normData} slot of chosen datasets}
#'  \item{Seurat method - Adds a named layer in chosen assay (V5), or update the
#'  \code{data} slot of the chosen assay (<=V4)}
#'  \item{\code{normalizePeak} - Updates the \code{normPeak} slot of chosen
#'  datasets.}
#' }
#' @rdname normalize
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
normalize <- function(object, ...) {
    UseMethod("normalize", object)
}

#' @rdname normalize
#' @export
#' @method normalize dgCMatrix
#' @param log Logical. Whether to do a \code{log(x + 1)} transform on the
#' normalized data. Default \code{TRUE}.
#' @param scaleFactor Numeric. Scale the normalized expression value by this
#' factor before transformation. \code{NULL} for not scaling. Default
#' \code{1e4}.
normalize.dgCMatrix <- function(
        object,
        log = FALSE,
        scaleFactor = NULL,
        ...
) {
    scaleFactor <- .checkArgLen(scaleFactor, ncol(object), repN = TRUE, class = "numeric")
    if (!is.null(scaleFactor) && any(scaleFactor <= 0)) {
        cli::cli_alert_danger("Invalid {.code scaleFactor} given. Setting to {.code NULL}.")
        scaleFactor <- NULL
    }
    normed <- object
    normed@x <- object@x / rep.int(Matrix::colSums(object), diff(object@p))
    if (!is.null(scaleFactor)) normed <- normed * scaleFactor
    if (isTRUE(log)) normed <- log1p(normed)
    return(normed)
}

#' @rdname normalize
#' @export
#' @param chunk Integer. Number of maximum number of cells in each chunk when
#' working on HDF5 file based ligerDataset. Default \code{1000}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @method normalize ligerDataset
normalize.ligerDataset <- function(
        object,
        chunk = 1000,
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    if (!isH5Liger(object)) {
        normData(object) <- normalize(rawData(object), ...)
    } else {
        # Initialize result
        results <- list(
            geneSumSq = rep(0, nrow(object)),
            geneMeans = rep(0, nrow(object))
        )
        h5file <- getH5File(object)
        resultH5Path <- "normData"
        # Use safe create here in practice
        ## This creates the CSC sparse non-zero value array
        safeH5Create(object = object, dataPath = resultH5Path,
                     dims = rawData(object)$dims, dtype = "double",
                     chunkSize = rawData(object)$chunk_dims)
        # Chunk run
        results <- H5Apply(
            object,
            function(chunk, sparseXIdx, cellIdx, values) {
                normChunk <- normalize(chunk)
                h5file[[resultH5Path]][sparseXIdx] <- normChunk@x
                row_sums <- rowSums(normChunk)
                values$geneSumSq <- values$geneSumSq +
                    rowSums(normChunk * normChunk)
                values$geneMeans <- values$geneMeans + row_sums
                return(values)
            },
            init = results, useData = "rawData",
            chunkSize = chunk, verbose = verbose
        )
        results$geneMeans <- results$geneMeans / ncol(object)
        featureMeta(object, check = FALSE)$geneMeans <- results$geneMeans
        featureMeta(object, check = FALSE)$geneSumSq <- results$geneSumSq
        normData(object, check = FALSE) <- h5file[[resultH5Path]]
        h5fileInfo(object, "normData", check = FALSE) <- resultH5Path
    }
    return(object)
}

#' @rdname normalize
#' @export
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to be normalized. Should specify ATACseq
#' datasets when using \code{normalizePeak}. Default \code{NULL} normalizes all
#' valid datasets.
#' @param format.type,remove.missing \bold{Deprecated}. The functionality of
#' these is covered through other parts of the whole workflow and is no long
#' needed. Will be ignored if specified.
#' @method normalize liger
normalize.liger <- function(
        object,
        useDatasets = NULL,
        verbose = getOption("ligerVerbose", TRUE),
        format.type = NULL,
        remove.missing = NULL,
        ...
) {
    .deprecateArgs(defunct = c("format.type", "remove.missing"))
    .checkObjVersion(object)
    useDatasets <- .checkUseDatasets(object, useDatasets)
    object <- recordCommand(object, ..., dependencies = "hdf5r")
    for (d in useDatasets) {
        if (isTRUE(verbose)) cliID <- cli::cli_process_start("Normalizing datasets {.val {d}}")
        # `d` is the name of each dataset
        ld <- dataset(object, d)
        ld <- normalize(ld, verbose = verbose, ...)
        datasets(object, check = FALSE)[[d]] <- ld
    }
    if (isTRUE(verbose)) cli::cli_process_done(id = cliID)
    object
}

#' @rdname normalize
#' @export
#' @param assay Name of assay to use. Default \code{NULL} uses current active
#' assay.
#' @param save For Seurat>=4.9.9, the name of layer to store normalized data.
#' Default \code{"ligerNormData"}. For older Seurat, stored to \code{data} slot.
#' @param layer Where the input raw counts should be from. Default
#' \code{"counts"}. For older Seurat, always retrieve from \code{counts} slot.
#' @method normalize Seurat
normalize.Seurat <- function(
        object,
        assay = NULL,
        layer = "counts",
        save = "ligerNormData",
        ...
) {
    raw <- .getSeuratData(object, layer = layer, slot = "counts",
                          assay = assay)
    if (!is.list(raw)) normed <- normalize(raw, ...)
    else normed <- lapply(raw, normalize, ...)
    object <- .setSeuratData(object, layer = layer, save = save, slot = "data",
                             value = normed, assay = assay)
    return(object)
}

#' @rdname normalize
#' @export
normalizePeak <- function(
        object,
        useDatasets = NULL,
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    useDatasets <- .checkUseDatasets(object, useDatasets, modal = "atac")
    object <- recordCommand(object, ..., dependencies = "hdf5r")
    for (d in useDatasets) {
        if (isTRUE(verbose)) cliID <- cli::cli_process_start("Normalizing peak of dataset: {.val {d}}")
        # `d` is the name of each dataset
        ld <- dataset(object, d)
        normPeak(ld, check = FALSE) <- normalize(rawPeak(ld), ...)
        datasets(object, check = FALSE)[[d]] <- ld
        if (isTRUE(verbose)) cli::cli_process_done(id = cliID)
    }
    object
}

############################### Select Genes ###################################

#' Select a subset of informative genes
#' @description This function identifies highly variable genes from each dataset
#' and combines these gene sets (either by union or intersection) for use in
#' downstream analysis. Assuming that gene expression approximately follows a
#' Poisson distribution, this function identifies genes with gene expression
#' variance above a given variance threshold (relative to mean gene expression).
#' Alternatively, we allow selecting a desired number of genes for each dataset
#' by ranking the relative variance, and then take the combination.
#' @export
#' @rdname selectGenes
#' @param object A \linkS4class{liger}, \linkS4class{ligerDataset} or
#' \code{Seurat} object, with normalized data available (no scale factor
#' multipled nor log transformed).
#' @param thresh Variance threshold used to identify variable genes. Higher
#' threshold results in fewer selected genes. Liger and Seurat S3 methods accept
#' a single value or a vector with specific threshold for each dataset in
#' \code{useDatasets}.* Default \code{0.1}.
#' @param nGenes Number of genes to find for each dataset. By setting this,
#' we optimize the threshold used for each dataset so that we get \code{nGenes}
#' selected features for each dataset. Accepts single value or a vector for
#' dataset specific setting matching \code{useDataset}.* Default \code{NULL}
#' does not optimize.
#' @param alpha Alpha threshold. Controls upper bound for expected mean gene
#' expression. Lower threshold means higher upper bound. Default \code{0.99}.
#' @param combine How to combine variable genes selected from all datasets.
#' Choose from \code{"union"} or \code{"intersection"}. Default \code{"union"}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param ... Arguments passed to other methods.
#' @return Updated object
#' \itemize{
#'  \item{liger method - Each involved dataset stored in
#'  \linkS4class{ligerDataset} is updated with its \code{\link{featureMeta}}
#'  slot and \code{varUnsharedFeatures} slot (if requested with
#'  \code{useUnsharedDatasets}), while \code{\link{varFeatures}(object)} will be
#'  updated with the final combined gene set.}
#'  \item{Seurat method - Final selection will be updated at
#'  \code{Seurat::VariableFeatures(object)}. Per-dataset information is
#'  stored in the \code{meta.features} slot of the chosen Assay.}
#' }
#' @examples
#' pbmc <- normalize(pbmc)
#' # Select basing on thresholding the relative variance
#' pbmc <- selectGenes(pbmc, thresh = .1)
#' # Select specified number for each dataset
#' pbmc <- selectGenes(pbmc, nGenes = c(60, 60))
selectGenes <- function(
        object,
        thresh = .1,
        nGenes = NULL,
        alpha = .99,
        ...
) {
    UseMethod("selectGenes", object)
}

#' @export
#' @rdname selectGenes
#' @method selectGenes liger
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to use for shared variable feature
#' selection. Default \code{NULL} uses all datasets.
#' @param useUnsharedDatasets A character vector of the names, a numeric or
#' logical vector of the index of the datasets to use for finding unshared
#' variable features. Default \code{NULL} does not attempt to find unshared
#' features.
#' @param unsharedThresh The same thing as \code{thresh} that is applied to test
#' unshared features. A single value for all datasets in
#' \code{useUnsharedDatasets} or a vector for dataset-specific setting.* Default
#' \code{0.1}.
#' @param chunk Integer. Number of maximum number of cells in each chunk, when
#' gene selection is applied to any HDF5 based dataset. Default \code{1000}.
#' @param var.thresh,alpha.thresh,num.genes,datasets.use,unshared.datasets,unshared.thresh \bold{Deprecated}.
#' These arguments are renamed and will be removed in the future. Please see
#' function usage for replacement.
#' @param tol,do.plot,cex.use,unshared \bold{Deprecated}. Gene variability
#' metric is now visualized with separated function
#' \code{\link{plotVarFeatures}}. Users can now set none-NULL
#' \code{useUnsharedDatasets} to select unshared genes, instead of having to
#' switch \code{unshared} on.
selectGenes.liger <- function(
        object,
        thresh = .1,
        nGenes = NULL,
        alpha = .99,
        useDatasets = NULL,
        useUnsharedDatasets = NULL,
        unsharedThresh = .1,
        combine = c("union", "intersection"),
        chunk = 1000,
        verbose = getOption("ligerVerbose", TRUE),
        var.thresh = thresh,
        alpha.thresh = alpha,
        num.genes = nGenes,
        datasets.use = useDatasets,
        unshared.datasets = useUnsharedDatasets,
        unshared.thresh = unsharedThresh,
        tol = NULL,
        do.plot = NULL,
        cex.use = NULL,
        unshared = NULL,
        ...
) {
    combine <- match.arg(combine)
    .deprecateArgs(replace = list(var.thresh = "thresh",
                                  alpha.thresh = "alpha",
                                  num.genes = "nGenes",
                                  datasets.use = "useDatasets",
                                  unshared.datasets = "useUnsharedDatasets",
                                  unshared.thresh = "unsharedThresh"),
                   defunct = c("tol", "do.plot", "cex.use"))
    object <- recordCommand(object, ...)
    datasetShared <- .checkUseDatasets(object, useDatasets)
    if (!is.null(useUnsharedDatasets))
        datasetUnshared <- .checkUseDatasets(object, useUnsharedDatasets)
    else datasetUnshared <- NULL
    useDatasets <- union(datasetShared, datasetUnshared)
    thresh <- .checkArgLen(thresh, length(datasetShared), class = "numeric")
    nGenes <- .checkArgLen(nGenes, length(datasetShared), class = "numeric")
    unsharedThresh <- .checkArgLen(unsharedThresh, length(datasetUnshared), class = "numeric")
    sharedFeature <- Reduce(intersect, lapply(datasets(object), rownames))
    selectList <- list()
    for (d in useDatasets) {
        if (isTRUE(verbose))
            cli::cli_alert_info("Selecting variable features for dataset {.val {d}}")
        ld <- dataset(object, d)
        thresh_i <- thresh[datasetShared == d]
        nGenes_i <- nGenes[datasetShared == d]
        unsharedThresh_i <- unsharedThresh[datasetUnshared == d]
        ld <- .selectGenes.ligerDataset(
            ld, sharedFeature = sharedFeature, thresh = thresh_i,
            nGenes = nGenes_i, unshared = d %in% datasetUnshared,
            unsharedThresh = unsharedThresh_i,
            nUMI = cellMeta(object, "nUMI", useDatasets = d),
            alpha = alpha, chunk = chunk, verbose = verbose
        )
        selectList[[d]] <- rownames(ld)[featureMeta(ld)$isVariable]
        datasets(object, check = FALSE)[[d]] <- ld
    }
    if (combine == "union") selected <- Reduce(union, selectList)
    else selected <- Reduce(intersect, selectList)
    if (length(selected) == 0) {
        cli::cli_alert_danger("No genes were selected. Lower {.code thresh} values or set {.code combine = 'union'}")
    } else {
        if (isTRUE(verbose))
            cli::cli_alert_success("Finally {length(selected)} shared variable feature{?s} are selected.")
    }
    varFeatures(object) <- selected
    for (d in names(object)) {
        ld <- dataset(object, d)
        featureMeta(ld, check = FALSE)$selected <- rownames(ld) %in% selected
        datasets(object, check = FALSE)[[d]] <- ld
    }
    return(object)
}

#' @param sharedFeature Character vector, the feature names that are common to
#' all datasets involved for the selection. Mostly set internally by the liger
#' S3 method. Default \code{NULL} tests with all available features.
#' @param unshared Logical, whether to select variable ones from unshared
#' features.
#' @param nUMI A vector of prior QC information, number of total counts per cell
#' derived with the raw count matrix of a dataset. Mostly set internally by the
#' liger S3 method.
#' @note
#' For \code{thresh}, \code{unsharedThresh} and \code{nGenes}, ligerDataset S3
#' method only accept a single value, which most of the time got passed from
#' upstream method chain.
#' @noRd
.selectGenes.ligerDataset <- function(
        object,
        nUMI,
        sharedFeature = NULL,
        thresh = .1,
        nGenes = NULL,
        unshared = FALSE,
        unsharedThresh = .1,
        alpha = .99,
        chunk = 1000,
        verbose = getOption("ligerVerbose", TRUE)
) {
    if (is.null(normData(object))) cli::cli_abort("Normalized data not available.")
    if (is.null(sharedFeature)) sharedFeature <- rownames(object)
    sharedFeature <- rownames(object) %in% sharedFeature
    unsharedFeature <- !sharedFeature
    selected.shared <- logical(sum(sharedFeature))
    selected.unshared <- logical(sum(unsharedFeature))
    if (isH5Liger(object)) {
        object <- calcGeneVars.H5(object, chunkSize = chunk, verbose = verbose)
    } else {
        featureMeta(object, check = FALSE)$geneMeans <-
            Matrix::rowMeans(normData(object))
        featureMeta(object, check = FALSE)$geneVars <-
            rowVars_sparse_rcpp(normData(object), featureMeta(object)$geneMeans)
    }
    selected.shared <- .selectGenes.withMetric(
        genes = rownames(object)[sharedFeature],
        means = featureMeta(object)$geneMeans[sharedFeature],
        vars = featureMeta(object)$geneVars[sharedFeature],
        nUMI = nUMI, dims = dim(object), thresh = thresh, alpha = alpha,
        nGenes = nGenes
    )
    featureMeta(object, check = FALSE)$isVariable <-
        rownames(object) %in% selected.shared
    if (isTRUE(verbose)) {
        cli::cli_alert_success("... {length(selected.shared)} feature{?s} selected out of {sum(sharedFeature)} shared feature{?s}.")
    }
    if (isTRUE(unshared) && length(unsharedFeature) > 0) {
        selected.unshared <- .selectGenes.withMetric(
            genes = rownames(object)[unsharedFeature],
            means = featureMeta(object)$geneMeans[unsharedFeature],
            vars = featureMeta(object)$geneVars[unsharedFeature],
            nUMI = nUMI, dims = dim(object), thresh = unsharedThresh,
            alpha = alpha, nGenes = nGenes
        )
        object@varUnsharedFeatures <- selected.unshared
        if (isTRUE(verbose)) {
            cli::cli_alert_success("... {length(selected.unshared)} feature{?s} selected out of {sum(unsharedFeature)} unshared feature{?s}.")
        }
    }
    return(object)
}

#' Calculate Gene Variance for ligerDataset object
#' @param object ligerDataset object
#' @param chunkSize Integer for the maximum number of cells in each chunk.
#' Default \code{1000}.
#' @param verbose Logical. Whether to show a progress bar. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @return The input \code{object} with calculated var updated in the H5 file.
#' @noRd
calcGeneVars.H5 <- function(object, chunkSize = 1000,
                            verbose = getOption("ligerVerbose", TRUE)) {
    geneVars <- rep(0, nrow(object))
    geneMeans <- featureMeta(object)$geneMeans
    geneVars <- H5Apply(
        object,
        function(chunk, sparseXIdx, cellIdx, values) {
            values + sumSquaredDeviations(chunk, geneMeans)
        },
        init = geneVars,
        useData = "normData",
        chunkSize = chunkSize,
        verbose = verbose
    )
    geneVars <- geneVars / (ncol(object) - 1)
    featureMeta(object, check = FALSE)$geneVars <- geneVars
    object
}

#' @export
#' @rdname selectGenes
#' @method selectGenes Seurat
#' @param layer Where the input normalized counts should be from. Default
#' \code{"ligerNormData"}. For older Seurat, always retrieve from \code{data}
#' slot.
#' @param assay Name of assay to use. Default \code{NULL} uses current active
#' assay.
#' @param datasetVar Metadata variable name that stores the dataset source
#' annotation. Default \code{"orig.ident"}.
selectGenes.Seurat <- function(
        object,
        thresh = .1,
        nGenes = NULL,
        alpha = .99,
        useDatasets = NULL,
        layer = "ligerNormData",
        assay = NULL,
        datasetVar = "orig.ident",
        combine = c("union", "intersection"),
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    combine <- match.arg(combine)
    assay <- assay %||% SeuratObject::DefaultAssay(object)
    matList <- .getSeuratData(object, layer = layer, slot = "data",
                              assay = assay)
    if (is.list(matList)) {
        # object contain split layers
        names(matList) <- gsub(paste0(layer, "."), "", names(matList))
        datasetVar <- factor(rep(names(matList), sapply(matList, ncol)),
                             levels = names(matList))
    } else {
        datasetVar <- object[[datasetVar, drop = TRUE]]
        if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
        datasetVar <- droplevels(datasetVar)
        matList <- splitRmMiss(matList, datasetVar)
    }
    useDatasets <- useDatasets %||% levels(datasetVar)
    matList <- matList[unique(useDatasets)]

    featureList <- lapply(matList, rownames)
    allshared <- Reduce(intersect, featureList)
    allFeatures <- SeuratObject::Features(object, assay = assay)
    thresh <- .checkArgLen(thresh, nlevels(datasetVar), class = "numeric")

    # Get nUMI metric into list
    nUMIVar <- paste0("nCount_", assay)
    nUMIAll <- object[[nUMIVar]][,1]
    nUMIList <- split(nUMIAll, datasetVar)

    selectList <- list()
    hvg.info <- data.frame(row.names = allFeatures)
    for (d in levels(datasetVar)) {
        if (isTRUE(verbose))
            cli::cli_alert_info("Selecting variable features for dataset: {.val {d}}")
        submat <- matList[[d]]
        # submat <- mat[, datasetVar == d, drop = FALSE]
        # submat <- submat[sort(expressed), , drop = FALSE]
        sharedFeature <- rownames(submat) %in% allshared
        means <- Matrix::rowMeans(submat)
        hvg.info[[paste0("liger.mean.", d)]] <- 0
        hvg.info[rownames(submat), paste0("liger.mean.", d)] <- means
        vars <- rowVars_sparse_rcpp(submat, means)
        hvg.info[[paste0("liger.variance.", d)]] <- 0
        hvg.info[rownames(submat), paste0("liger.variance.", d)] <- vars
        thresh_i <- thresh[levels(datasetVar) == d]
        selected <- .selectGenes.withMetric(
            genes = rownames(submat)[sharedFeature],
            means = means[sharedFeature],
            vars = vars[sharedFeature],
            nUMI = nUMIList[[d]],
            dims = dim(submat),
            thresh = thresh_i, alpha = alpha,
            nGenes = nGenes
        )
        if (isTRUE(verbose)) {
            cli::cli_alert_success("... {length(selected)} features selected out of {length(allshared)} shared features")
        }
        selectList[[d]] <- selected
    }
    if (combine == "union") selected <- Reduce(union, selectList)
    else selected <- Reduce(intersect, selectList)
    if (length(selected) == 0) {
        cli::cli_alert_danger("No genes were selected. Lower {.code thresh} values or set {.code combine = 'union'}")
    } else {
        if (isTRUE(verbose))
            cli::cli_alert_success("Finally {length(selected)} shared variable features selected.")
    }
    hvg.info$liger.variable <- allFeatures %in% selected
    assayObj <- Seurat::GetAssay(object, assay = assay)
    assayObj[[names(hvg.info)]] <- hvg.info
    object[[assay]] <- assayObj
    SeuratObject::VariableFeatures(object, assay = assay) <- selected
    return(object)
}

# returns selected gene names
.selectGenes.withMetric <- function(
        genes,
        means,
        vars,
        nUMI,
        dims,
        thresh = .1,
        alpha = .99,
        nGenes = NULL
) {
    nolan_constant <- mean((1 / nUMI))
    alphathresh.corrected <- alpha / dims[1]
    geneMeanUpper <- means +
        stats::qnorm(1 - alphathresh.corrected / 2) *
        sqrt(means * nolan_constant / dims[2])
    basegenelower <- log10(means * nolan_constant)
    if (!is.null(nGenes)) {
        preselect <- vars / nolan_constant > geneMeanUpper &
            log10(vars) > basegenelower
        relativeVar <- log10(vars) - basegenelower
        names(relativeVar) <- genes
        relativeVar <- relativeVar[preselect]
        relativeVar <- relativeVar[order(relativeVar, decreasing = TRUE)]
        selected <- names(relativeVar[seq(nGenes)])
    } else {
        selected <- genes[vars / nolan_constant > geneMeanUpper &
                              log10(vars) > basegenelower + thresh]
    }
    return(selected)
}

#' Plot the variance vs mean of feature expression
#' @description For each dataset where the feature variablitity is calculated,
#' a plot of log10 feature expression variance and log10 mean will be produced.
#' Features that are considered as variable would be highlighted in red.
#' @param object \linkS4class{liger} object. \code{\link{selectGenes}} needs to
#' be run in advance.
#' @param combinePlot Logical. If \code{TRUE}, sub-figures for all datasets will
#' be combined into one plot. if \code{FALSE}, a list of plots will be returned.
#' Default \code{TRUE}.
#' @param dotSize Controls the size of dots in the main plot. Default
#' \code{0.8}.
#' @param ... More theme setting parameters passed to
#' \code{\link{.ggplotLigerTheme}}.
#' @return \code{ggplot} object when \code{combinePlot = TRUE}, a list of
#' \code{ggplot} objects when \code{combinePlot = FALSE}
#' @export
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' plotVarFeatures(pbmc)
plotVarFeatures <- function(
        object,
        combinePlot = TRUE,
        dotSize = 1,
        ...
) {
    plotList <- list()
    maxVar <- max(sapply(datasets(object),
                         function(ld) max(log10(featureMeta(ld)$geneVars))))
    minVar <- min(sapply(datasets(object),
                         function(ld) min(log10(featureMeta(ld)$geneVars))))
    maxMean <- max(sapply(datasets(object),
                          function(ld) max(log10(featureMeta(ld)$geneMeans))))
    minMean <- min(sapply(datasets(object),
                          function(ld) min(log10(featureMeta(ld)$geneMeans))))
    for (d in names(object)) {
        ld <- dataset(object, d)
        trx_per_cell <- cellMeta(object, "nUMI", cellIdx = object$dataset == d)
        nolan_constant <- mean((1 / trx_per_cell))

        data <- .DataFrame.as.data.frame(featureMeta(ld))
        nSelect <- sum(data$isVariable)
        data$geneMeans <- log10(data$geneMeans)
        data$geneVars <- log10(data$geneVars)
        data$isVariable <- factor(data$isVariable,
                                  levels = c("TRUE", "FALSE"))
        p <- ggplot2::ggplot(
            data,
            ggplot2::aes(x = .data[["geneMeans"]],
                         y = .data[["geneVars"]],
                         color = .data[["isVariable"]])
        ) +
            ggplot2::geom_point(size = dotSize, stroke = 0) +
            ggplot2::geom_abline(intercept = log10(nolan_constant), slope = 1,
                                 color = "purple") +
            ggplot2::xlim(minMean, maxMean) +
            ggplot2::ylim(minVar, maxVar)
        p <- .ggplotLigerTheme(p, title = d,
                               subtitle = paste0(nSelect, " variable features"),
                               xlab = "Gene Expression Mean (log10)",
                               ylab = "Gene Expression Variance (log10)",
                               legendColorTitle = "Variable\nfeature",
                               colorLabels = c("TRUE", "FALSE"),
                               colorValues = c("RED", "BLACK"),
                               ...)
        plotList[[d]] <- p
    }
    if (isTRUE(combinePlot)) {
        suppressWarnings({
            legend <- cowplot::get_legend(plotList[[1]])
        })
        plotList <- lapply(plotList, function(x) {
            x + ggplot2::theme(legend.position = "none")
        })
        combined <- cowplot::plot_grid(plotlist = plotList,
                                       align = "hv",
                                       axis = "tblr")
        combined <- cowplot::plot_grid(combined, legend, rel_widths = c(5,1))
        return(combined)
    }
    else return(plotList)
}


#' Select variable genes from one dataset with Seurat VST method
#' @description
#' Seurat FindVariableFeatures VST method. This allows the selection of a fixed
#' number of variable features, but only applies to one dataset. No
#' normalization is needed in advance.
#' @param object A \linkS4class{liger} object.
#' @param useDataset The names, a numeric or logical index of the dataset to
#' be considered for selection.
#' @param n Number of variable features needed. Default \code{2000}.
#' @param loessSpan Loess span parameter used when fitting the variance-mean
#' relationship. Default \code{0.3}.
#' @param clipMax After standardization values larger than \code{clipMax} will
#' be set to \code{clipMax}. Default \code{"auto"} sets this value to the square
#' root of the number of cells.
#' @param useShared Logical. Whether to only select from genes shared by all
#' dataset. Default \code{TRUE}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @references Seurat::FindVariableFeatures.default(selection.method = "vst")
#' @export
#' @examples
#' pbmc <- selectGenesVST(pbmc, "ctrl", n = 50)
selectGenesVST <- function(
        object,
        useDataset,
        n = 2000,
        loessSpan = 0.3,
        clipMax = "auto",
        useShared = TRUE,
        verbose = getOption("ligerVerbose", TRUE))
{
    useDataset <- .checkUseDatasets(object, useDataset)
    useDataset <- .checkArgLen(useDataset, 1)
    if (isTRUE(verbose)) {
        cli::cli_alert_info("Selecting top {n} HVG{?s} with VST method for dataset: {.val {useDataset}}")
    }
    ld <- dataset(object, useDataset)
    data <- rawData(ld)
    if (isTRUE(useShared)) {
        useGenes <- lapply(datasets(object), rownames)
        useGenes <- Reduce(intersect, useGenes)
    } else {
        useGenes <- rownames(data)
    }
    if (isTRUE(verbose)) {
        cli::cli_alert_info("... Totally {length(useGenes)} {ifelse(useShared, 'shared', 'dataset specific')} genes to be selected from.")
    }
    if (clipMax == "auto") {
        clipMax <- sqrt(ncol(data))
    }
    hvf.info <- data.frame(mean = Matrix::rowMeans(data))
    hvf.info$variance <- rowVars_sparse_rcpp(data, hvf.info$mean)
    not.const <- hvf.info$variance > 0
    hvf.info$variance.expected <- 0
    fit <- stats::loess(formula = log10(variance) ~ log10(mean),
                        data = hvf.info[not.const, ],
                        span = loessSpan)
    hvf.info$variance.expected[not.const] <- 10^fit$fitted

    # `SparseRowVarStd` Rcpp function reimplemented with Armadillo,
    # Seurat original uses Eigen
    hvf.info$variance.standardized <- SparseRowVarStd(
        x = data,
        mu = hvf.info$mean,
        sd = sqrt(hvf.info$variance.expected),
        vmax = clipMax
    )
    colnames(hvf.info) <- c("mean", "var", "vst.variance.expected",
                            "vst.variance.standardized")
    rank <- order(hvf.info$vst.variance.standardized, decreasing = TRUE)
    selected <- rownames(ld)[rank][seq(n)]
    selected <- selected[selected %in% useGenes]
    varFeatures(object) <- selected
    fm <- featureMeta(ld)
    for (col in colnames(hvf.info)) fm[[col]] <- hvf.info[[col]]
    fm$selected <- FALSE
    fm$selected[rank[seq(n)]] <- TRUE
    featureMeta(ld, check = FALSE) <- fm
    datasets(object, check = FALSE)[[useDataset]] <- ld
    return(object)
}





############################# Scale Not Center #################################

#' Scale genes by root-mean-square across cells
#' @description This function scales normalized gene expression data after
#' variable genes have been selected. We do not mean-center the data before
#' scaling in order to address the non-negativity constraint of NMF.
#' Computation applied to each normalized dataset matrix can form the following
#' equation:
#'
#' \deqn{S_{i,j}=\frac{N_{i,j}}{\sqrt{\sum_{p}^{n}\frac{N_{i,p}^2}{n-1}}}}
#'
#' Where \eqn{N} denotes the normalized matrix for an individual dataset,
#' \eqn{S} is the output scaled matrix for this dataset, and \eqn{n} is the
#' number of cells in this dataset. \eqn{i, j} denotes the specific gene and
#' cell index, and \eqn{p} is the cell iterator.
#'
#' Please see detailed section below for explanation on methylation dataset.
#' @note
#' Since the scaling on genes is applied on a per dataset base, other scaling
#' methods that apply to a whole concatenated matrix of multiple datasets might
#' not be considered as equivalent alternatives, even if options like
#' \code{center} are set to \code{FALSE}. Hence we implemented an efficient
#' solution that works under such circumstance, provided with the Seurat S3
#' method.
#' @section Methylation dataset:
#' Because gene body mCH proportions are negatively correlated with gene
#' expression level in neurons, we need to reverse the direction of the
#' methylation data before performing the integration. We do this by simply
#' subtracting all values from the maximum methylation value. The resulting
#' values are positively correlated with gene expression. This will only be
#' applied to variable genes detected in prior. Please make sure that argument
#' \code{modal} is set accordingly when running \code{\link{createLiger}}. In
#' this way, this function can automatically detect it and take proper action.
#' If it is not set, users can still manually have the equivalent processing
#' done by doing \code{scaleNotCenter(lig, useDataset = c("other", "datasets"))},
#' and then \code{\link{reverseMethData}(lig, useDataset = c("meth", "datasets"))}.
#' @param object \linkS4class{liger} object, \linkS4class{ligerDataset} object,
#' \linkS4class{dgCMatrix}, or a Seurat object.
#' @param ... Arguments passed to other methods. The order goes by: "liger"
#' method calls "ligerDataset" method", which then calls "dgCMatrix" method.
#' "Seurat" method directly calls "dgCMatrix" method.
#' @return Updated \code{object}
#' \itemize{
#'  \item{dgCMatrix method - Returns scaled dgCMatrix object}
#'  \item{ligerDataset method - Updates the \code{scaleData} and
#'        \code{scaledUnsharedData} (if unshared variable feature available) slot
#'        of the object}
#'  \item{liger method - Updates the \code{scaleData} and
#'        \code{scaledUnsharedData} (if unshared variable feature available) slot
#'        of chosen datasets}
#'  \item{Seurat method - Adds a named layer in chosen assay (V5), or update the
#'  \code{scale.data} slot of the chosen assay (<=V4)}
#' }
#' @export
#' @rdname scaleNotCenter
#' @examples
#' pbmc <- normalize(pbmc)
#' pbmc <- selectGenes(pbmc)
#' pbmc <- scaleNotCenter(pbmc)
scaleNotCenter <- function(object, ...) {
    UseMethod("scaleNotCenter", object)
}

#' @export
#' @rdname scaleNotCenter
#' @method scaleNotCenter dgCMatrix
scaleNotCenter.dgCMatrix <- function(
        object,
        ...)
{
    if (nrow(object) == 0) return(object)
    scaled <- scaleNotCenter_byRow_rcpp(object)
    scaled@x[is.na(scaled@x)] <- 0 # Is this really happening?
    dimnames(scaled) <- dimnames(object)
    return(scaled)
}

#' @export
#' @rdname scaleNotCenter
#' @method scaleNotCenter ligerDataset
#' @param features Character, numeric or logical index that choose the variable
#' feature to be scaled. "liger" method by default uses
#' \code{\link{varFeatures}(object)}. "ligerDataset" method by default uses all
#' features. "Seurat" method by default uses
#' \code{Seurat::VariableFeatures(object)}.
#' @param chunk Integer. Number of maximum number of cells in each chunk, when
#' scaling is applied to any HDF5 based dataset. Default \code{1000}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
scaleNotCenter.ligerDataset <- function(
        object,
        features = NULL,
        chunk = 1000,
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    features <- .idxCheck(object, features, "feature")
    unsharedIdx <- .idxCheck(object, object@varUnsharedFeatures, "feature")
    if (!isH5Liger(object)) {
        scaleData(object) <- scaleNotCenter(
            normData(object)[features, , drop = FALSE]
        )
        if (length(unsharedIdx) > 0)
            scaleUnsharedData(object) <- scaleNotCenter(
                normData(object)[unsharedIdx, , drop = FALSE]
            )
    } else {
        object <- .scaleH5SpMatrix(object, features,
                                   resultH5Path = "scaleDataSparse",
                                   chunk = chunk, verbose = verbose)
        if (length(unsharedIdx) > 0)
            object <- .scaleH5SpMatrix(object, unsharedIdx,
                                       resultH5Path = "scaleUnsharedDataSparse",
                                       chunk = chunk, verbose = verbose)
    }
    return(object)
}

#' @export
#' @rdname scaleNotCenter
#' @method scaleNotCenter ligerMethDataset
scaleNotCenter.ligerMethDataset <- function(
        object,
        features = NULL,
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    raw <- rawData(object)
    scaled <- max(raw) - as.matrix(raw[features, , drop = FALSE])
    # Need more experiment to see if it's worth doing so
    scaled <- methods::as(scaled, "CsparseMatrix")
    scaleData(object, check = FALSE) <- scaled
    return(object)
}

#' @export
#' @rdname scaleNotCenter
#' @method scaleNotCenter liger
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to be scaled but not centered. Default
#' \code{NULL} applies to all datasets.
#' @param remove.missing \bold{Deprecated}. The functionality of this is covered
#' through other parts of the whole workflow and is no long needed. Will be
#' ignored if specified.
scaleNotCenter.liger <- function(
        object,
        useDatasets = NULL,
        features = varFeatures(object),
        verbose = getOption("ligerVerbose", TRUE),
        remove.missing = NULL,
        ...
) {
    .deprecateArgs(defunct = "remove.missing")
    .checkObjVersion(object)
    if (is.null(features) || length(features) == 0) {
        cli::cli_abort("No variable feature specified. Run {.fn selectGenes} first.")
    }
    useDatasets <- .checkUseDatasets(object, useDatasets)
    object <- recordCommand(object, ...,
                            dependencies = c("RcppArmadillo", "Rcpp"))

    for (d in useDatasets) {
        if (isTRUE(verbose)) cliID <- cli::cli_process_start("Scaling dataset {.val {d}}")
        ld <- dataset(object, d)
        ld <- scaleNotCenter(ld, features = features, verbose = verbose, ...)
        datasets(object, check = FALSE)[[d]] <- ld
        if (isTRUE(verbose)) cli::cli_process_done(id = cliID)
    }
    return(object)
}

#' @rdname scaleNotCenter
#' @export
#' @param assay Name of assay to use. Default \code{NULL} uses current active
#' assay.
#' @param save For Seurat>=4.9.9, the name of layer to store normalized data.
#' Default \code{"ligerScaleData"}. For older Seurat, stored to
#' \code{scale.data} slot.
#' @param layer For Seurat>=4.9.9, the name of layer to retrieve normalized
#' data. Default \code{"ligerNormData"}. For older Seurat, always retrieve from
#' \code{data} slot.
#' @param datasetVar Metadata variable name that stores the dataset source
#' annotation. Default \code{"orig.ident"}.
#' @method scaleNotCenter Seurat
scaleNotCenter.Seurat <- function(
        object,
        assay = NULL,
        layer = "ligerNormData",
        save = "ligerScaleData",
        datasetVar = "orig.ident",
        features = NULL,
        ...
) {
    normed <- .getSeuratData(object, layer = layer, slot = "data",
                             assay = assay)

    features <- features %||% SeuratObject::VariableFeatures(object)
    if (!length(features)) {
        cli::cli_abort("No variable feature specified. Run {.fn selectGenes} first")
    }

    if (is.list(normed)) {
        scaled <- lapply(normed, function(x) {
            scaleNotCenter(x[features, , drop = FALSE])
        })
    } else {
        # Condition for all batches in one matrix
        normed <- normed[features, , drop = FALSE]
        # the last [,1] converts data.frame to the vector/factor
        datasetVar <- object[[datasetVar, drop = TRUE]]
        if (!is.factor(datasetVar)) datasetVar <- factor(datasetVar)
        datasetVar <- droplevels(datasetVar)
        nlevel <- nlevels(datasetVar)
        datasetVar <- as.integer(datasetVar) - 1
        # efficient solution where we don't need to split and merge
        scaled <- scaleNotCenter_byRow_perDataset_rcpp(normed, datasetVar, nlevel)
        scaled@x[is.na(scaled@x)] <- 0
        dimnames(scaled) <- list(features, colnames(normed))
    }

    object <- .setSeuratData(object, layer = layer, save = save,
                             slot = "scale.data",
                             value = scaled, assay = assay,
                             denseIfNeeded = TRUE)
    return(object)
}

.scaleH5SpMatrix <- function(ld, featureIdx, resultH5Path, chunk, verbose) {
    features <- rownames(ld)[featureIdx]
    geneSumSq <- featureMeta(ld)$geneSumSq[featureIdx]
    nCells <- ncol(ld)
    geneRootMeanSumSq = sqrt(geneSumSq / (nCells - 1))
    h5file <- getH5File(ld)
    # Count the subset nnz first before creating data space
    nnz <- 0
    nnz <- H5Apply(
        ld, useData = "normData", chunkSize = chunk, verbose = FALSE,
        FUN = function(chunk, sparseXIdx, cellIdx, values) {
            chunk <- chunk[featureIdx, , drop = FALSE]
            values <- values + length(chunk@x)
        },
        init = nnz
    )
    # Create datasets
    dataPath <- paste0(resultH5Path, "/data")
    rowindPath <- paste0(resultH5Path, "/indices")
    colptrPath <- paste0(resultH5Path, "/indptr")
    safeH5Create(ld, dataPath = dataPath, dims = nnz,
                 dtype = "double", chunkSize = 2048)
    safeH5Create(ld, dataPath = rowindPath, dims = nnz,
                 dtype = "int", chunkSize = 2048)
    safeH5Create(ld, dataPath = colptrPath, dims = nCells + 1,
                 dtype = "int", chunkSize = 1024)
    # Process chunks of sparse normData, and write to sparse scaleData
    h5file[[colptrPath]][1] <- 0
    H5Apply(
        ld,
        useData = "normData",
        init = c(1, 0), # [1] record of nnz idx start [2] record of last colptr
        chunkSize = chunk,
        verbose = verbose,
        FUN = function(chunk, sparseXIdx, cellIdx, values) {
            # Subset variable features
            chunk <- chunk[featureIdx, , drop = FALSE]
            # Calculate scale not center
            chunk <- rowDivide_rcpp(chunk, geneRootMeanSumSq)
            chunk@x[is.na(chunk@x)] = 0
            # Write
            nnzRange <- seq(from = values[1], length.out = length(chunk@i))
            h5file[[rowindPath]][nnzRange] <- chunk@i
            h5file[[dataPath]][nnzRange] <- chunk@x
            values[1] <- values[1] + length(nnzRange)
            increColptr <- chunk@p + values[2]
            h5file[[colptrPath]][cellIdx + 1] =
                increColptr[seq(2, length(increColptr))]
            values[2] <- increColptr[length(increColptr)]
            return(values)
        }
    )
    safeH5Create(ld, dataPath = paste0(resultH5Path, "/featureIdx"),
                 dims = length(features), dtype = "int")
    h5file[[paste0(resultH5Path, "/featureIdx")]][1:length(featureIdx)] <-
        featureIdx
    h5fileInfo(ld, "scaleData", check = FALSE) <- resultH5Path
    return(ld)
}

#' Create "scaled data" for DNA methylation datasets
#' @description
#' Because gene body mCH proportions are negatively correlated with gene
#' expression level in neurons, we need to reverse the direction of the
#' methylation data. We do this by simply subtracting all values from the
#' maximum methylation value. The resulting values are positively correlated
#' with gene expression. This will only be applied to variable genes detected in
#' prior.
#' @param object A \linkS4class{liger} object, with variable genes identified.
#' @param useDatasets Required. A character vector of the names, a numeric or
#' logical vector of the index of the datasets that should be identified as
#' methylation data where the reversed data will be created.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @return The input \linkS4class{liger} object, where the \code{scaleData} slot
#' of the specified datasets will be updated with value as described in
#' Description.
#' @export
#' @examples
#' # Assuming the second dataset in example data "pbmc" is methylation data
#' pbmc <- normalize(pbmc, useDatasets = 1)
#' pbmc <- selectGenes(pbmc, datasets.use = 1)
#' pbmc <- scaleNotCenter(pbmc, useDatasets = 1)
#' pbmc <- reverseMethData(pbmc, useDatasets = 2)
reverseMethData <- function(object, useDatasets,
                            verbose = getOption("ligerVerbose", TRUE)) {
    useDatasets <- .checkUseDatasets(object, useDatasets)
    if (is.null(varFeatures(object)) || length(varFeatures(object)) == 0) {
        cli::cli_abort("No variable feature available. Run {.fn selectGenes} first.")
    }
    for (d in useDatasets) {
        ld <- dataset(object, d)
        raw <- rawData(ld)
        if (isTRUE(verbose)) cli::cli_alert_info("Substracting methylation data: {.val {d}}")
        scaleData(ld, check = FALSE) <- methods::as(
            max(raw) - raw[varFeatures(object), , drop = FALSE],
            "CsparseMatrix"
        )
        datasets(object, check = FALSE)[[d]] <- ld
    }
    return(object)
}

