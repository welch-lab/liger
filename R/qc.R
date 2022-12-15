#' General QC for liger object
#' @description Calculate number of UMIs, number of detected features and
#' percentage of feature subset (e.g. mito) expression per cell.
#' @param object liger object with \code{raw.data} available in
#' each ligerDataset embedded
#' @param features Feature names matching the feature subsets that users want to
#' calculate the expression percentage with. A vector for a single subset, or a
#' named list for multiple subset. Default \code{NULL}.
#' @param pattern Regex patterns for matching the feature subsets that users
#' want to calculate the expression percentage with. A vector for a single
#' subset, or a named list for multiple subset. Default \code{NULL}.
#' @param chunkSize Integer number of cells to include in a chunk when working
#' on HDF5 based dataset. Default \code{1000}
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{TRUE}.
#' @return A liger object with \code{nUMI}, \code{nGene} updated
#' in \code{cell.meta(object)}, as well as expression percentage value for each
#' feature subset
#' @export
runGeneralQC <- function(
    object,
    features = NULL,
    pattern = NULL,
    chunkSize = 1000,
    verbose = TRUE
) {
    if (!inherits(object, "liger")) {
        stop("Please use a liger object.")
    }

    # Process the the two arguments all into one named list of feature names
    # before exactly calculate the percentage
    featureSubsets <- list()
    allFeatures <- unique(unlist(lapply(datasets(object), rownames),
                                 use.names = FALSE))
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
            featureSubsets[["featureSubsets_pattern"]] <- pattern
        }
    }

    # Start calculation on each dataset
    newResultNames <- c("nUMI", "nGene", names(featureSubsets))
    for (nrn in newResultNames) object[[nrn]] <- 0
    for (d in names(object)) {
        ld <- dataset(object, d)
        if (isTRUE(verbose))
            message(date(), ' ... calculating QC for dataset "', d, '"')
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
        cell.meta(object)[object$dataset == d, newResultNames] <- results
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
        verbose = TRUE) {
    allFeatures <- rownames(object)
    # Initialize results
    results <- data.frame(row.names = colnames(object))
    results$nUMI <- 0
    results$nGene <- 0
    for (i in names(featureSubsets)) {
        results[[i]] <- 0
    }
    rowIndices <- lapply(featureSubsets, function(x) allFeatures %in% x)

    # Calculate in only one iteration
    results <- H5Apply(
        object,
        init = results,
        useData = "raw.data",
        chunkSize = chunkSize,
        verbose = verbose,
        FUN = function(chunk, sparseXIdx, cellIdx, values) {
            nUMI <- colSums(chunk)
            values$nUMI[cellIdx] <- nUMI
            values$nGene[cellIdx] <- colSums(chunk > 0)
            for (fs in names(rowIndices)) {
                values[[fs]][cellIdx] <-
                    colSums(chunk[rowIndices[[fs]],]) / nUMI * 100
            }
            return(values)
        }
    )
    return(results)
}

#' Calculate general QC on in memory matrix based ligerDataset object
#' @param object ligerDataset object
#' @param featureSubsets  Named list passed from \code{runGeneralQC}
#' @return data.frame
#' @noRd
runGeneralQC.Matrix <- function(
        object,
        featureSubsets = NULL,
        verbose = TRUE) {
    nUMI <- colSums(raw.data(object))
    nGene <- colSums(raw.data(object) > 0)
    results <- data.frame(nUMI = nUMI, nGene = nGene,
                         row.names = colnames(object))
    if (length(featureSubsets) > 0) {
        percentages <- lapply(featureSubsets, function(x) {
            row.idx <- rownames(object) %in% x
            colSums(raw.data(object)[row.idx,]) / colSums(raw.data(object)) * 100
        })
        results <- cbind(results, as.data.frame(percentages))
    }
    results
}
