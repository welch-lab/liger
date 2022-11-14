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
#' @return A liger object with \code{nUMI}, \code{nGene} updated
#' in \code{cell.meta(object)}, as well as expression percentage value for each
#' feature subset
#' @export
runGeneralQC <- function(object, features = NULL, pattern = NULL) {
    if (!inherits(object, "liger")) {
        stop("Please use a liger object.")
    }
    # Quick process for nUMI and nGene first
    nUMI <- unlist(lapply(datasets(object), function(x) {
        # x is ligerDataset object
        colSums(raw.data(x))
    }), use.names = FALSE)
    nGene <- unlist(lapply(datasets(object), function(x) {
        colSums(raw.data(x) > 0)
    }), use.names = FALSE)
    object$nUMI <- nUMI
    object$nGene <- nGene
    # Process the next two arguments all into one named list of feature names
    # before exactly calculate the percentage
    featureSubsets <- list()
    allFeatures <- unique(unlist(lapply(datasets(object), function(x) {
        rownames(raw.data(x))
    }), use.names = FALSE))
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
    # Finally calculate them altogether
    if (!is.null(pattern) | !is.null(features)) {
        percentages <- lapply(featureSubsets, function(x) {
            unlist(lapply(datasets(object), function(d) {
                row.idx <- rownames(raw.data(d)) %in% x
                colSums(raw.data(d)[row.idx,]) / colSums(raw.data(d)) * 100
            }), use.names = FALSE)
        })
        names(percentages) <- paste0(names(percentages), "_pct")
        for (i in seq_along(percentages)) {
            object[[names(percentages)[i]]] <- percentages[[i]]
        }
    }
    return(object)
}
