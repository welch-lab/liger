#' `r lifecycle::badge("experimental")` Batch-aware highly variable gene selection
#' @rdname selectBatchHVG
#' @description
#' Method to select HVGs based on mean dispersions of genes that are highly
#' variable genes in all batches. Using a the top target_genes per batch by
#' average normalize dispersion. If target genes still hasn't been reached,
#' then HVGs in all but one batches are used to fill up. This is continued
#' until HVGs in a single batch are considered.
#'
#' This is an *rliger* implementation of the method originally published in
#' [SCIB](https://scib.readthedocs.io/en/latest/api/scib.preprocessing.hvg_batch.html).
#' We found the potential that it can improve integration under some
#' circumstances, and is currently testing it.
#'
#' This function currently only works for shared features across all datasets.
#' For selection from only part of the datasets and selection for
#' dataset-specific unshared features, please use \code{\link{selectGenes}()}.
#' @references
#' Luecken, M.D., Büttner, M., Chaichoompu, K. et al. (2022), Benchmarking
#' atlas-level data integration in single-cell genomics. *Nat Methods*, 19,
#' 41–50. https://doi.org/10.1038/s41592-021-01336-8.
#' @seealso [selectGenes()]
#' @param object A \code{\linkS4class{liger}} object,
#' \code{\linkS4class{ligerDataset}} object or a sparse/dense matrix. The liger
#' objects must have raw counts available. A direct matrix input is preferably
#' log-1p transformed from CPM normalized counts in cell per column orientation.
#' @param nGenes Integer number of target genes to select. Default \code{2000}.
#' @param features For ligerDataset method, the feature subset to limit the
#' selection to, due to liger downstream requires non-zero features to be
#' considered. Default \code{NULL} selects from all features.
#' @param verbose Logical. Whether to show a progress bar. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param returnStats Logical, for dgCMatrix-method, whether to return a data
#' frame of statistics for all features, or by default \code{FALSE} just return
#' a character vector of selected features.
#' @param ... Arguments passed to S3 methods.
#' @return
#' \itemize{
#' \item{liger-method: Returns the input liger object with the selected genes
#' updated in \code{varFeatures} slot, which can be accessed with
#' \code{varFeatures(object)}. Additionally, the statistics are updated in
#' the \code{featureMeta} slot of each ligerDataset object within the
#' \code{datasets} slot of the \code{object}.}
#' \item{ligerDataset-method: Returns the input ligerDataset object with the
#' statistics updated in the \code{featureMeta} slot.}
#' \item{dgCMatrix-method: By default returns a character vector of selected
#' variable features. If \code{returnStats = TRUE}, returns a data.frame of the
#' statistics.}
#' }
#' @export
#' @examples
#' pbmc <- selectBatchHVG(pbmc, nGenes = 10)
#' varFeatures(pbmc)
selectBatchHVG <- function(
        object,
        ...
) {
    UseMethod("selectBatchHVG", object)
}

#' @export
#' @method selectBatchHVG liger
#' @rdname selectBatchHVG
selectBatchHVG.liger <- function(
        object,
        nGenes = 2000,
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    object <- recordCommand(object, ...)
    sharedFeature <- Reduce(intersect, lapply(datasets(object), rownames))
    hvgDFList <- list()
    for (d in names(object)) {
        if (isTRUE(verbose))
            cli::cli_alert_info("Selecting variable features for dataset {.val {d}}")
        ld <- dataset(object, d)
        ld <- selectBatchHVG(ld, nGenes = nGenes, features = sharedFeature,
                             verbose = verbose, ...)
        object@datasets[[d]] <- ld
        fm <- featureMeta(ld)[sharedFeature, , drop = FALSE]
        fm$dataset <- d
        fm$feature <- rownames(fm)
        hvgDFList[[d]] <- as.data.frame(fm)
    }
    hvgDFAgg <- Reduce(rbind, hvgDFList)
    hvgDFAgg <- hvgDFAgg %>%
        dplyr::mutate(feature = factor(.data[['feature']], levels = sharedFeature)) %>%
        dplyr::group_by(.data[['feature']]) %>%
        dplyr::summarise(
            means = mean(.data[['means']]),
            dispersion = mean(.data[['dispersion']]),
            dispersion_norm = mean(.data[['dispersion_norm']]),
            n_batch = sum(.data[['highly_variable']])
        )
    res <- character()
    enough <- FALSE
    nSelected <- 0
    nNeeded <- nGenes
    n_batch_look <- length(object)
    while (!enough && n_batch_look > 0) {
        hvgDFSub <- hvgDFAgg %>%
            dplyr::filter(.data[['n_batch']] == n_batch_look) %>%
            dplyr::arrange(-.data[['dispersion_norm']])
        if (nrow(hvgDFSub) < nNeeded) {
            if (isTRUE(verbose)) {
                cli::cli_alert_info(
                    "Selected {nrow(hvgDFSub)} gene{?s} that {?is/are} highly variable in {n_batch_look} batch{?es}"
                )
            }
            nSelected <- nSelected + nrow(hvgDFSub)
            nNeeded <- nGenes - nSelected
            n_batch_look <- n_batch_look - 1
        } else {
            if (isTRUE(verbose)) {
                cli::cli_alert_info(
                    "Selected {nNeeded} gene{?s} that {?is/are} highly variable in {n_batch_look} batch{?es}"
                )
            }
            hvgDFSub <- hvgDFSub %>% dplyr::slice_head(n = nNeeded)
            enough <- TRUE
        }
        res <- c(res, as.character(hvgDFSub$feature))
    }
    if (isTRUE(verbose)) {
        if (!enough) {
            cli::cli_alert_warning("Only {length(res)} gene{?s} {?is/are} selected while {nGenes} {?is/are} requested.")
        } else {
            cli::cli_alert_success("Totally {length(res)} gene{?s} {?is/are} selected.")
        }
    }
    varFeatures(object) <- res
    return(object)
}

#' @export
#' @method selectBatchHVG ligerDataset
#' @rdname selectBatchHVG
selectBatchHVG.ligerDataset <- function(
        object,
        nGenes = 2000,
        features = NULL,
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    features <- .idxCheck(object, features, orient = "feature")
    mat <- rawData(object)
    mat <- normalize(mat, log = TRUE, scaleFactor = 1e6)
    fm <- featureMeta(object)
    idx <- rep(FALSE, nrow(object))
    idx[features] <- TRUE
    # idx <- idx & (fm$nCell > 0)
    mat <- mat[idx, , drop = FALSE]
    stats <- selectBatchHVG(mat, nGenes = nGenes, verbose = verbose,
                            returnStats = TRUE, ...)
    fm$means <- 0
    fm$dispersion <- 0
    fm$dispersion_norm <- 0
    fm$highly_variable <- FALSE
    fm[stats$feature, "means"] <- stats$means
    fm[stats$feature, "dispersion"] <- stats$dispersion
    fm[stats$feature, "dispersion_norm"] <- stats$dispersion_norm
    fm[stats$feature, "highly_variable"] <- stats$highly_variable
    featureMeta(object) <- fm
    return(object)
}

#' @export
#' @method selectBatchHVG dgCMatrix
#' @rdname selectBatchHVG
selectBatchHVG.dgCMatrix <- function(
        object,
        nGenes = 2000,
        returnStats = FALSE,
        verbose = getOption("ligerVerbose", TRUE),
        ...
) {
    # Get mean and var
    idx <- Matrix::rowSums(object > 0) > 0
    object <- object[idx, , drop = FALSE]
    means <- Matrix::rowMeans(object)
    vars <- rowVars_sparse_rcpp(object, means)
    stats <- .selectBatchHVG.by.metric(
        feature = rownames(object),
        means = means,
        vars = vars,
        nGenes = nGenes,
        verbose = verbose
    )
    if (isTRUE(returnStats)) {
        return(stats)
    } else {
        stats %>%
            dplyr::arrange(-.data[['dispersion_norm']]) %>%
            dplyr::filter(.data[['highly_variable']]) %>%
            dplyr::pull(.data[['feature']])
    }
}

.selectBatchHVG.by.metric <- function(
        feature,
        means,
        vars,
        nGenes,
        verbose
) {
    means[means == 0] <- 1e-12
    dispersion <- vars / means
    df <- data.frame(
        feature = feature,
        means = means,
        dispersion = dispersion
    )
    bins <- c(-Inf, stats::quantile(means, seq(0.1, 1, 0.05)), Inf)
    df$mean_bin <- cut(means, breaks = bins)
    df <- df %>%
        dplyr::group_by(.data[['mean_bin']]) %>%
        dplyr::mutate(
            dispersion_norm = (.data[['dispersion']] - stats::median(.data[['dispersion']]))/.mad(.data[['dispersion']])
        )
    if (nGenes > nrow(df)) nGenes <- nrow(df)
    df$highly_variable <- FALSE
    df$highly_variable[order(df$dispersion_norm, decreasing = TRUE)[seq(nGenes)]] <- TRUE
    df
}

.mad <- function(x) {
    stats::median(abs(x - stats::median(x)))/0.6744897501960817
}
