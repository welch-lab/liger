#' @title Find DEG between two groups
#' @description Find DEG between two groups. Two methods are supported:
#' \code{"wilcoxon"} and \code{"pseudoBulk"}. Wilcoxon rank sum test is
#' performed on single-cell level, while pseudo-bulk method aggregates cells
#' basing on biological replicates and calls bulk RNAseq DE methods, DESeq2 wald
#' test. When real biological replicates are not available, pseudo replicates
#' can be generated. Please see below for detailed scenario usage.
#' @section Pairwise DEG Scenarios:
#' Users can select classes of cells from a variable in \code{cellMeta}.
#' \code{variable1} and \code{variable2} are used to specify a column in
#' \code{cellMeta}, and \code{groupTest} and \code{groupCtrl} are used to specify
#' existing classes from \code{variable1} and \code{variable2}, respectively.
#' When \code{variable2} is missing, \code{groupCtrl} will be considered from
#' \code{variable1}.
#'
#' For example, when \code{variable1 = "celltype"} and \code{variable2 = NULL},
#' \code{groupTest} and \code{groupCtrl} should be valid cell types in
#' \code{object$celltype}.
#'
#' When \code{variable1} is "celltype" and \code{variable2} is "gender",
#' \code{groupTest} should be a valid cell type from \code{object$celltype} and
#' \code{groupCtrl} should be a valid class from \code{object$gender}.
#'
#' When both \code{variable1} and \code{variable2} are missing, \code{groupTest}
#' and \code{groupCtrl} should be valid index of cells in \code{object}.
#' @param object A \linkS4class{liger} object, with normalized data available
#' @param groupTest,groupCtrl,variable1,variable2 Condition specification. See
#' \code{?runPairwiseDEG} section \bold{Pairwise DEG Scenarios} for detail.
#' @param method DEG test method to use. Choose from \code{"wilcoxon"} or
#' \code{"pseudoBulk"}. Default \code{"wilcoxon"}
#' @param usePeak Logical. Whether to use peak count instead of gene count.
#' Only supported when ATAC datasets are involved. Default \code{FALSE}.
#' @param useReplicate \code{cellMeta} variable of biological replicate
#' annotation. Only used with \code{method = "pseudoBulk"}. Default \code{NULL}
#' will create \code{nPsdRep} pseudo replicates per group.
#' @param nPsdRep Number of pseudo replicates to create. Only used when
#' \code{method = "pseudoBulk", useReplicate = NULL}. Default \code{5}.
#' @param minCellPerRep Numeric, will not make pseudo-bulk for replicate with
#' less than this number of cells. Default \code{10}.
#' @param seed Random seed to use for pseudo-replicate generation. Default
#' \code{1}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @return A data.frame with DEG information
#' @rdname liger-DEG
#' @export
#' @examples
#' # Compare between cluster "0" and cluster "1"
#' degStats <- runPairwiseDEG(pbmcPlot, groupTest = 0, groupCtrl = 1,
#'                            variable1 = "leiden_cluster")
#' # Compare between all cells from cluster "5" and
#' # all cells from dataset "stim"
#' degStats <- runPairwiseDEG(pbmcPlot, groupTest = "5", groupCtrl = "stim",
#'                            variable1 = "leiden_cluster",
#'                            variable2 = "dataset")
runPairwiseDEG <- function(
        object,
        groupTest,
        groupCtrl,
        variable1 = NULL,
        variable2 = NULL,
        method = c("wilcoxon", "pseudoBulk"),
        usePeak = FALSE,
        useReplicate = NULL,
        nPsdRep = 5,
        minCellPerRep = 10,
        seed = 1,
        verbose = getOption("ligerVerbose", TRUE)
) {
    method <- match.arg(method)
    if (is.null(variable1) && is.null(variable2)) {
        # Directly using cell index
        groups <- list(
            .idxCheck(object, groupTest, "cell"),
            .idxCheck(object, groupCtrl, "cell")
        )
        group1Name <- "test"
        group2Name <- "control"
        names(groups) <- c("test", "control")
    } else if (!is.null(variable1)) {
        var1 <- .fetchCellMetaVar(object, variable1,
                                  checkCategorical = TRUE, drop = TRUE,
                                  droplevels = TRUE)
        group1Idx <- which(var1 %in% groupTest)
        group1Name <- paste(groupTest, collapse = ".")
        if (is.null(variable2)) {
            variable2 <- variable1
            var2 <- var1
        } else {
            var2 <- .fetchCellMetaVar(object, variable2,
                                           checkCategorical = TRUE, drop = TRUE,
                                           droplevels = TRUE)
        }
        group2Idx <- which(var2 %in% groupCtrl)
        group2Name <- paste(groupCtrl, collapse = ".")
        groups <- list(group1Idx, group2Idx)
        names(groups) <- c(group1Name, group2Name)
    } else {
        cli::cli_abort("Please see {.code ?runPairwiseDEG} for usage.")
    }
    result <- .runDEG(object, groups = groups, method = method,
                      usePeak = usePeak, useReplicate = useReplicate,
                      nPsdRep = nPsdRep,
                      minCellPerRep = minCellPerRep, seed = seed,
                      verbose = verbose)
    result <- result[result$group == group1Name,]
    attributes(result)$meta <- list(
        groupTest = groupTest,
        variable1 = variable1,
        groupCtrl = groupCtrl,
        variable2 = variable2
    )
    return(result)
}

#' @rdname liger-DEG
#' @export
#' @param conditionBy \code{cellMeta} variable(s). Marker detection will be
#' performed for each level of this variable. Multiple variables will be
#' combined. Default \code{NULL} uses default cluster.
#' @param splitBy Split data by \code{cellMeta} variable(s) here and identify
#' markers for \code{conditionBy} within each chunk. Default \code{NULL}.
#' @param useDatasets Datasets to perform marker detection within. Default
#' \code{NULL} will use all datasets.
#' @section Marker Detection Scenarios:
#' Marker detection is generally performed in a one vs. rest manner. The
#' grouping of such condition is specified by \code{conditionBy}, which should
#' be a column name in \code{cellMeta}. When \code{splitBy} is specified as
#' another variable name in \code{cellMeta}, the marker detection will be
#' iteratively done for within each level of \code{splitBy} variable.
#'
#' For example, when \code{conditionBy = "celltype"} and \code{splitBy = NULL},
#' marker detection will be performed by comparing all cells of "celltype_i"
#' against all other cells, and etc.
#'
#' When \code{conditionBy = "celltype"} and \code{splitBy = "gender"}, marker
#' detection will be performed by comparing "celltype_i" cells from "gender_j"
#' against other cells from "gender_j", and etc.
#' @examples
#' # Identify markers for each cluster. Equivalent to old version
#' # `runWilcoxon(method = "cluster")`
#' markerStats <- runMarkerDEG(pbmcPlot, conditionBy = "leiden_cluster")
#' # Identify dataset markers within each cluster. Equivalent to old version
#' # `runWilcoxon(method = "dataset")`.
#' markerStatsList <- runMarkerDEG(pbmcPlot, conditionBy = "dataset",
#'                                 splitBy = "leiden_cluster")
runMarkerDEG <- function(
        object,
        conditionBy = NULL,
        splitBy = NULL, # The previous by dataset strategy
        method = c("wilcoxon", "pseudoBulk"),
        useDatasets = NULL,
        usePeak = FALSE,
        useReplicate = NULL,
        nPsdRep = 5,
        minCellPerRep = 10,
        seed = 1,
        verbose = getOption("ligerVerbose", TRUE)
) {
    useDatasets <- .checkUseDatasets(object, useDatasets)
    allCellIdx <- seq(ncol(object))[object$dataset %in% useDatasets]
    conditionBy <- conditionBy %||% object@uns$defaultCluster
    if (is.null(conditionBy)) {
        cli::cli_abort("No {.var conditionBy} given or default cluster not set.")
    }
    conditionBy <- .fetchCellMetaVar(
        object, conditionBy, cellIdx = allCellIdx,
        checkCategorical = TRUE, drop = FALSE, droplevels = TRUE
    )
    conditionBy <- interaction(conditionBy, drop = TRUE)
    splitBy <- .fetchCellMetaVar(
        object, splitBy, cellIdx = allCellIdx,
        checkCategorical = TRUE, drop = FALSE, droplevels = TRUE
    )
    splitBy <- interaction(splitBy, drop = TRUE)
    if (nlevels(splitBy) <= 1) {
        groups <- split(allCellIdx, conditionBy)
        result <- .runDEG(object, groups = groups, method = method,
                          usePeak = usePeak, useReplicate = useReplicate,
                          nPsdRep = nPsdRep,
                          minCellPerRep = minCellPerRep,
                          seed = seed, verbose = verbose)
    } else {
        result <- list()
        for (i in seq_along(levels(splitBy))) {
            subIdx <- splitBy == levels(splitBy)[i]
            subCellIdx <- allCellIdx[subIdx]
            groups <- split(subCellIdx, conditionBy[subIdx])
            result[[levels(splitBy)[i]]] <- .runDEG(
                object, groups = groups, method = method, usePeak = usePeak,
                useReplicate = useReplicate, nPsdRep = nPsdRep,
                minCellPerRep = minCellPerRep, seed = seed, verbose = verbose
            )
        }
    }

    return(result)
}

#' @rdname liger-DEG
#' @export
#' @param data.use Same as \code{useDatasets}.
#' @param compare.method Choose from \code{"clusters"} (default) or
#' \code{"datasets"}. \code{"clusters"} compares each cluster against all other
#' cells, while \code{"datasets"} run within each cluster and compare each
#' dataset against all other datasets.
runWilcoxon <- function(
        object,
        data.use = NULL,
        compare.method = c("clusters", "datasets")
) {
    lifecycle::deprecate_warn(
        "1.99.0", "runWilcoxon()",
        details = "Please use `runMarkerDEG()` with `method = 'wilcoxon'` instead."
    )
    compare.method <- match.arg(compare.method)
    if (compare.method == "clusters") {
        res <- runMarkerDEG(object, conditionBy = object@uns$defaultCluster,
                            splitBy = NULL, method = "wilcoxon")
    } else if (compare.method == "datasets") {
        res <- runMarkerDEG(object, conditionBy = "dataset",
                            splitBy = object@uns$defaultCluster,
                            method = "wilcoxon")
    }
    return(res)
}

# groups - As underlying function, this must be organized into list of numeric
# cell index vectors.
.runDEG <- function(
        object,
        groups,
        method = c("wilcoxon", "pseudoBulk"),
        # byDataset = FALSE,
        usePeak = FALSE,
        useReplicate = NULL,
        nPsdRep = 5,
        minCellPerRep = 10,
        seed = 1,
        verbose = getOption("ligerVerbose", TRUE)
) {
    method <- match.arg(method)
    allCellIdx <- unlist(groups)
    if (length(allCellIdx) == 0)
        cli::cli_abort(c(x = "No cell selected"))
    allCellBC <- colnames(object)[allCellIdx]
    datasetInvolve <- levels(object$dataset[allCellIdx, drop = TRUE])
    var <- factor(rep(names(groups), lengths(groups)), levels = names(groups))
    if (isTRUE(usePeak)) {
        useDatasets <- .checkUseDatasets(object, useDatasets = datasetInvolve,
                                         modal = "atac")
    } else {
        useDatasets <- .checkUseDatasets(object, useDatasets = datasetInvolve)
    }
    slot <- .DE.checkDataAvail(object, datasetInvolve, method, usePeak)
    dataList <- getMatrix(object, slot, datasetInvolve, returnList = TRUE)
    features <- Reduce(intersect, lapply(dataList, rownames))
    dataList <- lapply(dataList, function(x) x[features, , drop = FALSE])

    mat <- Reduce(cbind, dataList)
    mat <- mat[, allCellBC, drop = FALSE]
    if (method == "wilcoxon") {
        cliID <- cli::cli_process_start("Running Wilcoxon rank-sum test")
        mat <- log1p(1e10*mat)
        result <- wilcoxauc(mat, var)
        cli::cli_process_done(id = cliID)
    } else if (method == "pseudoBulk") {
        if (is.null(useReplicate)) {
            replicateAnn <- setupPseudoRep(var, nRep = nPsdRep,
                                            seed = seed)
        } else {
            replicateAnn <- .fetchCellMetaVar(
                object, useReplicate,
                cellIdx = allCellIdx,
                drop = FALSE,
                checkCategorical = TRUE,
                droplevels = TRUE
            )
            replicateAnn$groups <- var
        }

        pbs <- makePseudoBulk(mat, replicateAnn,
                               minCellPerRep = minCellPerRep,
                               verbose = verbose)
        pb <- pbs[[1]]
        replicateAnn <- pbs[[2]]
        var <- sapply(levels(replicateAnn$groups), function(x) {
            nlevels(interaction(replicateAnn[replicateAnn$groups == x, , drop = FALSE],
                                drop = TRUE))
        })
        var <- factor(rep(names(var), var), levels = names(var))
        result <- .callDESeq2(pb, var, verbose)
    }
    return(result)
}

.DE.checkDataAvail <- function(object, useDatasets, method, usePeak) {
    if (isH5Liger(object, useDatasets)) { # nocov start
        cli::cli_abort(
            c("HDF5 based datasets detected but is not supported. ",
              "i" = "Try {.code object.sub <- downsample(object, useSlot = 'normData')} to create another object with in memory data")
        )
    } # nocov end
    if (method == "wilcoxon") {
        slot <- ifelse(usePeak, "normPeak", "normData")
    } else if (method == "pseudoBulk") {
        if (!requireNamespace("DESeq2", quietly = TRUE)) # nocov start
            cli::cli_abort(
                "Package {.pkg DESeq2} is needed for this function to work.
                Please install it by command:
                {.code BiocManager::install('DESeq2')}"
            ) # nocov end
        slot <- ifelse(usePeak, "rawPeak", "rawData")
    }
    allAvail <- all(sapply(useDatasets, function(d) {
        ld <- dataset(object, d)
        !is.null(methods::slot(ld, slot))
    }))
    if (!allAvail)
        cli::cli_abort(
            c("{.field {slot}} not all available for involved datasets: {.val {useDatasets}}",
              "i" = "{.code method = '{method}'}; {.code usePeak = {usePeak}}")
        )
    return(slot)
}

###################### Pseudo-bulk Method helper ###############################

setupPseudoRep <- function(groups, nRep = 3, seed = 1) {
    # The output data.frame should be cell per row by variable per col
    set.seed(seed)
    psdRep <- c()
    for (i in seq_along(levels(groups))) {
        groupSize <- sum(groups == levels(groups)[i])
        repVar <- sample(seq_len(groupSize) %% nRep) + 1 + (i - 1)*nRep
        psdRep <- c(psdRep, repVar)
    }
    return(data.frame(
        groups = groups,
        pseudoRep = factor(psdRep)
    ))
}

makePseudoBulk <- function(mat, replicateAnn, minCellPerRep, verbose = TRUE) {
    # mat - Extracted and contatenated matrix. intersection of genes by
    #       c(groupTest, groupCtrl) cells
    # groups - list of groups
    # replicateAnn - data.frame of replicate annotation, with rownames as
    #                barcodes and columns as variables

    # Check whether enough replicates per condition
    for (gr in levels(replicateAnn$groups)) {
        subrep <- replicateAnn[replicateAnn$groups == gr,]
        splitLabel <- interaction(subrep, drop = TRUE)
        if (nlevels(splitLabel) < 2) {
            cli::cli_abort(
                c("Too few replicates for condition {.val {gr}}. Cannot create pseudo-bulks.",
                  "i" = "Please consider creating pseudo-replicates or using {.code method = 'wilcoxon'} instead.")
            )
        }
    }
    splitLabel <- interaction(replicateAnn, drop = TRUE)

    labelCounts <- table(splitLabel)
    ignored <- names(labelCounts)[labelCounts < minCellPerRep]
    keep <- names(labelCounts)[labelCounts >= minCellPerRep]
    idx <- splitLabel %in% keep
    splitLabel <- splitLabel[idx, drop = TRUE]
    mat <- mat[, idx, drop = FALSE]
    replicateAnn <- replicateAnn[idx, , drop = FALSE]
    if (verbose) {
        if (length(ignored) > 0) cli::cli_alert_warning("Ignoring replicates with too few cells: {.val {ignored}}")
        cli::cli_alert_info("Replicate sizes:")
        print(table(splitLabel))
    }
    pseudoBulks <- colAggregateSums_sparse(mat, as.integer(splitLabel) - 1,
                                           nlevels(splitLabel))
    dimnames(pseudoBulks) <- list(rownames(mat), levels(splitLabel))
    pseudoBulks <- pseudoBulks[rowSums(pseudoBulks) > 0,]
    return(list(pseudoBulks, replicateAnn))
}

.callDESeq2 <- function(pseudoBulks, groups,
                         verbose = getOption("ligerVerbose", TRUE)) {
    # DESeq2 workflow
    if (isTRUE(verbose)) cliID <- cli::cli_process_start("Calling DESeq2 Wald test")
    ## NOTE: DESeq2 wishes that the contrast/control group is the first level
    ## whereas we required it as the second in upstream input. So we need to
    ## reverse it here.
    groups <- stats::relevel(groups, ref = levels(groups)[2])
    ## Now levels(groups)[1] becomes control and levels(groups)[2] becomes
    ## the test group
    des <- DESeq2::DESeqDataSetFromMatrix(
        countData = pseudoBulks,
        colData = data.frame(groups = groups),
        design = stats::formula("~groups")
    )
    des <- DESeq2::DESeq(des, test = "Wald", quiet = !verbose)
    res <- DESeq2::results(des, contrast = c("groups", levels(groups)[2],
                                             levels(groups)[1]))
    res <- .DataFrame.as.data.frame(res)
    res$feature <- rownames(res)
    rownames(res) <- NULL
    res$group <- levels(groups)[2]
    res <- res[, c(7, 8, 2, 5, 6)]
    colnames(res) <- c("feature", "group", "logFC", "pval", "padj")
    if (isTRUE(verbose)) cli::cli_process_done(id = cliID)
    return(res)
}


####################### Wilcoxon rank-sum test helper ##########################

# X: matrix of data to be tested
# y: grouping label of columns of X
# Rcpp source code located in src/wilcoxon.cpp
wilcoxauc <- function(x, clusterVar) {
    if (methods::is(x, 'dgTMatrix')) x <- methods::as(x, 'CsparseMatrix') # nocov start
    if (methods::is(x, 'TsparseMatrix')) x <- methods::as(x, 'CsparseMatrix')
    if (is.null(row.names(x))) {
        rownames(x) <- paste0('Feature', seq(nrow(x)))
    } # nocov end
    if (!is.factor(clusterVar)) clusterVar <- factor(clusterVar)
    clusterVar <- droplevels(clusterVar)
    groupSize <- as.numeric(table(clusterVar))

    ## Compute primary statistics
    n1n2 <- groupSize * (ncol(x) - groupSize)
    # rankRes - list(X_ranked, ties), where X_ranked is obs x feature
    xRanked <- Matrix::t(x)
    # This computes the ranking of non-zero values and the ties
    ties <- cpp_rank_matrix_dgc(xRanked@x, xRanked@p,
                                nrow(xRanked), ncol(xRanked))
    # ranksRes <- list(X_ranked = xT, ties = ties)

    # rankRes <- colRanking(x)
    ustat <- computeUstat(xRanked, clusterVar, n1n2, groupSize)
    auc <- t(ustat / n1n2)
    pvals <- computePval(ustat, ties, ncol(x), n1n2)
    fdr <- apply(pvals, 2, function(p) stats::p.adjust(p, 'BH'))

    ### Auxiliary Statistics (AvgExpr, PctIn, LFC, etc)
    groupSums <- colAggregateSum_sparse(x, as.integer(clusterVar) - 1, length(unique(clusterVar)))
    # groupSums <- colAggregateSum(x, clusterVar)
    group_nnz <- colNNZAggr_sparse(x, as.integer(clusterVar) - 1, length(unique(clusterVar)))
    # group_nnz <- colNNZAggr(x, clusterVar)
    group_pct <- t(sweep(group_nnz, 1, as.numeric(table(clusterVar)), "/"))

    group_pct_out <- sweep(-group_nnz, 2, colSums(group_nnz), "+")
    group_pct_out <- sweep(group_pct_out, 1,
                           as.numeric(length(clusterVar) - table(clusterVar)),
                           "/")
    group_pct_out <- t(group_pct_out)

    groupMeans <- t(sweep(groupSums, 1, as.numeric(table(clusterVar)), "/"))

    cs <- colSums(groupSums)
    gs <- as.numeric(table(clusterVar))
    lfc <- Reduce(cbind, lapply(seq_along(levels(clusterVar)), function(g) {
        groupMeans[, g] - (cs - groupSums[g, ])/(length(clusterVar) - gs[g])
    }))

    data.frame(
        feature = rep(row.names(x), times = length(levels(clusterVar))),
        group = factor(rep(levels(clusterVar), each = nrow(x)),
                       levels = levels(clusterVar)),
        avgExpr = as.numeric(groupMeans),
        logFC = as.numeric(lfc),
        statistic = as.numeric(t(ustat)),
        auc = as.numeric(auc),
        pval = as.numeric(pvals),
        padj = as.numeric(fdr),
        pct_in = as.numeric(100 * group_pct),
        pct_out = as.numeric(100 * group_pct_out)
    )
}

computeUstat <- function(Xr, cols, n1n2, groupSize) {
    grs <- rowAggregateSum_sparse(Xr, as.integer(cols) - 1, length(unique(cols)))
    # grs <- rowAggregateSum(Xr, cols)

    # if (inherits(Xr, 'dgCMatrix')) {
    # With the ranking of only non-zero features, here the tie-ranking of
    # zeros need to be added.
    nnz <- rowNNZAggr_sparse(Xr, as.integer(cols) - 1, length(unique(cols)))
    gnz <- groupSize - nnz
    zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
    ustat <- t((t(gnz) * zero.ranks)) + grs - groupSize*(groupSize + 1)/2
    # } else {
    #     ustat <- grs - groupSize * (groupSize + 1) / 2
    # }
    return(ustat)
}

computePval <- function(ustat, ties, N, n1n2) {
    z <- ustat - .5 * n1n2
    z <- z - sign(z) * .5
    .x1 <- N ^ 3 - N
    .x2 <- 1 / (12 * (N ^ 2 - N))
    rhs <- unlist(lapply(ties, function(tvals) {
        (.x1 - sum(tvals ^ 3 - tvals)) * .x2
    }))
    usigma <- sqrt(matrix(n1n2, ncol = 1) %*% matrix(rhs, nrow = 1))
    z <- t(z / usigma)
    pvals <- matrix(2 * stats::pnorm(-abs(as.numeric(z))), ncol = ncol(z))
    return(pvals)
}





######################## Visualization #########################################

#' Create heatmap for showing top marker expression in conditions
#' @export
#' @param object A \linkS4class{liger} object, with normalized data and metadata
#' to annotate available.
#' @param result The data.frame returned by \code{\link{runMarkerDEG}}.
#' @param topN Number of top features to be plot for each group. Default
#' \code{5}.
#' @param lfcThresh Hard threshold on logFC value. Default \code{1}.
#' @param padjThresh Hard threshold on adjusted P-value. Default \code{0.05}.
#' @param pctInThresh,pctOutThresh Threshold on expression percentage. These
#' mean that a feature will only pass the filter if it is expressed in more than
#' \code{pctInThresh} percent of cells in the corresponding cluster. Similarly
#' for \code{pctOutThresh}. Default \code{50} percent for both.
#' @param dedupBy When ranking by padj and logFC and a feature is ranked as top
#' for multiple clusters, assign this feature as the marker of a cluster when
#' it has the largest \code{"logFC"} in the cluster or has the lowest
#' \code{"padj"}. Default \code{"logFC"}.
#' @param groupBy Cell metadata variable names for cell grouping. Downsample
#' balancing will also be aware of this. Default \code{c("dataset",
#' "leiden_cluster")}.
#' @param groupSize Maximum number of cells in each group to be downsampled for
#' plotting. Default \code{50}.
#' @param column_title Title on the column. Default \code{NULL}.
#' @param ... Parameter passed to wrapped functions in the inheritance order:
#' \code{\link{plotGeneHeatmap}}, \code{\link{.plotHeatmap}},
#' \code{ComplexHeatmap::\link[ComplexHeatmap]{Heatmap}}
#' @examples
#' markerTable <- runMarkerDEG(pbmcPlot)
#' plotMarkerHeatmap(pbmcPlot, markerTable)
plotMarkerHeatmap <- function(
        object,
        result,
        topN = 5,
        lfcThresh = 1,
        padjThresh = 0.05,
        pctInThresh = 50,
        pctOutThresh = 50,
        dedupBy = c("logFC", "padj"),
        groupBy = NULL,
        groupSize = 50,
        column_title = NULL,
        ...
) {
    dedupBy <- match.arg(dedupBy)
    if (dedupBy == "logFC") {
        result <- result[order(result[[dedupBy]], decreasing = TRUE), ]
    } else if (dedupBy == "padj") {
        result <- result[order(result[[dedupBy]], decreasing = FALSE), ]
    }
    pctInThresh <- pctInThresh %||% 0
    pctOutThresh <- pctOutThresh %||% 100
    groupBy <- groupBy %||% c("dataset", object@uns$defaultCluster)
    result <- result[!duplicated(result$feature), ]
    result <- result %>% dplyr::filter(.data$logFC > lfcThresh,
                                       .data$padj < padjThresh,
                                       .data$pct_in > pctInThresh,
                                       .data$pct_out < pctOutThresh) %>%
        dplyr::group_by(.data[["group"]]) %>%
        dplyr::arrange(.data[["padj"]], -.data[["logFC"]], .by_group = TRUE) %>%
        dplyr::filter(dplyr::row_number() %in% seq(topN)) %>%
        as.data.frame()
    cellIdx <- downsample(object, maxCells = groupSize, balance = groupBy,
                          returnIndex = TRUE)
    featureAnn <- result[, "group", drop = FALSE]

    rownames(featureAnn) <- result$feature
    colnames(featureAnn) <- "marker"
    plotGeneHeatmap(object, features = result$feature,
                    cellIdx = cellIdx,
                    useCellMeta = groupBy,
                    featureAnnotation = featureAnn,
                    cellSplitBy = rev(groupBy),
                    featureSplitBy = "marker",
                    showFeatureLegend = FALSE,
                    cluster_columns = FALSE,
                    cluster_column_slices = FALSE,
                    cluster_rows = FALSE,
                    cluster_row_slices = FALSE,
                    column_title = column_title,
                    ...)
}
