#' Perform Wilcoxon rank-sum test
#' @description Perform Wilcoxon rank-sum tests on specified dataset using
#' given method.
#' @param object A \linkS4class{liger} object with cluster assignment available.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets to be normalized. Default
#' \code{NULL} uses all datasets.
#' @param method Choose from \code{"clusters"} or \code{"datasets"}. Default
#' \code{"clusters"} compares between clusters across all datasets, while
#' \code{"datasets"} compares between datasets within each cluster.
#' @param useCluster The name of the column in \code{cellMeta} slot storing the
#' cluster assignment variable. Default \code{"leiden_cluster"}
#' @param usePeak Logical, whether to test peak counts instead of gene
#' expression. Requires presence of ATAC modility datasets. Default
#' \code{FALSE}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param data.use,compare.method \bold{Deprecated}. See Usage section for
#' replacement.
#' @return A 10-columns data.frame with test results.
#' @export
#' @examples
#' library(dplyr)
#' result <- runWilcoxon(pbmcPlot)
#' result %>% group_by(group) %>% top_n(2, logFC)
runWilcoxon <- function(
        object,
        useDatasets = NULL,
        method = c("clusters", "datasets"),
        useCluster = "leiden_cluster",
        usePeak = FALSE,
        verbose = getOption("ligerVerbose"),
        # Deprecated coding style,
        data.use = useDatasets,
        compare.method = method
) {
    .deprecateArgs(list(data.use = "useDatasets", compare.method = "method"))
    method <- match.arg(method)
    # Input checks ####
    useDatasets <- .checkUseDatasets(object, useDatasets,
                                     modal = ifelse(usePeak, "atac", "default"))
    if (!isTRUE(usePeak)) {
        if (method == "datasets" & length(useDatasets) < 2)
            stop("Should have at least 2 datasets as input ",
                 "when compare between datasets")
        if (isH5Liger(object, useDatasets)) {
            stop("HDF5 based datasets detected but is not supported. \n",
                 "Try `object.sub <- downsample(object, useSlot = ",
                 "'normData')` to create ANOTHER object with in memory data.")
        }
        allNormed <- all(sapply(datasets(object),
                                function(ld) !is.null(normData(ld))))
        if (!allNormed)
            stop("All datasets being involved has to be normalized")

        ## get all shared genes of datasets involved
        normDataList <- getMatrix(object, "normData", dataset = useDatasets,
                                  returnList = TRUE)
        features <- Reduce(intersect, lapply(normDataList, rownames))
        normDataList <- lapply(normDataList, function(x) x[features,])
        featureMatrix <- Reduce(cbind, normDataList)
    } else {
        if (method == "datasets" || length(useDatasets) != 1)
            stop("For wilcoxon test on peak counts, can only use ",
                 "\"cluster\" method on one dataset.")
        normPeakList <- lapply(useDatasets, function(d) normPeak(object, d))
        features <- Reduce(intersect, lapply(normPeakList, rownames))
        featureMatrix <- Reduce(cbind, normPeakList)
        #featureMatrix <- normPeak(object, useDatasets)
        if (is.null(featureMatrix))
            stop("Peak counts of specified dataset has to be normalized. ",
                 "Please try `normalizePeak(object, useDatasets = '",
                 useDatasets, "')`.")
        #features <- rownames(featureMatrix)
    }

    ## Subset metadata involved
    cellIdx <- object$dataset %in% useDatasets
    cellSource <- object$dataset[cellIdx]
    clusters <- .fetchCellMetaVar(object, useCluster, cellIdx = cellIdx,
                                  checkCategorical = TRUE)

    if (isTRUE(verbose))
        .log("Performing Wilcoxon test on ", length(useDatasets), " datasets: ",
             paste(useDatasets, collapse = ", "))
    # perform wilcoxon test ####
    if (method == "clusters") {
        # compare between clusters across datasets
        nfeatures <- length(features)
        if (nfeatures > 100000) {
            if (isTRUE(verbose)) .log("Calculating Large-scale Input...")
            results <- Reduce(rbind, lapply(
                suppressWarnings(split(seq(nfeatures),
                                       seq(nfeatures / 100000))),
                function(index) {
                    wilcoxauc(log(featureMatrix[index, ] + 1e-10), clusters)
                }))
        } else {
            # TODO: If we add log-transformation to normalization method in the
            # future, remember to have conditions here.
            results <- wilcoxauc(log(featureMatrix + 1e-10), clusters)
        }
    } else {
        # compare between datasets within each cluster
        results <- Reduce(rbind, lapply(levels(clusters), function(cluster) {
            clusterIdx <- clusters == cluster
            subLabel <- paste0(cluster, "-", cellSource[clusterIdx])
            if (length(unique(subLabel)) == 1) {
                # if cluster has only 1 data source
                warning("Skipped Cluster ", cluster,
                        " since it has only one dataset source.")
                return()
            } else {
                subMatrix <- featureMatrix[, clusterIdx]
                return(wilcoxauc(log(subMatrix + 1e-10), subLabel))
            }
        }))
    }
    return(results)
}

#' Find shared and dataset-specific markers
#' @description Applies various filters to genes on the shared (\eqn{W}) and
#' dataset-specific (\eqn{V}) components of the factorization, before selecting
#' those which load most significantly on each factor (in a shared or
#' dataset-specific way).
#' @param object \linkS4class{liger} object with factorization results.
#' @param dataset1 Name of first dataset. Required.
#' @param dataset2 Name of second dataset. Required
#' @param factorShareThresh Numeric. Only factors with a dataset specificity
#' less than or equal to this threshold will be used. Default \code{10}.
#' @param datasetSpecificity Numeric vector. Pre-calculated dataset specificity
#' if available. Length should match number of all factors available. Default
#' \code{NULL} automatically calculates with
#' \code{\link{calcDatasetSpecificity}}.
#' @param logFCThresh Numeric. Lower log-fold change threshold for differential
#' expression in markers. Default \code{1}.
#' @param pvalThresh Numeric. Upper p-value threshold for Wilcoxon rank test for
#' gene expression. Default \code{0.05}.
#' @param nGenes Integer. Max number of genes to report for each dataset.
#' Default \code{30}.
#' @param printGenes Logical. Whether to print ordered markers passing logFC,
#' UMI and frac thresholds, when \code{verbose = TRUE}. Default \code{FALSE}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} which is \code{TRUE} if users have not set.
#' @param factor.share.thresh,dataset.specificity,log.fc.thresh,pval.thresh,num.genes,print.genes
#' \bold{Deprecated}. See Usage section for replacement.
#' @return A list object consisting of the following entries:
#' \item{[value of `dataset1`]}{data.frame of dataset1-specific markers}
#' \item{shared}{data.frame of shared markers}
#' \item{[value of `dataset1`]}{data.frame of dataset2-specific markers}
#' \item{num_factors_V1}{A frequency table indicating the number of factors each
#' marker appears, in dataset1}
#' \item{num_factors_V2}{A frequency table indicating the number of factors each
#' marker appears, in dataset2}
#' @export
#' @examples
#' library(dplyr)
#' result <- getFactorMarkers(pbmcPlot, dataset1 = "ctrl", dataset2 = "stim")
#' print(class(result))
#' print(names(result))
#' result$shared %>% group_by(factor_num) %>% top_n(2, log2fc)
getFactorMarkers <- function(
        object,
        dataset1,
        dataset2,
        factorShareThresh = 10,
        datasetSpecificity = NULL,
        logFCThresh = 1,
        pvalThresh = 0.05,
        nGenes = 30,
        printGenes = FALSE,
        verbose = getOption("ligerVerbose"),
        # Deprecated coding style
        factor.share.thresh = factorShareThresh,
        dataset.specificity = datasetSpecificity,
        log.fc.thresh = logFCThresh,
        pval.thresh = pvalThresh,
        num.genes = nGenes,
        print.genes = printGenes
) {
    .deprecateArgs(list(factor.share.thresh = "factorShareThresh",
                        dataset.specificity = "datasetSpecificity",
                        log.fc.thresh = "logFCThresh",
                        pval.thresh = "pvalThresh", num.genes = "nGenes",
                        print.genes = "printGenes"))
    dataset1 <- .checkUseDatasets(object, useDatasets = dataset1)
    dataset2 <- .checkUseDatasets(object, useDatasets = dataset2)
    if (any(isH5Liger(object, dataset = c(dataset1, dataset2))))
        stop("Please use in-memory liger object for this analysis.`")
    if (is.null(nGenes)) {
        nGenes <- length(varFeatures(object))
    }
    if (is.null(datasetSpecificity)) {
        datasetSpecificity <- calcDatasetSpecificity(object,
                                                     dataset1 = dataset1,
                                                     dataset2 = dataset2,
                                                     doPlot = FALSE)[[3]]
    }
    useFactors <- which(abs(datasetSpecificity) <= factorShareThresh)
    if (length(useFactors) == 0) {
        stop("No factor passed the dataset specificity threshold, ",
             "please try a larger `factorShareThresh`.")
    }
    if (length(useFactors) == 1 && isTRUE(verbose)) {
        warning("Only 1 factor passed the dataset specificity threshold.")
    }

    H <- getMatrix(object, "H", dataset = c(dataset1, dataset2))
    H_scaled <- lapply(H, function(x) {
        scale(t(x), scale = TRUE, center = TRUE)
    })
    labels <- list()
    for (i in seq_along(H_scaled)) {
        idx <- apply(H_scaled[[i]][, useFactors, drop = FALSE], 1, which.max)
        labels[[i]] <- useFactors[idx]
    }
    names(labels) <- names(H_scaled)

    V1_matrices <- list()
    V2_matrices <- list()
    W_matrices <- list()
    vargene <- varFeatures(object)
    if (isTRUE(verbose)) {
        .log("Performing wilcoxon test between datasets \"", dataset1,
             "\" and \"", dataset2, "\", \nbasing on factor loading.")
        if (!isTRUE(printGenes))
            pb <- utils::txtProgressBar(0, length(useFactors), style = 3)
    }
    for (j in seq_along(useFactors)) {
        i <- useFactors[j]

        W <- getMatrix(object, "W")
        V1 <- getMatrix(object, "V", dataset = dataset1)
        V2 <- getMatrix(object, "V", dataset = dataset2)

        # if not max factor for any cell in either dataset
        if (sum(labels[[dataset1]] == i) <= 1 ||
            sum(labels[[dataset2]] == i) <= 1) {
            warning("Factor ", i, " did not appear as max in ",
                    "any cell in either dataset", immediate. = TRUE)
            next
        }

        # filter genes by gene_count and cell_frac thresholds
        normDataList <- getMatrix(object, "normData",
                                  dataset = c(dataset1, dataset2))
        normData <- cbind(normDataList[[1]][vargene, labels[[1]] == i],
                          normDataList[[2]][vargene, labels[[2]] == i])
        cellLabel <- rep(c(dataset1, dataset2),
                         c(sum(labels[[1]] == i), sum(labels[[2]] == i)))
        wilcoxResult <- wilcoxauc(log(normData + 1e-10), cellLabel)

        log2fc <- wilcoxResult[wilcoxResult$group == dataset1, ]$logFC
        names(log2fc) <- wilcoxResult[wilcoxResult$group == dataset1, ]$feature
        filteredGenesV1 <- wilcoxResult[wilcoxResult$logFC > logFCThresh &
                                            wilcoxResult$pval < pvalThresh,
                                        "feature"]
        filteredGenesV2 <- wilcoxResult[-wilcoxResult$logFC > logFCThresh &
                                            wilcoxResult$pval < pvalThresh,
                                        "feature"]

        W <- pmin(W + V1, W + V2)
        V1 <- V1[filteredGenesV1, , drop = FALSE]
        V2 <- V2[filteredGenesV2, , drop = FALSE]

        topGenesV1 <- character(0)
        if (length(filteredGenesV1) > 0) {
            topGeneIdx1 <- order(V1[, i], decreasing = TRUE)[seq(nGenes)]
            topGenesV1 <- row.names(V1)[topGeneIdx1]
            topGenesV1 <- topGenesV1[!is.na(topGenesV1)]
            topGenesV1 <- topGenesV1[V1[topGenesV1, i] > 0]
        }
        topGenesV2 <- character(0)
        if (length(filteredGenesV2) > 0) {
            topGenesIdx2 <- order(V2[, i], decreasing = TRUE)[seq(nGenes)]
            topGenesV2 <- row.names(V2)[topGenesIdx2]
            topGenesV2 <- topGenesV2[!is.na(topGenesV2)]
            topGenesV2 <- topGenesV2[V2[topGenesV2, i] > 0]
        }
        topGeneIdxW <- order(W[, i], decreasing = TRUE)[seq(nGenes)]
        topGenesW <- row.names(W)[topGeneIdxW]
        topGenesW <- topGenesW[!is.na(topGenesW)]
        topGenesW <- topGenesW[W[topGenesW, i] > 0]

        if (isTRUE(verbose)) {
            if (isTRUE(printGenes)) {
                .log("Factor ", i)
                message("Dataset 1:\n",
                        paste(topGenesV1, collapse = ", "),
                        "\nShared:\n",
                        paste(topGenesW, collapse = ", "),
                        "\nDataset 2\n",
                        paste(topGenesV2, collapse = ", "), "\n")
            } else {
                utils::setTxtProgressBar(pb, j)
            }
        }

        # order is V1, V2, W
        topGenes <- list(V1 = topGenesV1, V2 = topGenesV2, W = topGenesW)
        pvals <- lapply(topGenes, function(tg) {
            wilcoxResult[wilcoxResult$feature %in% tg &
                             wilcoxResult$group == dataset1, "pval"]
        })
        # bind values in matrices
        V1_matrices[[j]] <- Reduce(cbind, list(
            rep(i, length(topGenesV1)), topGenesV1, log2fc[topGenesV1], pvals$V1
        ))
        V2_matrices[[j]] <- Reduce(cbind, list(
            rep(i, length(topGenesV2)), topGenesV2, log2fc[topGenesV2], pvals$V2
        ))
        W_matrices[[j]] <- Reduce(cbind, list(
            rep(i, length(topGenesW)), topGenesW, log2fc[topGenesW], pvals$W
        ))
    }
    V1_genes <- data.frame(Reduce(rbind, V1_matrices), stringsAsFactors = FALSE)
    V2_genes <- data.frame(Reduce(rbind, V2_matrices), stringsAsFactors = FALSE)
    W_genes <- data.frame(Reduce(rbind, W_matrices), stringsAsFactors = FALSE)
    df_cols <- c("factor_num", "gene", "log2fc", "p_value")
    outputList <- list(V1_genes, W_genes, V2_genes)
    outputList <- lapply(seq_along(outputList), function(x) {
        df <- outputList[[x]]
        colnames(df) <- df_cols
        df <- transform(df, factor_num = as.numeric(df$factor_num),
                        gene = as.character(df$gene),
                        log2fc = as.numeric(df$log2fc),
                        p_value = as.numeric(df$p_value))
        # Cutoff only applies to dataset-specific dfs
        if (x != 2) {
            df[which(df$p_value < pvalThresh), ]
        } else {
            df
        }
    })
    names(outputList) <- c(dataset1, "shared", dataset2)
    outputList[["num_factors_V1"]] <- table(outputList[[dataset1]]$gene)
    outputList[["num_factors_V2"]] <- table(outputList[[dataset2]]$gene)
    return(outputList)
}

#' Calculate a dataset-specificity score for each factor
#' @description This score represents the relative magnitude of the
#' dataset-specific components of each factor's gene loadings compared to the
#' shared components for two datasets. First, for each dataset we calculate the
#' norm of the sum of each factor's shared loadings (\eqn{W}) and
#' dataset-specific loadings (\eqn{V}). We then determine the ratio of these two
#' values and subtract from 1... TODO: finish description.
#' @param object \linkS4class{liger} object with factorization results.
#' @param dataset1 Name of first dataset. Required.
#' @param dataset2 Name of second dataset. Required.
#' @param doPlot Logical. Whether to display a barplot of dataset specificity
#' scores (by factor). Default \code{FALSE}.
#' @param do.plot \bold{Deprecated}. Use \code{doPlot} instead.
#' @return List containing three elements.
#' \item{pct1}{Vector of the norm of each metagene factor for dataset1.}
#' \item{pct2}{Vector of the norm of each metagene factor for dataset2.}
#' \item{pctSpec}{Vector of dataset specificity scores.}
#' @export
calcDatasetSpecificity <- function(
        object,
        dataset1,
        dataset2,
        doPlot = FALSE,
        do.plot = doPlot
) {
    .deprecateArgs(list(do.plot = "doPlot"))
    H1 <- getMatrix(object, slot = "H", dataset = 1)
    # V: List of two g x k matrices
    V <- getMatrix(object, slot = "V", dataset = c(dataset1, dataset2))
    W <- getMatrix(object, slot = "W")
    k <- nrow(H1)
    pct1 <- rep(0, k)
    pct2 <- rep(0, k)
    for (i in seq(k)) {
        pct1[i] <- norm(as.matrix(V[[dataset1]][,i] + W[,i]), "F")
        pct2[i] <- norm(as.matrix(V[[dataset2]][,i] + W[,i]), "F")
    }
    if (isTRUE(doPlot)) {
        graphics::barplot(
            100 * (1 - (pct1 / pct2)),
            xlab = "Factor",
            ylab = "Percent Specificity",
            main = "Dataset Specificity of Factors",
            names.arg = seq(k),
            cex.names = 0.75,
            mgp = c(2, 0.5, 0)
        ) # or possibly abs(pct1-pct2)
    }
    result <- list(pct1 = pct1,
                   pct2 = pct2,
                   pctSpec = 100 * (1 - (pct1 / pct2)))
    return(result)
}

############################# For fast Wilcoxon test ###########################
# helper function for wilcoxon tests on general variables like matrix and
# dgCMatrix
# X: matrix of data to be tested
# y: grouping label of columns of X
wilcoxauc <- function(X,
                      y,
                      groups_use = NULL,
                      verbose = getOption("ligerVerbose")) {
    ## Check and possibly correct input values
    if (methods::is(X, 'dgeMatrix'))
        X <- as.matrix(X)
    if (methods::is(X, 'data.frame'))
        X <- as.matrix(X)
    # if (methods::is(X, 'DataFrame')) X <- as.matrix(X)
    # if (methods::is(X, 'data.table')) X <- as.matrix(X)
    if (methods::is(X, 'dgTMatrix'))
        X <- methods::as(X, 'dgCMatrix')
    if (methods::is(X, 'TsparseMatrix'))
        X <- methods::as(X, 'dgCMatrix')
    if (ncol(X) != length(y))
        stop("number of columns of X does not match length of y")
    if (!is.null(groups_use)) {
        idx_use <- y %in% groups_use
        y <- y[idx_use]
        X <- X[, idx_use]
    }

    if (!is.factor(y)) y <- factor(y)
    else y <- droplevels(y)
    idx_use <- which(!is.na(y))
    if (length(idx_use) < length(y)) {
        y <- y[idx_use]
        X <- X[, idx_use]
        if (isTRUE(verbose)) message('Removing NA values from labels')
    }

    group.size <- as.numeric(table(y))
    if (length(group.size[group.size > 0]) < 2) {
        stop('Must have at least 2 groups defined.')
    }

    #     features_use <- which(apply(!is.na(X), 1, all))
    #     if (verbose & length(features_use) < nrow(X)) {
    #         message('Removing features with NA values')
    #     }
    #     X <- X[features_use, ]
    if (is.null(row.names(X))) {
        rownames(X) <- paste0('Feature', seq(nrow(X)))
    }

    ## Compute primary statistics
    n1n2 <- group.size * (ncol(X) - group.size)
    if (methods::is(X, 'dgCMatrix')) {
        rank_res <- rank_matrix(Matrix::t(X))
    } else {
        rank_res <- rank_matrix(X)
    }
    ustat <- compute_ustat(rank_res$X_ranked, y, n1n2, group.size)
    auc <- t(ustat / n1n2)
    pvals <- compute_pval(ustat, rank_res$ties, ncol(X), n1n2)
    fdr <- apply(pvals, 2, function(x)
        stats::p.adjust(x, 'BH'))

    ### Auxiliary Statistics (AvgExpr, PctIn, LFC, etc)
    group_sums <- sumGroups(X, y, 1)
    group_nnz <- nnzeroGroups(X, y, 1)
    group_pct <-
        sweep(group_nnz, 1, as.numeric(table(y)), "/") %>% t()
    group_pct_out <- -group_nnz %>%
        sweep(2, colSums(group_nnz) , "+") %>%
        sweep(1, as.numeric(length(y) - table(y)), "/") %>% t()
    group_means <-
        sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
    cs <- colSums(group_sums)
    gs <- as.numeric(table(y))
    lfc <-
        Reduce(cbind, lapply(seq_len(length(levels(
            y
        ))), function(g) {
            group_means[, g] - ((cs - group_sums[g,]) / (length(y) - gs[g]))
        }))

    data.frame(
        feature = rep(row.names(X), times = length(levels(y))),
        group = factor(rep(levels(y), each = nrow(X)), levels = levels(y)),
        avgExpr = as.numeric(group_means),
        logFC = as.numeric(lfc),
        statistic = as.numeric(t(ustat)),
        auc = as.numeric(auc),
        pval = as.numeric(pvals),
        padj = as.numeric(fdr),
        pct_in = as.numeric(100 * group_pct),
        pct_out = as.numeric(100 * group_pct_out)
    )
}

compute_ustat <- function(Xr, cols, n1n2, group.size) {
    grs <- sumGroups(Xr, cols)

    if (methods::is(Xr, 'dgCMatrix')) {
        gnz <- (group.size - nnzeroGroups(Xr, cols))
        zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
        ustat <- t((t(gnz) * zero.ranks)) + grs - group.size *
            (group.size + 1) / 2
    } else {
        ustat <- grs - group.size * (group.size + 1) / 2
    }
    return(ustat)
}

compute_pval <- function(ustat, ties, N, n1n2) {
    z <- ustat - .5 * n1n2
    z <- z - sign(z) * .5
    .x1 <- N ^ 3 - N
    .x2 <- 1 / (12 * (N ^ 2 - N))
    rhs <- lapply(ties, function(tvals) {
        (.x1 - sum(tvals ^ 3 - tvals)) * .x2
    }) %>% unlist
    usigma <- sqrt(matrix(n1n2, ncol = 1) %*% matrix(rhs, nrow = 1))
    z <- t(z / usigma)
    pvals <- matrix(2 * stats::pnorm(-abs(as.numeric(z))), ncol = ncol(z))
    return(pvals)
}


#' rank_matrix
#'
#' Utility function to rank columns of matrix
#'
#' @param X feature by observation matrix.
#'
#' @return List with 2 items
#' @noRd
rank_matrix <- function(X) {
    UseMethod('rank_matrix')
}

rank_matrix.dgCMatrix <- function(X) {
    Xr <- Matrix::Matrix(X, sparse = TRUE)
    ties <- cpp_rank_matrix_dgc(Xr@x, Xr@p, nrow(Xr), ncol(Xr))
    return(list(X_ranked = Xr, ties = ties))
}

rank_matrix.matrix <- function(X) {
    cpp_rank_matrix_dense(X)
}

#' sumGroups
#'
#' Utility function to sum over group labels
#'
#' @param X matrix
#' @param y group labels
#' @param MARGIN whether observations are rows (=2) or columns (=1)
#'
#' @return Matrix of groups by features
#' @noRd
sumGroups <- function(X, y, MARGIN = 2) {
    if (MARGIN == 2 & nrow(X) != length(y)) {
        stop('wrong dims')
    } else if (MARGIN == 1 & ncol(X) != length(y)) {
        stop('wrong dims')
    }
    UseMethod('sumGroups')
}

sumGroups.dgCMatrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_sumGroups_dgc_T(X@x,
                            X@p,
                            X@i,
                            ncol(X),
                            nrow(X),
                            as.integer(y) - 1,
                            length(unique(y)))
    } else {
        cpp_sumGroups_dgc(X@x, X@p, X@i, ncol(X), as.integer(y) - 1,
                          length(unique(y)))
    }
}

sumGroups.matrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_sumGroups_dense_T(X, as.integer(y) - 1, length(unique(y)))
    } else {
        cpp_sumGroups_dense(X, as.integer(y) - 1, length(unique(y)))
    }
}



#' nnzeroGroups
#'
#' Utility function to compute number of zeros-per-feature within group
#'
#' @param X matrix
#' @param y group labels
#' @param MARGIN whether observations are rows (=2) or columns (=1)
#'
#' @return Matrix of groups by features
#' @noRd
nnzeroGroups <- function(X, y, MARGIN = 2) {
    if (MARGIN == 2 & nrow(X) != length(y)) {
        stop('wrong dims')
    } else if (MARGIN == 1 & ncol(X) != length(y)) {
        stop('wrong dims')
    }
    UseMethod('nnzeroGroups')
}

nnzeroGroups.dgCMatrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_nnzeroGroups_dgc_T(X@p,
                               X@i,
                               ncol(X),
                               nrow(X),
                               as.integer(y) - 1,
                               length(unique(y)))
    } else {
        cpp_nnzeroGroups_dgc(X@p, X@i, ncol(X), as.integer(y) - 1,
                             length(unique(y)))
    }
}

nnzeroGroups.matrix <- function(X, y, MARGIN = 2) {
    if (MARGIN == 1) {
        cpp_nnzeroGroups_dense_T(X, as.integer(y) - 1, length(unique(y)))
    } else {
        cpp_nnzeroGroups_dense(X, as.integer(y) - 1, length(unique(y)))
    }
}


################################################################################
# Visualization ####
################################################################################

plotMarkerHeatmap <- function(
        object,
        result,
        topN = 5,
        dedupBy = c("logFC", "padj"),
        groupBy = c("dataset", "leiden_cluster"),
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
    result <- result[-which(duplicated(result$feature)), ]
    # TODO
    result <- result %>%
        dplyr::group_by(.data[["group"]]) %>%
        dplyr::top_n(topN, .data[["logFC"]]) %>%
        dplyr::arrange(.data[["group"]]) %>%
        as.data.frame()

    object <- downsample(object, maxCells = groupSize,
                         balance = groupBy,
                         useSlot = "normData")
    featureAnn <- result[, "group", drop = FALSE]
    rownames(featureAnn) <- result$feature
    plotGeneHeatmap(object, features = result$feature,
                    useCellMeta = groupBy,
                    featureAnnotation = featureAnn,
                    cellSplitBy = rev(groupBy),
                    featureSplitBy = "group",
                    showFeatureLegend = FALSE,
                    cluster_columns = FALSE,
                    cluster_column_slices = FALSE,
                    cluster_rows = FALSE,
                    cluster_row_slices = FALSE,
                    column_title = column_title,
                    ...)
}

