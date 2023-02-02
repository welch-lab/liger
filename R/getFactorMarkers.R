#' Find shared and dataset-specific markers
#'
#' Applies various filters to genes on the shared (W) and dataset-specific (V) components of the
#' factorization, before selecting those which load most significantly on each factor (in a shared
#' or dataset-specific way).
#'
#' @param object \code{liger} object. Should call optimizeALS before calling.
#' @param dataset1 Name of first dataset (default first dataset by order)
#' @param dataset2 Name of second dataset (default second dataset by order)
#' @param factorShareThresh Use only factors with a dataset specificity less than or equalt to
#'   threshold (default 10).
#' @param datasetSpecificity Pre-calculated dataset specificity if available. Will calculate if not
#'   available.
#' @param logFCThresh Lower log-fold change threshold for differential expression in markers
#'   (default 1).
#' @param pvalThresh Upper p-value threshold for Wilcoxon rank test for gene expression
#'   (default 0.05).
#' @param nGenes Max number of genes to report for each dataset (default 30).
#' @param printGenes Print ordered markers passing logfc, umi and frac thresholds (default FALSE).
#' @param verbose Print messages (TRUE by default)
#' @param factor.share.thresh,dataset.specificity,log.fc.thresh,pval.thresh,num.genes,print.genes Deprecated
#' @return List of shared and specific factors. First three elements are dataframes of dataset1-
#'   specific, shared, and dataset2-specific markers. Last two elements are tables indicating the
#'   number of factors in which marker appears.
#' @export
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
        verbose = TRUE,
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
                        pval.thresh = "pvalThresh",
                        num.genes = "nGenes",
                        print.genes = "printGenes"),
                   call = rlang::call_args(match.call()))

    if (any(isH5Liger(object, dataset = c(dataset1, dataset2))))
        stop("Please use in-memory liger object for this analysis.`")
    if (is.null(nGenes)) {
        nGenes <- length(var.features(object))
    }
    if (is.null(datasetSpecificity)) {
        datasetSpecificity <- calcDatasetSpecificity(object,
                                                     dataset1 = dataset1,
                                                     dataset2 = dataset2,
                                                     do.plot = FALSE)[[3]]
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
        labels[[i]] <- useFactors[apply(H_scaled[[i]][, useFactors],
                                        1, which.max)]
    }
    names(labels) <- names(H_scaled)

    V1_matrices <- list()
    V2_matrices <- list()
    W_matrices <- list()
    vargene <- var.features(object)
    for (j in 1:length(useFactors)) {
        i <- useFactors[j]

        W <- getMatrix(object, "W")
        V1 <- getMatrix(object, "V", dataset = dataset1)
        V2 <- getMatrix(object, "V", dataset = dataset2)

        # if not max factor for any cell in either dataset
        if (sum(labels[[dataset1]] == i) <= 1 |
            sum(labels[[dataset2]] == i) <= 1) {
            warning("Factor ", i, " did not appear as max in ",
                    "any cell in either dataset", immediate. = TRUE)
            next
        }

        # filter genes by gene_count and cell_frac thresholds
        normDataList <- getMatrix(object, "norm.data",
                                  dataset = c(dataset1, dataset2))
        normData <- cbind(normDataList[[1]][vargene, labels[[1]] == i],
                          normDataList[[2]][vargene, labels[[2]] == i])
        cellLabel <- rep(c(dataset1, dataset2),
                         c(sum(labels[[1]] == i), sum(labels[[2]] == i)))
        wilcoxResult <- wilcoxauc(log(normData + 1e-10), cellLabel)

        log2fc <- wilcoxResult[wilcoxResult$group == dataset1, ]$logFC
        names(log2fc) = wilcoxResult[wilcoxResult$group == dataset1, ]$feature
        filteredGenesV1 = wilcoxResult[wilcoxResult$logFC > logFCThresh &
                                           wilcoxResult$pval < pvalThresh,
                                       "feature"]
        filteredGenesV2 = wilcoxResult[-wilcoxResult$logFC > logFCThresh &
                                           wilcoxResult$pval < pvalThresh,
                                       "feature"]

        W <- pmin(W + V1, W + V2)
        V1 <- V1[filteredGenesV1, , drop = FALSE]
        V2 <- V2[filteredGenesV2, , drop = FALSE]

        if (length(filteredGenesV1) == 0) {
            topGenesV1 <- character(0)
        } else {
            topGeneIdx1 <- order(V1[, i], decreasing = TRUE)[1:nGenes]
            topGenesV1 <- row.names(V1)[topGeneIdx1]
            topGenesV1 <- topGenesV1[!is.na(topGenesV1)]
            topGenesV1 <- topGenesV1[V1[topGenesV1, i] > 0]
        }
        if (length(filteredGenesV2) == 0) {
            topGenesV2 <- character(0)
        } else {
            topGenesIdx2 <- order(V2[, i], decreasing = TRUE)[1:num.genes]
            topGenesV2 <- row.names(V2)[topGenesIdx2]
            topGenesV2 <- topGenesV2[!is.na(topGenesV2)]
            topGenesV2 <- topGenesV2[V2[topGenesV2, i] > 0]
        }
        topGeneIdxW <- order(W[, i], decreasing = TRUE)[1:num.genes]
        topGenesW <- row.names(W)[topGeneIdxW]
        topGenesW <- topGenesW[!is.na(topGenesW)]
        topGenesW <- topGenesW[W[topGenesW, i] > 0]

        if (isTRUE(verbose) && isTRUE(print.genes)) {
            .log("Factor ", i)
            message("Dataset 1:\n",
                 paste(topGenesV1, collapse = ", "),
                 "\nShared:\n",
                 paste(topGenesW, collapse = ", "),
                 "\nDataset 2\n",
                 paste(topGenesV2, collapse = ", "), "\n")
        }

        # order is V1, V2, W
        topGenes <- list(V1 = topGenesV1, V2 = topGenesV2, W = topGenesW)
        pvals <- lapply(topGenes, function(tg) {
            wilcoxResult[wilcoxResult$feature %in% tg &
                             wilcoxResult$group == dataset1, "pval"]
        })
        # bind values in matrices
        V1_matrices[[j]] <- Reduce(cbind, list(
            rep(i, length(topGenesV1)), topGenesV1,
            log2fc[topGenesV1], pvals$V1
        ))
        V2_matrices[[j]] <- Reduce(cbind, list(
            rep(i, length(topGenesV2)), topGenesV2,
            log2fc[topGenesV2], pvals$V2
        ))
        W_matrices[[j]] <- Reduce(cbind, list(
            rep(i, length(topGenesW)), topGenesW,
            log2fc[topGenesW], pvals$W
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
            df[which(df$p_value < pval.thresh), ]
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
#'
#' This score represents the relative magnitude of the dataset-specific components of each factor's
#' gene loadings compared to the shared components for two datasets. First, for each dataset we
#' calculate the norm of the sum of each factor's shared loadings (W) and dataset-specific loadings
#' (V). We then determine the ratio of these two values and subtract from 1... TODO: finish
#' description.
#'
#' @param object \code{liger} object. Should run optimizeALS before calling.
#' @param dataset1 Name of first dataset (by default takes first two datasets for dataset1 and 2)
#' @param dataset2 Name of second dataset
#' @param do.plot Display barplot of dataset specificity scores (by factor) (default TRUE).
#' @return List containing three elements. First two elements are the norm of each metagene factor
#' for each dataset. Last element is the vector of dataset specificity scores.
#' @export
calcDatasetSpecificity <- function(
        object,
        dataset1,
        dataset2,
        do.plot = FALSE
) {
    H1 <- getMatrix(object, slot = "H", dataset = 1)
    # V: List of two g x k matrices
    V <- getMatrix(object, slot = "V", dataset = c(dataset1, dataset2))
    W <- getMatrix(object, slot = "W")
    k <- nrow(H1)
    pct1 <- rep(0, k)
    pct2 <- rep(0, k)
    for (i in 1:k) {
        pct1[i] <- norm(as.matrix(V[[dataset1]][,i] + W[,i]), "F")
        pct2[i] <- norm(as.matrix(V[[dataset2]][,i] + W[,i]), "F")
    }
    if (do.plot) {
        graphics::barplot(
            100 * (1 - (pct1 / pct2)),
            xlab = "Factor",
            ylab = "Percent Specificity",
            main = "Dataset Specificity of Factors",
            names.arg = 1:k,
            cex.names = 0.75,
            mgp = c(2, 0.5, 0)
        ) # or possibly abs(pct1-pct2)
    }
    result <- list(pct1 = pct1,
                   pct2 = pct2,
                   pctSpec = 100 * (1 - (pct1 / pct2)))
    return(result)
}
