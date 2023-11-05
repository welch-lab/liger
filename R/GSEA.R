#' Analyze biological interpretations of metagene
#' @description Identify the biological pathways (gene sets from Reactome) that
#' each metagene (factor) might belongs to.
#' @param object A \linkS4class{liger} object with valid factorization result.
#' @param genesets Character vector of the Reactome gene sets names to be
#' tested. Default \code{NULL} uses all the gene sets from the Reactome.
#' @param useW Logical, whether to use the shared factor loadings (\eqn{W}).
#' Default \code{TRUE}.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets where the \eqn{V} matrices will be
#' included for analysis. Default \code{NULL} uses all datasets.
#' @param customGenesets A named list of character vectors of entrez gene ids.
#' Default \code{NULL} uses all the gene symbols from the input matrix.
#' @param gene_sets,mat_w,mat_v,custom_gene_sets \bold{Deprecated}. See Usage
#' section for replacement.
#' @return A list of matrices with GSEA analysis for each factor
#' @export
#' @examples
#' runGSEA(pbmcPlot)
runGSEA <- function(
        object,
        genesets = NULL,
        useW = TRUE,
        useDatasets = NULL,
        customGenesets = NULL,
        # Deprecated coding style
        gene_sets = genesets,
        mat_w = useW,
        mat_v = useDatasets,
        custom_gene_sets = customGenesets
) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
        stop("Package \"org.Hs.eg.db\" needed for this function to work. ",
             "Please install it by command:\n",
             "BiocManager::install('org.Hs.eg.db')",
             call. = FALSE)

    if (!requireNamespace("reactome.db", quietly = TRUE))
        stop("Package \"reactome.db\" needed for this function to work. ",
             "Please install it by command:\n",
             "BiocManager::install('reactome.db')",
             call. = FALSE)

    if (!requireNamespace("fgsea", quietly = TRUE))
        stop("Package \"fgsea\" needed for this function to work. ",
             "Please install it by command:\n",
             "BiocManager::install('fgsea')",
             call. = FALSE)

    .deprecateArgs(list(gene_sets = "genesets",
                        mat_w = "useW",
                        mat_v = "useDatasets",
                        custom_gene_sets = "customGenesets"))
    useDatasets <- .checkUseDatasets(object, useDatasets = useDatasets)
    .checkValidFactorResult(object, useDatasets)

    # list of V matrices: gene x k
    Vs <- getMatrix(object, "V", dataset = useDatasets, returnList = TRUE)
    # Get gene ranks in each factor
    geneLoading <- Reduce("+", Vs)
    if (isTRUE(useW)) geneLoading <- geneLoading + getMatrix(object, "W")
    geneRanks <- t(apply(geneLoading, MARGIN = 2, FUN = rank))
    colnames(geneRanks) <- toupper(colnames(geneRanks))
    # Rename gene with entrez ID
    geneID <- as.character(
        AnnotationDbi::mapIds(
            org.Hs.eg.db::org.Hs.eg.db,
            colnames(geneRanks),
            "ENTREZID",
            "SYMBOL"
        )
    )
    geneRanks <- geneRanks[, !is.na(geneID)]
    geneID <- geneID[!is.na(geneID)]
    colnames(geneRanks) <- geneID

    # Prepare `pathways` which should be a list, where each
    # element is a character vector of entrez IDs.
    if (is.null(customGenesets)) {
        pathways <- fgsea::reactomePathways(colnames(geneRanks))
        if (!is.null(genesets)) {
            pathways <- pathways[intersect(genesets, names(pathways))]
        }
    } else {
        if (inherits((customGenesets)[1], "tbl_df")) {
            pathways <- split(customGenesets,
                              x = customGenesets$entrez_gene,
                              f = customGenesets$gs_name)
            pathways <- lapply(pathways, as.character)
        } else {
            pathways <- customGenesets
        }
    }

    gsea <- apply(geneRanks, MARGIN = 1, function(x) {
        fgsea::fgsea(
            pathways,
            x,
            minSize = 15,
            maxSize = 500,
            nperm = 10000
        )
    })
    gsea <- lapply(gsea, function(x) {
        as.data.frame(x[order(x$padj), ])
    })
    names(gsea) <- paste0("Factor_", seq_along(gsea))
    return(gsea)
}

#' Run Gene Ontology enrichment analysis on metagenes
#' @description
#' This function forms genesets basing on the differential expression result,
#' and calls gene ontology (GO) analysis method provided by gprofiler2.
#' @param result Data frame of unfiltered output from \code{\link{runWilcoxon}}.
#' @param group Selection of one group available from \code{result$group}.
#' Default \code{NULL} uses all groups involved in DE \code{result} table.
#' @param useBg Logical, whether to set all genes involved in DE analysis
#' (before threshold filtering) as a domain background of GO analysis. Default
#' \code{TRUE}.
#' @param orderBy Name of DE statistics metric to order the gene list for each
#' group. Choose from \code{"logFC"} (default), \code{"pval"} or \code{"padj"}.
#' Or set \code{NULL} to turn off ranked mode.
#' @param logFCThresh The log2FC threshold above which the genes will be used.
#' Default \code{1}.
#' @param padjThresh The adjusted p-value threshold less than which the genes
#' will be used. Default \code{0.01}.
#' @param ... Additional arguments passed to \code{\link[gprofiler2]{gost}}.
#' @references Kolberg, L. et al, 2020 and Raudvere, U. et al, 2019
#' @return A list object with the following entries
#' \item{result}{data.frame of main GO analysis result.}
#' \item{meta}{Meta information for the query.}
#'
#' See \code{\link[gprofiler2]{gost}}. for detailed explanation.
#' @export
#' @examples
#' res <- runWilcoxon(pbmcPlot)
#' # Setting `significant = FALSE` because it's hard for a gene list obtained
#' # from small test dataset to represent real-life biology.
#' go <- runGOEnrich(res, group = 0, significant = FALSE)
runGOEnrich <- function(
        result,
        group = NULL,
        useBg = TRUE,
        orderBy = "logFC",
        logFCThresh = 1,
        padjThresh = 0.01,
        ...
) {
    if (!requireNamespace("gprofiler2", quietly = TRUE))
        stop("Package \"gprofiler2\" needed for this function to work. ",
             "Please install it by command:\n",
             "install.packages('gprofiler2')",
             call. = FALSE)
    if (is.null(group)) group <- unique(result$group)
    if (any(!group %in% result$group)) {
        stop("Selected groups not available `result`: ",
             paste(group[!group %in% result$group], collapse = ", "))
    }
    bg <- NULL
    domain_scope <- "annotated" # gprofiler2 default
    if (isTRUE(useBg)) {
        bg <- unique(result$feature)
        domain_scope <- "custom"
    }
    filter <- result$group %in% group &
        result$logFC > logFCThresh &
        result$padj < padjThresh
    result <- result[filter, ]

    ordered_query <- FALSE
    if (!is.null(orderBy)) {
        ordered_query <- TRUE
        if (length(orderBy) > 1) stop("Only one `orderBy` metric allowed")
        if (!orderBy %in% c("logFC", "pval", "padj")) {
            stop("`orderBy` should be one of 'logFC', 'pval' or 'padj'.")
        }
        if (orderBy == "logFC") {
            result <- result[order(result$logFC, decreasing = TRUE),]
        } else {
            result <- result[order(result[[orderBy]], decreasing = FALSE),]
        }
    }

    gsLists <- split(result$feature, droplevels(result$group))
    output <- gprofiler2::gost(
        query = gsLists, custom_bg = bg, domain_scope = domain_scope,
        ordered_query = ordered_query, ...
    )

    return(output)
}
