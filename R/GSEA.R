#' Analyze biological interpretations of metagene
#' @description Identify the biological pathways (gene sets from Reactome) that
#' each metagene (factor) might belongs to.
#' @param object A \linkS4class{liger} object with valid factorization result.
#' @param genesets Character vector of the Reactome gene sets names to be
#' tested. Default \code{NULL} uses all the gene sets from the Reactome.
#' @param useW Logical, whether to use the shared factor loadings (\eqn{W}).
#' Default \code{TRUE}.
#' @param useV A character vector of the names, a numeric or logical
#' vector of the index of the datasets where the \eqn{V} matrices will be
#' included for analysis. Default \code{NULL} uses all datasets.
#' @param customGenesets A named list of character vectors of entrez gene ids.
#' Default \code{NULL} uses all the gene symbols from the input matrix.
#' @param gene_sets,mat_w,mat_v,custom_gene_sets \bold{Deprecated}. See Usage
#' section for replacement.
#' @return A list of matrices with GSEA analysis for each factor
#' @export
#' @examples
#' \donttest{
#' runGSEA(pbmcPlot)
#' }
runGSEA <- function(
        object,
        genesets = NULL,
        useW = TRUE,
        useV = NULL,
        customGenesets = NULL,
        # Deprecated coding style
        gene_sets = genesets,
        mat_w = useW,
        mat_v = useV,
        custom_gene_sets = customGenesets
) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) # nocov start
        cli::cli_abort(
            "Package {.pkg org.Hs.eg.db} is needed for this function to work.
            Please install it by command:
            {.code BiocManager::install('org.Hs.eg.db')}")

    if (!requireNamespace("reactome.db", quietly = TRUE))
        cli::cli_abort(
            "Package {.pkg reactome.db} is needed for this function to work.
            Please install it by command:
            {.code BiocManager::install('reactome.db')}")

    if (!requireNamespace("fgsea", quietly = TRUE))
        cli::cli_abort(
            "Package {.pkg fgsea} is needed for this function to work.
            Please install it by command:
            {.code BiocManager::install('fgsea')}") # nocov end

    .deprecateArgs(list(gene_sets = "genesets",
                        mat_w = "useW",
                        mat_v = "useV",
                        custom_gene_sets = "customGenesets"))
    useV <- .checkUseDatasets(object, useDatasets = useV)
    .checkValidFactorResult(object, useV)

    # list of V matrices: gene x k
    Vs <- getMatrix(object, "V", dataset = useV, returnList = TRUE)
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

#' Run Gene Ontology enrichment analysis on differentially expressed genes.
#' @description
#' This function forms genesets basing on the differential expression result,
#' and calls gene ontology (GO) analysis method provided by gprofiler2.
#' @param result Data frame of unfiltered output from \code{\link{runMarkerDEG}}
#' or \code{\link{runPairwiseDEG}}.
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
#' will be used. Default \code{0.05}.
#' @param splitReg Whether to have queries of both up-regulated and
#' down-regulated genes for each group. Default \code{FALSE} only queries
#' up-regulated genes and should be preferred when \code{result} comes from
#' marker detection test. When \code{result} comes from group-to-group DE test,
#' it is recommended to set \code{splitReg = TRUE}.
#' @param ... Additional arguments passed to
#' \code{gprofiler2::\link[gprofiler2]{gost}}. Arguments \code{query},
#' \code{custom_bg}, \code{domain_scope}, and \code{ordered_query} are
#' pre-specified by this wrapper function. Users must set
#' \code{organism = "mmusculus"} when working on mouse data.
#' @references Kolberg, L. et al, 2020 and Raudvere, U. et al, 2019
#' @return A list object where each element is a result list for a group. Each
#' result list contains two elements:
#' \item{result}{data.frame of main GO analysis result.}
#' \item{meta}{Meta information for the query.}
#'
#' See \code{\link[gprofiler2]{gost}}. for detailed explanation.
#' @export
#' @examples
#' res <- runMarkerDEG(pbmcPlot)
#' # Setting `significant = FALSE` because it's hard for a gene list obtained
#' # from small test dataset to represent real-life biology.
#' \donttest{
#' go <- runGOEnrich(res, group = 0, significant = FALSE)
#' }
runGOEnrich <- function(
        result,
        group = NULL,
        useBg = TRUE,
        orderBy = "padj",
        logFCThresh = 1,
        padjThresh = 0.05,
        splitReg = FALSE,
        ...
) {
    if (!requireNamespace("gprofiler2", quietly = TRUE)) # nocov start
        cli::cli_abort(
            "Package {.pkg gprofiler2} is needed for this function to work.
            Please install it by command:
            {.code install.packages('gprofiler2')}") # nocov end

    group <- group %||% unique(result$group)
    if (any(!group %in% result$group)) {
        cli::cli_abort("Selected groups not available in {.code result$group}: {.val {group[!group %in% result$group]}}")
    }
    bg <- NULL
    domain_scope <- "annotated" # gprofiler2 default
    if (isTRUE(useBg)) {
        bg <- unique(result$feature)
        domain_scope <- "custom"
    }
    filter <- result$group %in% group &
        abs(result$logFC) > logFCThresh &
        result$padj < padjThresh
    filter[is.na(filter)] <- FALSE
    result <- result[filter, , drop = FALSE]
    resultUp <- result[result$logFC > 0,]
    if (isTRUE(splitReg)) resultDown <- result[result$logFC < 0,]

    ordered_query <- FALSE
    if (!is.null(orderBy)) {
        ordered_query <- TRUE
        if (length(orderBy) > 1) cli::cli_abort("Only one {.code orderBy} metric allowed")
        if (!orderBy %in% c("logFC", "pval", "padj")) {
            cli::cli_abort("{.code orderBy} should be one of {.val logFC}, {.val pval} or {.val padj}.")
        }
        if (orderBy == "logFC") {
            resultUp <- resultUp[order(resultUp$logFC, decreasing = TRUE),]
            if (isTRUE(splitReg))
                resultDown <- resultDown[order(resultDown$logFC),]
        } else {
            resultUp <- resultUp[order(resultUp[[orderBy]]),]
            if (isTRUE(splitReg))
                resultDown <- resultDown[order(resultDown[[orderBy]]),]
        }
    }
    resultList <- list()
    for (g in unique(result$group)) {
        query <- list(Up = resultUp$feature[resultUp$group == g])
        if (splitReg) query$Down <- resultDown$feature[resultDown$group == g]
        resultList[[g]] <- gprofiler2::gost(
            query = query, custom_bg = bg, domain_scope = domain_scope,
            ordered_query = ordered_query, ...
        )
    }

    return(resultList)
}

