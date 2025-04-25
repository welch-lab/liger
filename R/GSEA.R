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
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
#'     requireNamespace("reactome.db", quietly = TRUE) &&
#'     requireNamespace("fgsea", quietly = TRUE) &&
#'     requireNamespace("AnnotationDbi", quietly = TRUE)) {
#'     runGSEA(pbmcPlot)
#' }
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
            {.code BiocManager::install('fgsea')}")

    if (!requireNamespace("AnnotationDbi", quietly = TRUE))
        cli::cli_abort(
            "Package {.pkg AnnotationDbi} is needed for this function to work.
            Please install it by command:
            {.code BiocManager::install('AnnotationDbi')}")  # nocov end

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
#' \code{TRUE}. Otherwise use all annotated genes from gprofiler2 database.
#' @param orderBy Name of DE statistics metric to order the gene list for each
#' group. Choose from \code{"logFC"}, \code{"pval"} or \code{"padj"} to enable
#' ranked mode. Default \code{NULL} to use two-list mode.
#' @param logFCThresh The absolute valued log2FC threshold above which the
#' genes will be used. Default \code{1}.
#' @param padjThresh The adjusted p-value threshold less than which the genes
#' will be used. Default \code{0.05}.
#' @param splitReg Whether to have queries of both up-regulated and
#' down-regulated genes for each group. Default \code{FALSE} only queries
#' up-regulated genes and should be preferred when \code{result} comes from
#' marker detection test. When \code{result} comes from group-to-group DE test,
#' it is recommended to set \code{splitReg = TRUE}.
#' @param ... Additional arguments passed to \code{gprofiler2::gost()}. Useful
#' ones are:
#' \describe{
#' \item{\code{organism}}{The organism to be used for the analysis. "hsapiens"
#' for human, "mmusculus" for mouse.}
#' \item{\code{evcodes}}{Whether to include overlapping genes for each term.
#' Default \code{FALSE}.}
#' \item{\code{significant}}{Whether to filter out non-significant terms.
#' Default \code{TRUE}.}
#' }
#' Arguments \code{query}, \code{custom_bg}, \code{domain_scope}, and
#' \code{ordered_query} are pre-specified by this wrapper function.
#' @references Kolberg, L. et al, 2020 and Raudvere, U. et al, 2019
#' @details
#' GO term enrichment test often goes with two modes: two-list mode and ranked
#' mode.
#'
#' Two-list mode comes with a query gene set and a background gene set.
#' A query gene set contains the filtered DEGs in this analysis. A background
#' can be all the genes involved in the DEG test (default, \code{useBg = TRUE}),
#' or use all annotated genes in the gprofiler2 database (\code{useBg = FALSE}).
#'
#' Ranked mode comes with only one query gene set, which is sorted. It should
#' contain the whole domain background genes while significant genes are
#' supposed to come first. Set \code{orderBy} to one of the DE statistics metric
#' to enable this mode. \code{useBg} will be ignored in this mode.
#'
#' @return A list object where each element is a result list for a group. Each
#' result list contains two elements:
#' \item{result}{data.frame of main GO analysis result.}
#' \item{meta}{Meta information for the query.}
#'
#' See \code{gprofiler2::gost()}. for detailed explanation.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("gprofiler2", quietly = TRUE)) {
#'     go <- runGOEnrich(deg.pw)
#' }
#' }
runGOEnrich <- function(
        result,
        group = NULL,
        useBg = TRUE,
        orderBy = NULL,
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

    group <- group %||% as.character(unique(result$group))
    group <- rlang::arg_match(
        arg = group,
        values =  as.character(unique(result$group)),
        multiple = TRUE
    )
    resultList <- list()
    for (g in group) {
        groupResult <- result %>% dplyr::filter(.data[['group']] == g)
        query <- list()
        if (!is.null(orderBy)) {
            # Ranked mode
            orderBy <- rlang::arg_match(arg = orderBy,
                                        values = c("logFC", "pval", "padj"),
                                        multiple = FALSE)
            if (orderBy == 'logFC') {
                query$Up <- groupResult %>%
                    dplyr::arrange(dplyr::desc(.data[['logFC']])) %>%
                    dplyr::pull(.data[['feature']])
                if (isTRUE(splitReg)) {
                    query$Down <- groupResult %>%
                        dplyr::arrange(.data[['logFC']]) %>%
                        dplyr::pull(.data[['feature']])
                }
            } else {
                if (isTRUE(splitReg)) {
                    cli::cli_alert_warning(
                        "Unable to split for up- and down-regulated genes when ranking by {orderBy}."
                    )
                }
                query <- groupResult %>%
                    dplyr::arrange(.data[[orderBy]]) %>%
                    dplyr::pull(.data[['feature']])
            }
            res <- gprofiler2::gost(
                query = query,
                ordered_query = TRUE,
                ...
            )
            res$result <- res$result %>%
                dplyr::mutate(
                    fold_enrichment =
                        (.data[['intersection_size']] / .data[['term_size']]) /
                        (.data[['query_size']] / .data[['effective_domain_size']])
                )
        } else {
            # Two-list mode
            query$Up <- groupResult %>%
                dplyr::filter(.data[['logFC']] > logFCThresh,
                              .data[['padj']] < padjThresh) %>%
                dplyr::pull(.data[['feature']])
            if (isTRUE(splitReg)) {
                query$Down <- groupResult %>%
                    dplyr::filter(.data[['logFC']] < -logFCThresh,
                                  .data[['padj']] < padjThresh) %>%
                    dplyr::pull(.data[['feature']])
            }
            bg <- NULL
            domainScope <- 'annotated'
            if (isTRUE(useBg)) {
                bg <- unique(groupResult$feature)
                domainScope <- 'custom'
            }
            res <- gprofiler2::gost(
                query = query,
                custom_bg = bg,
                domain_scope = domainScope,
                ...
            )
            res$result <- res$result %>%
                dplyr::mutate(
                    fold_enrichment =
                        (.data[['intersection_size']] / .data[['term_size']]) /
                        (.data[['query_size']] / .data[['effective_domain_size']])
                )
        }
        resultList[[g]] <- res
    }
    return(resultList)
}


#' Visualize GO enrichment test result in dot plot
#' @param result Returned list object from \code{\link{runGOEnrich}}.
#' @param group Character vector of group names, must be available in
#' \code{names(result)}. Default \code{NULL} make plots for all groups.
#' @param query A single string selecting from which query to show the result.
#' Choose from \code{"Up"} for results using up-regulated genes, \code{"Down"}
#' for down-regulated genes. Default NULL use anything available.
#' @param pvalThresh Numeric scalar, cutoff for p-value where smaller values are
#' considered as significant. Default \code{0.05}.
#' @param n Number of top terms to be shown, ranked by p-value. Default
#' \code{20}.
#' @param minDotSize The size of the dot representing the minimum gene count.
#' Default \code{3}.
#' @param maxDotSize The size of the dot representing the maximum gene count.
#' @param termIDMatch Regular expression pattern to match the term ID. Default
#' \code{"^GO"} for only using GO terms from returned results.
#' @param colorPalette,colorDirection Viridis palette options. Default
#' \code{"E"} and \code{1}.
#' @inheritDotParams .ggplotLigerTheme title subtitle legendColorTitle legendSizeTitle showLegend legendPosition baseSize titleSize subtitleSize xTextSize xTitleSize yTextSize yTitleSize legendTextSize legendTitleSize plotly
#' @return A ggplot object if only one group or a list of ggplot objects.
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("gprofiler2", quietly = TRUE)) {
#'    go <- runGOEnrich(deg.pw)
#'    plotGODot(go)
#' }
#' }
plotGODot <- function(
        result,
        group = NULL,
        query = NULL,
        pvalThresh = 0.05,
        n = 20,
        minDotSize = 3,
        maxDotSize = 7,
        termIDMatch = "^GO",
        colorPalette = "E",
        colorDirection = -1,
        ...
) {
    group <- group %||% names(result)
    group <- rlang::arg_match(
        arg = group,
        values = names(result),
        multiple = TRUE
    )
    query <- query %||% as.character(unique(result[[group]]$result$query))
    query <- rlang::arg_match(
        arg = query,
        values = as.character(unique(result[[group]]$result$query)),
        multiple = FALSE
    )
    plotList <- list()
    for (i in seq_along(group)) {
        gname <- group[i]
        resdf <- result[[gname]]$result
        # Can't use dplyr call here because the environment object name 'query'
        # will be wrongly understood as a column.
        resdf <- resdf[resdf$query == query, , drop = FALSE]
        if (is.null(resdf) || nrow(resdf) == 0) {
            cli::cli_alert_warning(
                "No result returned for group {.val {gname}} and query {.val {query}}."
            )
            next
        }
        resdf %<>% dplyr::filter(
            grepl(termIDMatch, .data[['term_id']]),
            .data[['p_value']] <= pvalThresh
        )
        if (nrow(resdf) == 0) {
            cli::cli_alert_warning(
                "No enough matching terms ({.val termIDMatch}) nor significant terms (p-value <= {.val {pvalThresh}}) for group {.val {gname}}."
            )
            next
        }
        g <- resdf %>%
            dplyr::select(dplyr::all_of(c(
                'term_name',
                'p_value',
                'intersection_size',
                'fold_enrichment'
            ))) %>%
            dplyr::arrange(.data[['p_value']]) %>%
            dplyr::slice_head(n = n) %>%
            dplyr::mutate(
                term_name = factor(
                    .data[['term_name']],
                    levels = rev(.data[['term_name']])
                )
            ) %>%
            ggplot2::ggplot(ggplot2::aes(
                x = .data[['fold_enrichment']],
                y = .data[['term_name']],
                size = .data[['intersection_size']],
                color = -log10(.data[['p_value']])
            )) +
            ggplot2::geom_segment(
                mapping = ggplot2::aes(
                    x = 0,
                    y = .data[['term_name']],
                    xend = .data[['fold_enrichment']],
                    yend = .data[['term_name']],
                    color = -log10(.data[['p_value']])
                ),
                size = 0.5
            ) +
            ggplot2::geom_point() +
            ggplot2::scale_size_continuous(range = c(minDotSize, maxDotSize)) +
            ggplot2::labs(x = 'Fold Enrichment',
                          y = NULL,
                          color = bquote(~-log[10](p-value)),
                          size = "Gene Count")
        plotList[[gname]] <- .ggplotLigerTheme(
            plot = g,
            colorPalette = colorPalette,
            colorDirection = colorDirection,
            panelBorder = TRUE,
            ...
        ) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
            ggplot2::theme(
                legend.title.position = 'left',
                legend.title = ggplot2::element_text(angle = 270, hjust = 0)
            )
    }
    if (length(plotList) == 1) {
        return(plotList[[1]])
    } else {
        return(plotList)
    }
}
