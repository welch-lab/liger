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
#' This function forms genesets basing on the non-zero gene loading in each
#' factor, and calls gene ontology (GO) analysis method provided by gprofiler2.
#' @param object A \linkS4class{liger} object with valid factorization result.
#' @param useW Logical, whether to consider shared gene loading value (i.e.
#' \eqn{W} matrix) when forming the genesets. Default \code{TRUE}.
#' @param useDatasets A character vector of the names, a numeric or logical
#' vector of the index of the datasets where the gene loading need to be
#' considered. Default \code{NULL} considers all datasets.
#' @param sumLoading Logical, whether to sum up the gene loading from all
#' considered datasets for each factor, or to query a geneset per factor and
#' dataset. Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link[gprofiler2]{gost}}.
#' @references Kolberg, L. et al, 2020 and Raudvere, U. et al, 2019
#' @return A list object with the following entries
#' \item{result}{data.frame of main GO analysis result.}
#' \item{meta}{Meta information for the query.}
#' @export
#' @examples
#' go <- runFactorGeneGO(pbmcPlot)
#' head(go$result)
runFactorGeneGO <- function(
        object,
        useW = TRUE,
        useDatasets = NULL,
        sumLoading = TRUE,
        ...
) {
    if (!requireNamespace("gprofiler2", quietly = TRUE))
        stop("Package \"gprofiler2\" needed for this function to work. ",
             "Please install it by command:\n",
             "install.packages('gprofiler2')",
             call. = FALSE)
    useDatasets <- .checkUseDatasets(object, useDatasets = useDatasets)
    .checkValidFactorResult(object, useDatasets)

    # list of V matrices: gene x k
    geneLoading <- getMatrix(object, "V", dataset = useDatasets,
                             returnList = TRUE)
    k <- ncol(geneLoading[[1]])
    # Get gene ranks in each factor
    if (isTRUE(sumLoading)) {
        geneLoading <- list(Reduce("+", geneLoading))
    }

    if (isTRUE(useW)) {
        geneLoading <- lapply(geneLoading, function(v) {
            v + getMatrix(object, "W")
        })
    }
    genes <- rownames(geneLoading)
    gsLists <- lapply(geneLoading, function(v) {
        genes <- rownames(v)
        gs <- lapply(seq_len(ncol(v)), function(j) {
            genes <- genes[order(v[, j], decreasing = TRUE)]
            genes[v[genes, j] > 0]
        })
        names(gs) <- colnames(v)
        return(gs)
    })
    if (is.null(names(gsLists))) {
        # Summing up loadings from all used datasets
        gsLists <- gsLists[[1]]
    } else {
        prefix <- rep(names(gsLists), each = k)
        gsLists <- Reduce(c, gsLists)
        names(gsLists) <- paste0(prefix, "_", names(gsLists))
    }
    output <- gprofiler2::gost(query = gsLists, ...)
    output$meta$query_metadata$useW <- useW
    output$meta$query_metadata$sumLoading <- sumLoading
    return(output)
}

