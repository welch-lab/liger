#' Analyze biological interpretations of metagene
#'
#' Identify the biological pathways (gene sets from Reactome) that each metagene (factor) might belongs to.
#' @param object \code{liger} object.
#' @param gene_sets A list of the Reactome gene sets names to be tested. If not specified,
#' this function will use all the gene sets from the Reactome by default
#' @param mat_w This indicates whether to use the shared factor loadings 'W' (default TRUE)
#' @param mat_v This indicates which V matrix to be added to the analysis. It can be a numeric number or a list
#' of the numerics.
#' @param custom_gene_sets A named list of character vectors of entrez gene ids. If not specified,
#' this function will use all the gene symbols from the input matrix by default
#' @return A list of matrices with GSEA analysis for each factor
#' @export
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
    Vs <- getMatrix(object, "V", dataset = useDatasets)
    if (length(useDatasets) == 1) {
        # getMatrix directly returns the matrix when only one dataset
        Vs <- list(Vs)
        names(Vs) <- useDatasets
    }
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

    # TODO: bugs for using certain pathways, still don't know why
    # TODO: Try other GSEA packages
    gsea <- apply(geneRanks, MARGIN = 1, function(x) {
        fgsea::fgsea(
            pathways,
            x,
            minSize = 15,
            maxSize = 500,
            nperm = 10000
        )
        #scgsea
    })
    return(gsea)
    gsea <- lapply(gsea, function(x) {
        as.matrix(x[order(x$pval), ])
    })
    return(gsea)
}
