#' Test all factors for enrichment in a gene set
#' @description
#' This function takes the factorized \eqn{W} matrix, with gene loading in
#' factors, to get the ranked gene list for each factor. Then it runs simply
#' implemented GSEA against given gene sets. So if genes in the given gene set
#' are top loaded in a factor, this function will return high positive
#' enrichment score (ES) as well as significant p-value.
#'
#' For the returned result object, use \code{print()} or \code{summary()} to
#' show concise results, and use \code{plot()} to visualize the GSEA statistics.
#'
#' This function can be useful in various scenarios:
#'
#' For example, when clusters with strong cell cycle activity are detected,
#' users can apply this function with cell cycle gene sets to identify if any
#' factor is enriched with such genes. Then in the downstream when aligning the
#' iNMF factor loadings, users can simply opt to exclude these factors so the
#' variation in cell cycle is regressed out. Objects \code{cc.gene.human} and
#' \code{cc.gene.mouse} are deliverered in package for convenience.
#'
#' In other cases, this function can also be used to understand the biological
#' meaning of each cluster. Since the downstream clustering result is largely
#' determined by the top loaded factor in each cell, understanding what
#' genes are loaded in the top factor helps understand the identity and activity
#' of the cell. This will require users to have there own gene sets prepared.
#'
#' @param object A \linkS4class{liger} object with factorized \eqn{W} matrix
#' available.
#' @param geneSet A character vector for a single gene set, or a list of
#' character vectors for multiple gene sets.
#' @param nPerm Integer number for number of permutations to estimate p-value.
#' Default \code{1000}.
#' @param seed Integer number for random seed. Default \code{1}. Set to
#' \code{NULL} to not set seed.
#' @param verbose Logical, whether to print progress bar. Default
#' \code{getOptions('ligerVerbose')} otherwise \code{TRUE}.
#' @return If \code{geneSet} is a single character vector, returns a data frame
#' with enrichment score (ES), normalized enrichment score (NES), and p-value
#' for the test in each factor. If \code{geneSet} is a list, returns a list of
#' such data frames.
#' @export
#' @examples
#' \dontrun{
#' pbmc <- pbmc %>%
#'     selectBatchHVG() %>%
#'     scaleNotCenter() %>%
#'     runINMF()
#' factorGSEAres <- factorGSEA(pbmc, ccGeneHuman)
#' # Print summary of significant results
#' print(factorGSEAres)
#' summary(factorGSEAres)
#' # Make GSEA plot for certain gene set and factor
#' plot(factorGSEAres, 'g2m.genes', 'Factor_1')
#' }
#'
factorGSEA <- function(
        object,
        geneSet,
        nPerm = 1000,
        seed = 1,
        verbose = getOption('ligerVerbose', TRUE)
) {
    # Input checks
    W <- getMatrix(object, 'W')
    if (is.null(W)) {
        cli::cli_abort('{.field W} matrix is not available. Please run {.fn runIntegration} first.')
    }

    if (is.character(geneSet)) {
        # Make the structure uniform
        geneSet <- list(geneset = geneSet)
    } else if (is.list(geneSet)) {
        for (gs in geneSet) {
            if (!is.character(gs)) {
                cli::cli_abort('Elements in {.field geneSet} must all be character vectors.')
            }
        }
    } else {
        cli::cli_abort('{.field geneSet} must be either a character vector or a list of character vectors.')
    }

    # Start processing
    # Get the list of genes for each factor,
    # ranked by the gene loading in each factor
    rankedGenes <- apply(W, MARGIN = 2, FUN = function(x) {
        names(sort(x, decreasing = TRUE))
    })
    results <- list()
    cliID <- NULL
    if (isTRUE(verbose)) {
        cliID <- cli::cli_progress_bar(
            name = 'Running GSEA on factors',
            total = length(geneSet) * nPerm * ncol(rankedGenes)
        )
    }
    if (!is.null(seed)) set.seed(seed)
    for (i in seq_along(geneSet)) {
        gsname <- names(geneSet)[i]
        if (is.null(gsname) || nchar(gsname) == 0) {
            gsname <- paste0('geneset', i)
        }
        res <- list()
        for (j in seq_len(ncol(rankedGenes))) {
            res[[j]] <- gseaSingleGeneSet(
                rankedGenes = rankedGenes[, j],
                gs = geneSet[[i]],
                nPerm = nPerm,
                cliID = cliID
            )
            res[[j]]$loading <- sort(W[, j], decreasing = TRUE)
        }
        names(res) <- colnames(rankedGenes)
        # res <- Reduce(rbind, res)
        # rownames(res) <- colnames(rankedGenes)
        # res$sig <- ifelse(res$pval < 0.05 & res$NES > 0, '*', '')
        # res$sig <- ifelse(res$pval < 0.01 & res$NES > 0, '**', res$sig)
        results[[gsname]] <- res
    }
    class(results) <- 'factorGSEA'
    return(results)
}

gseaSingleGeneSet <- function(rankedGenes, gs, nPerm = 1000, cliID) {
    N <- length(rankedGenes)
    inSet <- rankedGenes %in% gs
    Nh <- sum(inSet)
    if (Nh == 0) {
        return(data.frame(ES = NA, NES = NA, pval = NA))
    }

    # No weight, uniform increment for hits
    hits <- as.numeric(inSet) / Nh
    misses <- as.numeric(!inSet) / (N - Nh)
    runningSum <- cumsum(hits - misses)
    ES <- max(runningSum)
    ES_sign <- ifelse(which.max(runningSum) <= which.min(runningSum), 1, -1)
    ES <- ES * ES_sign

    # Permutations for p-value estimation
    permutedES <- numeric(nPerm)
    for (i in seq_len(nPerm)) {
        permutedInSet <- sample(inSet)
        hitsPerm <- as.numeric(permutedInSet) / Nh
        missesPerm <- as.numeric(!permutedInSet) / (N - Nh)
        runningSumPerm <- cumsum(hitsPerm - missesPerm)
        permutedES[i] <- max(runningSumPerm) * ES_sign
        if (!is.null(cliID)) {
            cli::cli_progress_update(inc = 1, id = cliID)
        }
    }

    pval <- mean(abs(permutedES) >= abs(ES))
    NES <- ES / mean(abs(permutedES))
    list(
        runningSum = runningSum,
        hits = inSet,
        ES = ES,
        NES = NES,
        pval = pval
    )
}

#' Show significant results from factorGSEA
#' @param object A \code{factorGSEA} object.
#' @param ... S3 method convention, not used for now.
#' @return A data frame of significant tests with gene set names, factor names
#' and other GSEA statistics.
#' @method summary factorGSEA
#' @export
summary.factorGSEA <- function(object, ...) {
    rows <- list()
    for (i in seq_along(object)) {
        # at gene set level
        for (j in seq_along(object[[i]])) {
            # at factor level
            if (object[[i]][[j]]$pval < 0.05 && object[[i]][[j]]$NES > 0) {
                rows <- c(rows, list(
                    data.frame(
                        geneSet = names(object)[i],
                        factor = names(object[[i]])[j],
                        ES = object[[i]][[j]]$ES,
                        NES = object[[i]][[j]]$NES,
                        pval = object[[i]][[j]]$pval
                    )
                ))
            }
        }
    }
    Reduce(rbind, rows)
}

#' Show information about factorGSEA object
#' @param x A \code{factorGSEA} object.
#' @param ... S3 method convention, not used for now.
#' @method print factorGSEA
#' @export
print.factorGSEA <- function(x, ...) {
    cat(cli::format_inline(
        'A factorGSEA object tested {length(x)} gene set{?s} over {length(x[[1]])} factor{?s}\n\n'
    ))
    cat('Significat test results:\n')
    print(summary(x))
    return(invisible(NULL))
}

#' GSEA plot for specific gene set and factor using factorGSEA results
#' @param x A \code{factorGSEA} object.
#' @param y Not used, for S3 method convention.
#' @param geneSetName A character string for the gene set name to plot.
#' @param useFactor A character string (e.g. 'Factor_1') or just numeric index
#' for the factor name to plot.
#' @param xTitleSize,yTitleSize Numeric, size for x or y axis titles,
#' respectively. Default \code{10}.
#' @param xTextSize,yTextSize Numeric, size for x or y axis text,
#' respectively. Default \code{8}.
#' @param titleSize Numeric, size for the main plot title. Default \code{12}.
#' @param captionTextSize Numeric, size for the caption text. Default \code{8}.
#' @param ESLineColor Color for the enrichment score line. Default
#' \code{'green'}.
#' @param ESLinewidth Numeric, line width for the enrichment score line.
#' Default \code{1}.
#' @param hitsLineColor Color for the hits line. Default \code{'black'}.
#' @param hitsLinewidth Numeric, line width for the hits line. Default
#' \code{0.5}.
#' @param loadingBarColor Color for the loading bar. Default \code{'grey'}.
#' @param ... Not used.
#' @return ggplot object
#' @method plot factorGSEA
#' @export
plot.factorGSEA <- function(
        x,
        y,
        geneSetName,
        useFactor,
        xTitleSize = 10,
        xTextSize = 8,
        yTitleSize = 10,
        yTextSize = 8,
        titleSize = 12,
        captionTextSize = 8,
        ESLineColor = 'green',
        ESLinewidth = 1,
        hitsLineColor = 'black',
        hitsLinewidth = 0.5,
        loadingBarColor = 'grey',
        ...
) {
    if (length(geneSetName) != 1) {
        cli::cli_abort('Only one gene set is allowed at a time while {length(geneSetName)} were given to {.field geneSetName}.')
    }
    geneSetName <- rlang::arg_match(
        arg = geneSetName, values = names(x)
    )

    if (length(useFactor) != 1) {
        cli::cli_abort('Only one factor is allowed at a time while {length(useFactor)} were given to {.field useFactor}.')
    }
    if (is.character(useFactor)) {
        useFactor <- rlang::arg_match(arg = useFactor, values = names(x[[geneSetName]]))
    } else if (is.numeric(useFactor) || is.integer(useFactor)) {
        useFactor <- names(x[[geneSetName]])[useFactor]
        if (is.na(useFactor)) {
            cli::cli_abort('{.field useFactor} is out of range. Only {length(x[[geneSetName]])} factors are available.')
        }
    }

    element <- x[[geneSetName]][[useFactor]]
    plotDF <- data.frame(
        rank = seq_along(element$runningSum),
        ES = element$runningSum,
        hits = element$hits,
        loading = element$loading
    )
    p1 <- ggplot2::ggplot(
        plotDF,
        ggplot2::aes(
            x = .data[['rank']],
            y = .data[['ES']]
        )) +
        ggplot2::geom_line(color = ESLineColor, linewidth = ESLinewidth) +
        ggplot2::theme_linedraw() +
        ggplot2::labs(
            title = sprintf('GSEA test of gene set %s in factor %s',
                            geneSetName, useFactor),
            y = 'Enrichment Score'
        ) +
        ggplot2::coord_cartesian(expand = FALSE) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = xTitleSize, face = 'bold', hjust = 0.5),
            axis.title.y = ggplot2::element_text(size = yTitleSize),
            axis.text.y = ggplot2::element_text(size = yTextSize),
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(t = 10, l = 10, r = 10, b = 0, unit = 'pt')
        )
    p2 <- ggplot2::ggplot(plotDF[plotDF$hits,]) +
        ggplot2::geom_vline(
            mapping = ggplot2::aes(xintercept = .data[['rank']]),
            color = hitsLineColor,
            linewidth = hitsLinewidth
        ) +
        ggplot2::theme_linedraw() +
        ggplot2::scale_x_continuous(limits = c(0, nrow(plotDF))) +
        ggplot2::coord_cartesian(expand = FALSE) +
        ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(t = 0, l = 10, r = 10, b = 0, unit = 'pt')
        )
    p3 <- ggplot2::ggplot(
        plotDF,
        ggplot2::aes(
            x = .data[['rank']],
            y = log1p(.data[['loading']])
        )
    ) +
        ggplot2::geom_bar(stat = 'identity', fill = loadingBarColor) +
        ggplot2::theme_linedraw() +
        ggplot2::labs(
            x = 'Rank',
            y = 'log(Gene loading + 1)',
            caption = sprintf(
                '# Hitting genes: %d\nNormalized enrichment score: %.3f\nP-value: %.3f',
                sum(plotDF$hits),
                element$NES,
                element$pval
            )
        ) +
        ggplot2::coord_cartesian(expand = FALSE) +
        ggplot2::theme(
            plot.margin = ggplot2::margin(t = 0, l = 10, r = 10, b = 10, unit = 'pt'),
            axis.title.x = ggplot2::element_text(size = xTitleSize),
            axis.text.x = ggplot2::element_text(size = xTextSize),
            axis.title.y = ggplot2::element_text(size = yTitleSize),
            axis.text.y = ggplot2::element_text(size = yTextSize),
            plot.caption = ggplot2::element_text(size = captionTextSize)

        )
    cowplot::plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(3, 1, 2), align = 'v', axis = 'lr')
}

