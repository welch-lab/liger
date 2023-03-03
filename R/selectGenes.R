#' Select a subset of informative genes
#' @description This function identifies highly variable genes from each dataset
#' and combines these gene sets (either by union or intersection) for use in
#' downstream analysis. Assuming that gene expression approximately follows a
#' Poisson distribution, this function identifies genes with gene expression
#' variance above a given variance threshold (relative to mean gene expression).
#' @param object \linkS4class{liger} object. Normalization has to be performed
#' in advance. See \code{\link{normalize}}.
#' @param var.thresh Variance threshold. Main threshold used to identify
#' variable genes. Genes with expression variance greater than threshold
#' (relative to mean) are selected. Higher threshold results in fewer selected
#' genes. Accepts single value or vector with specific threshold for each
#' dataset in \code{datasets.use}. Default \code{0.1}.
#' @param alpha.thresh Alpha threshold. Controls upper bound for expected mean
#' gene expression. Lower threshold means higher upper bound. Default
#' \code{0.99}.
#' @param num.genes Number of genes to find for each dataset. Optimizes the
#' value of \code{var.thresh} for each dataset to get this number of genes.
#' Accepts single value or vector with same length as number of
#' \code{datasets.use}. Default \code{NULL} does not optimize.
#' @param tol Tolerance to use for optimization if \code{num.genes} is
#' specified. Default \code{0.0001}.
#' @param capitalize Capitalize gene names to match homologous genes (i.e.
#' across species). Default \code{FALSE}.
#' @param combine How to combine variable genes across datasets. Choose from
#' \code{"union"} or \code{"intersection"}. Default \code{"union"}.
#' @param datasets.use A character vector of the names, a numeric or logical
#' vector of the index of the datasets to be included for selection. Default
#' \code{NULL} considers all datasets.
#' @param unshared Logical. Whether to consider unshared features. Default
#' \code{FALSE}.
#' @param unshared.datasets The dataset to be included for considering unshared
#' features. Similar to how \code{datasets.use} should be specified. Default
#' \code{NULL}.
#' @param unshared.thresh The threshold to be applied for each dataset in
#' \code{unshared.datasets}. Similar to how \code{var.thresh} should be
#' specified. Default \code{NULL}.
#' @param chunkSize Integer. Number of maximum number of cells in each chunk,
#' when gene selection is applied to any HDF5 based dataset. Default
#' \code{1000}.
#' @param verbose Logical. Whether to show the progress. Default \code{TRUE}.
#' @return The input \code{object}, with the \code{var.features} slot updated.
#' The \code{var.unshared.features} slot of the \linkS4class{ligerDataset}
#' object in the \code{datasets} slot will be updated if \code{unshared = TRUE}.
#' @seealso \code{\link{plotVarFeatures}}
#' @useDynLib rliger, .registration = TRUE
#' @export
selectGenes <- function(
        object,
        var.thresh = 0.1,
        alpha.thresh = 0.99,
        num.genes = NULL,
        tol = 0.0001,
        capitalize = FALSE,
        combine = c("union", "intersection"),
        datasets.use = NULL,
        unshared = FALSE,
        unshared.datasets = NULL,
        unshared.thresh = NULL,
        chunkSize = 1000,
        verbose = TRUE
) {
    # A bunch of input checks at first ####
    combine <- match.arg(combine)
    datasets.use <- .checkUseDatasets(object, datasets.use)
    object <- recordCommand(object, dependencies = "hdf5r")

    if (length(var.thresh) != 1 & length(var.thresh) != length(datasets.use))
        stop("Wrong length of `var.thresh`. Use 1 or `length(object)` values.")
    if (length(var.thresh) == 1)
        var.thresh <- rep(var.thresh, length(datasets.use))
    datasets.involved <- datasets.use

    ## Checks for the same thing for unshared variable features ####
    if (isTRUE(unshared)) {
        unshared.datasets <- .checkUseDatasets(object, unshared.datasets)
        if (length(unshared.thresh) != 1 &
            length(unshared.thresh) != length(unshared.datasets))
            stop("Wrong length of `unshared.thresh`. ",
                 "Use 1 or `length(object)` values.")
        if (length(unshared.thresh) == 1)
            unshared.thresh <- rep(unshared.thresh, length(unshared.datasets))

        datasets.involved <- unique(c(datasets.involved, unshared.datasets))
    } else unshared.datasets <- NULL

    # Identify all variable features for datasets being involved ####
    shared.features <- Reduce(intersect, lapply(datasets(object), rownames))
    selected.shared <- list()
    selected.unshared <- list()
    for (d in datasets.involved) {
        i <- datasets.use == d
        j <- unshared.datasets == d
        if (isTRUE(verbose))
            .log("Selecting HVG for dataset: ", d)
        ld <- dataset(object, d)
        ## Make sure that all required feature meta values exist ####
        if (isH5Liger(ld)) {
            ld <- calcGeneVars.H5(ld, chunkSize = chunkSize,
                                  verbose = verbose)
        } else {
            feature.meta(ld, check = FALSE)$geneMeans <-
                rowMeansFast(norm.data(ld))
            feature.meta(ld, check = FALSE)$geneVars <-
                rowVarsFast(norm.data(ld), feature.meta(ld)$geneMeans)
        }
        datasets(object, check = FALSE)[[d]] <- ld
        ## The real calculation starts here ####
        geneMeans <- feature.meta(ld)$geneMeans
        geneVars <- feature.meta(ld)$geneVars
        trx_per_cell <- cell.meta(object, "nUMI", cellIdx = object$dataset == d)
        nolan_constant <- mean((1 / trx_per_cell))
        alphathresh.corrected <- alpha.thresh / nrow(ld)
        geneMeanUpper <- geneMeans +
            stats::qnorm(1 - alphathresh.corrected / 2) *
            sqrt(geneMeans * nolan_constant / ncol(ld))
        basegenelower <- log10(geneMeans * nolan_constant)
        num_varGenes <- function(x, num.genes.des) {
            # This function returns the difference between the desired number of
            # genes and the number actually obtained when thresholded on x
            y <- length(which(geneVars / nolan_constant > geneMeanUpper &
                                  log10(geneVars) > basegenelower + x))
            return(abs(num.genes.des - y))
        }
        if (d %in% datasets.use) {
            ## Make the selection for shared var features ####
            if (!is.null(num.genes)) {
                # Optimize to find value of x which gives the desired number of
                # genes for this dataset if very small number of genes
                # requested, `var.thresh` may need to exceed 1
                optimized <- stats::optimize(num_varGenes, c(0, 1.5), tol = tol,
                                             num.genes.des = num.genes[i])
                var.thresh[i] <- optimized$minimum
                if (is.na(optimized$objective))
                    warning("Cannot optimize the number of selected genes for ",
                            "dataset \"", d, "\"")
                else
                    if (optimized$objective > 1)
                        warning("Returned number of genes for dataset ", d,
                                " differs from requested by ",
                                optimized$objective, ". Lower `tol` or ",
                                "`alpha.thresh` for better results.")
            }
            selected <- rownames(ld)[
                geneVars / nolan_constant > geneMeanUpper &
                    log10(geneVars) > basegenelower + var.thresh[i]
            ]
            feature.meta(ld, check = FALSE)$is_variable <-
                rownames(ld) %in% selected
            datasets(object, check = FALSE)[[d]] <- ld
            selected <- selected[selected %in% shared.features]
            if (isTRUE(verbose))
                .log("  ", length(selected),
                     " shared variable features selected")
            selected.shared[[d]] <- selected
        }
        if (d %in% unshared.datasets) {
            ## Make the selection for unshared var features ####
            if (nrow(ld) == length(shared.features)) {
                warning('No unshared feature exists for dataset "', d, '"')
            }
            selected <- rownames(ld)[
                geneVars / nolan_constant > geneMeanUpper &
                    log10(geneVars) > basegenelower + unshared.thresh[j]
            ]
            selected <- selected[!selected %in% shared.features]
            if (isTRUE(verbose))
                .log("  ", length(selected), " unshared features selected")
            selected.unshared[[d]] <- selected
        }
    }
    # Final process for shared variable features ####
    if (combine == "intersection")
        selected.shared <- Reduce(intersect, selected.shared)
    else selected.shared <- Reduce(union, selected.shared)

    if (length(selected.shared) == 0) {
        warning("No genes were selected. Lower `var.thresh` values or set ",
                '`combine = "union"`', immediate. = TRUE)
    } else {
        if (isTRUE(verbose))
            .log("Finally ", length(selected.shared),
                    " shared variable features selected.")
    }
    var.features(object, check = FALSE) <- selected.shared

    # Final process for unshared variable features ####
    for (d in unshared.datasets) {
        datasets(object, check = FALSE)[[d]]@var.unshared.features <-
            selected.unshared[[d]]
    }

    return(object)
}


#' Calculate Gene Variance for ligerDataset object
#' @param object ligerDataset object
#' @param chunkSize Integer for the maximum number of cells in each chunk.
#' Default \code{1000}.
#' @param verbose Logical. Whether to show a progress bar. Default \code{TRUE}.
#' @return The input \code{object} with calculated var updated in the H5 file.
#' @noRd
#' @useDynLib rliger, .registration = TRUE
calcGeneVars.H5 <- function(object, chunkSize = 1000, verbose = TRUE,
                            rerun = FALSE) {
    h5file <- getH5File(object)
    geneVars <- rep(0, nrow(object))
    geneMeans <- h5file[["gene_means"]][]
    safeH5Create(
        object = object,
        dataPath = "gene_vars",
        dims = nrow(object),
        dtype = "double"
    )
    geneVars <- H5Apply(
        object,
        function(chunk, sparseXIdx, cellIdx, values) {
            values + sumSquaredDeviations(chunk, geneMeans)
        },
        init = geneVars,
        useData = "norm.data",
        chunkSize = chunkSize,
        verbose = verbose
    )
    geneVars <- geneVars / (ncol(object) - 1)
    h5file[["gene_vars"]][seq_along(geneVars)] <- geneVars
    feature.meta(object, check = FALSE)$geneVars <- geneVars
    object
}

#' Plot the variance vs mean of feature expression
#' @description For each dataset where the feature variablitity is calculated,
#' a plot of log10 feature expression variance and log10 mean will be produced.
#' Features that are considered as variable would be highlighted in red.
#' @param object \linkS4class{liger} object. \code{\link{selectGenes}} need to
#' run in advance.
#' @param combinePlot Logical. If \code{TRUE}, sub-figures for all datasets will
#' be combined into one plot. if \code{FALSE}, a list of plots will be returned.
#' Default \code{TRUE}.
#' @param dotSize Controls the size of dots in the main plot. Default
#' \code{0.8}.
#' @param ... More theme setting parameters passed to
#' \code{\link{.ggplotLigerTheme}}.
#' @return \code{ggplot} object when \code{combinePlot = TRUE}, a list of
#' \code{ggplot} objects when \code{combinePlot = FALSE}
#' @export
plotVarFeatures <- function(
        object,
        combinePlot = TRUE,
        dotSize = 1,
        ...
) {
    plotList <- list()
    for (d in names(object)) {
        ld <- dataset(object, d)
        trx_per_cell <- cell.meta(object, "nUMI", cellIdx = object$dataset == d)
        nolan_constant <- mean((1 / trx_per_cell))

        data <- as.data.frame(feature.meta(ld))
        nSelect <- sum(data$is_variable)
        data$geneMeans <- log10(data$geneMeans)
        data$geneVars <- log10(data$geneVars)
        data$is_variable <- factor(data$is_variable,
                                   levels = c("TRUE", "FALSE"))
        p <- ggplot2::ggplot(data, ggplot2::aes_string(x = "geneMeans",
                                                       y = "geneVars",
                                                       color = "is_variable")) +
            ggplot2::geom_point(size = dotSize, stroke = 0) +
            ggplot2::scale_colour_discrete(name = "Variable\nfeature",
                                           type = c(`TRUE` = "red",
                                                    `FALSE` = "black")) +
            ggplot2::geom_abline(intercept = log10(nolan_constant), slope = 1,
                                 color = "purple")
        p <- .ggplotLigerTheme(p, title = d,
                               subtitle = paste0(nSelect, " variable features"),
                               xlab = "Gene Expression Mean (log10)",
                               ylab = "Gene Expression Variance (log10)",
                               ...)
        plotList[[d]] <- p
    }
    if (isTRUE(combinePlot)) {
        legend <- cowplot::get_legend(plotList[[1]])
        plotList <- lapply(plotList, function(x) {
            x + ggplot2::theme(legend.position = "none")
        })
        combined <- cowplot::plot_grid(plotlist = plotList,
                                       align = "hv",
                                       axis = "tblr")
        combined <- cowplot::plot_grid(combined, legend, rel_widths = c(5,1))
        return(combined)
    }
    else return(plotList)
}

