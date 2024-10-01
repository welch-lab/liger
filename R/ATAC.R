#' Impute the peak counts from gene expression data referring to an ATAC dataset
#' after integration
#' @description
#' This function is designed for creating peak data for a dataset with only gene
#' expression. This function uses aligned cell factor loading to find nearest
#' neighbors between cells from the queried dataset (without peak) and cells
#' from reference dataset (with peak). And then impute the peak for the former
#' basing on the weight. Therefore, the reference dataset selected must be of
#' "atac" modality setting.
#' @param object \linkS4class{liger} object with aligned factor loading computed
#' in advance.
#' @param nNeighbors The maximum number of nearest neighbors to search. Default
#' \code{20}.
#' @param reference Name of a dataset containing peak data to impute into query
#' dataset(s).
#' @param queries Names of datasets to be augmented by imputation. Should not
#' include \code{reference}. Default \code{NULL} uses all datasets except the
#' reference.
#' @param weight Logical. Whether to use KNN distances as weight matrix. Default
#' \code{FALSE}.
#' @param norm Logical. Whether to normalize the imputed data. Default
#' \code{TRUE}.
#' @param scale Logical. Whether to scale but not center the imputed data.
#' Default \code{TRUE}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param ... Optional arguments to be passed to \code{\link{normalize}} when
#' \code{norm = TRUE}.
#' @param knn_k \bold{Deprecated}. See Usage section for replacement.
#' @return The input \code{object} where queried \linkS4class{ligerDataset}
#' objects in \code{datasets} slot are replaced. These datasets will all be
#' converted to \linkS4class{ligerATACDataset} class with an additional slot
#' \code{rawPeak} to store the imputed peak counts, and \code{normPeak} for
#' normalized imputed peak counts if \code{norm = TRUE}.
#' @export
#' @examples
#' bmmc <- normalize(bmmc)
#' bmmc <- selectGenes(bmmc, datasets.use = "rna")
#' bmmc <- scaleNotCenter(bmmc)
#' if (requireNamespace("RcppPlanc", quietly = TRUE)) {
#'     bmmc <- runINMF(bmmc, k = 20)
#'     bmmc <- alignFactors(bmmc)
#'     bmmc <- normalizePeak(bmmc)
#'     bmmc <- imputeKNN(bmmc, reference = "atac", queries = "rna")
#' }
imputeKNN <- function(
        object,
        reference,
        queries = NULL,
        nNeighbors = 20,
        weight = TRUE,
        norm = TRUE,
        scale = FALSE,
        verbose = getOption("ligerVerbose", TRUE),
        ...,
        # Deprecated coding style,
        knn_k = nNeighbors
) {
    .deprecateArgs(list(knn_k = "nNeighbors"), defunct = "scale")
    if (is.null(getMatrix(object, "H.norm")))
        cli::cli_abort(
            "Aligned factor loading has to be available for imputation.
            Please run {.fn alignFactors} in advance.")
    reference <- .checkArgLen(reference, n = 1)
    reference <- .checkUseDatasets(object, reference)#, modal = "atac")
    queries <- .checkUseDatasets(object, queries)
    if (any(queries %in% reference)) {
        cli::cli_alert_warning(
            "Reference dataset cannot be inclued in the query datasets."
        )
        cli::cli_alert_warning(
            "Removed from query list: {.val {queries[queries %in% reference]}}"
        )
        queries <- queries[!queries %in% reference]
    }
    object <- recordCommand(object, ..., dependencies = c("RANN", "Matrix"))
    if (isTRUE(verbose)) {
        cli::cli_alert_info(
            "Imputing {length(queries)} query dataset{?s}: {.val {queries}}"
        )
        cli::cli_alert_info("from reference dataset: {.val {reference}}")
    }

    referenceCells <- colnames(dataset(object, reference))
    for (i in seq_along(queries)) {
        query <- queries[i]
        queryLD <- dataset(object, query)
        queryCells <- colnames(queryLD)

        # creating a (reference cell numbers X query cell numbers) weights
        # matrix for knn weights and unit weights
        # knn <- FNN::get.knnx(object@H.norm[referenceCells, ],
        #                       object@H.norm[queryCells, ],
        #                       k = nNeighbors,
        #                       algorithm = "CR")
        knn <- RANN::nn2(object@H.norm[referenceCells, ],
                         object@H.norm[queryCells, ],
                         k = nNeighbors)

        weights <- Matrix::Matrix(0, nrow = length(referenceCells),
                                  ncol = nrow(knn$nn.idx), sparse = TRUE)
        if (isTRUE(weight)) {
            # for weighted situation
            # find nearest neighbors for query cell in normed ref datasets
            for (n in seq_len(nrow(knn$nn.idx))) {
                # record ref-query cell-cell distances
                weights[knn$nn.idx[n, ], n] <-
                    exp(-knn$nn.dists[n, ]) / sum(exp(-knn$nn.dists[n, ]))
            }
        } else{
            # for unweighted situation
            for (n in seq_len(nrow(knn$nn.idx))) {
                # simply count the mean
                weights[knn$nn.idx[n, ], n] <- 1 / nNeighbors
            }
        }

        # (genes by ref cell num) multiply by the weight matrix
        # (ref cell num by query cell num)
        referenceRawData <- rawPeak(object, reference)
        imputed_vals <- referenceRawData %*% weights
        # assigning dimnames
        colnames(imputed_vals) <- queryCells
        rownames(imputed_vals) <- rownames(referenceRawData)

        if (!inherits(imputed_vals, "dgCMatrix"))
            imputed_vals <- methods::as(imputed_vals, "dgCMatrix")
        newQueryLD <- as.ligerDataset(queryLD, modal = "atac")
        rawPeak(newQueryLD) <- imputed_vals
        datasets(object, check = FALSE)[[query]] <- newQueryLD
    }

    if (isTRUE(norm)) {
        object <- normalizePeak(object, useDatasets = queries,
                                verbose = verbose, ...)
    }

    return(object)
}

#' Linking genes to putative regulatory elements
#'
#' @description Evaluate the relationships between pairs of genes and peaks
#' based on specified distance metric. Usually used for inferring the
#' correlation between gene expression and imputed peak counts for datasets
#' without the modality originally (i.e. applied to \code{\link{imputeKNN}}
#' result).
#' @param object A \linkS4class{liger} object, with datasets that is of
#' \linkS4class{ligerATACDataset} class in the \code{datasets} slot.
#' @param pathToCoords Path tothe gene coordinates file, usually a BED file.
#' @param useDataset Name of one dataset, with both normalized gene expression
#' and normalized peak counts available.
#' @param useGenes Character vector of gene names to be tested. Default
#' \code{NULL} uses all genes available in \code{useDataset}.
#' @param method Choose the type of correlation to calculate, from
#' \code{"spearman"}, \code{"pearson"} and \code{"kendall"}. Default
#' \code{"spearman"}
#' @param alpha Numeric, significance threshold for correlation p-value.
#' Peak-gene correlations with p-values below this threshold are considered
#' significant. Default \code{0.05}.
#' @param verbose Logical. Whether to show information of the progress. Default
#' \code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.
#' @param path_to_coords,genes.list,dist \bold{Deprecated}. See Usage section
#' for replacement.
#' @return A sparse matrix with peak names as rows and gene names as columns,
#' with each element indicating the correlation between peak i and gene j, 0 if
#' the gene and peak are not significantly linked.
#' @seealso \code{\link{imputeKNN}}
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("RcppPlanc", quietly = TRUE) &&
#'     requireNamespace("GenomicRanges", quietly = TRUE) &&
#'     requireNamespace("IRanges", quietly = TRUE) &&
#'     requireNamespace("psych", quietly = TRUE)) {
#'     bmmc <- normalize(bmmc)
#'     bmmc <- selectGenes(bmmc)
#'     bmmc <- scaleNotCenter(bmmc)
#'     bmmc <- runINMF(bmmc, miniBatchSize = 100)
#'     bmmc <- alignFactors(bmmc)
#'     bmmc <- normalizePeak(bmmc)
#'     bmmc <- imputeKNN(bmmc, reference = "atac", queries = "rna")
#'     corr <- linkGenesAndPeaks(
#'         bmmc, useDataset = "rna",
#'         pathToCoords = system.file("extdata/hg19_genes.bed", package = "rliger")
#'     )
#' }
#' }
linkGenesAndPeaks <- function(
        object,
        useDataset,
        pathToCoords,
        useGenes = NULL,
        method = c("spearman", "pearson", "kendall"),
        alpha = 0.05,
        verbose = getOption("ligerVerbose", TRUE),
        # Deprecated coding style
        path_to_coords = pathToCoords,
        genes.list = useGenes,
        dist = method
) {
    ## check dependency
    if (!requireNamespace("GenomicRanges", quietly = TRUE))
        cli::cli_abort(
            "Package {.pkg GenomicRanges} is needed for this function to work.
            Please install it by command:
            {.code BiocManager::install('GenomicRanges')}")

    if (!requireNamespace("IRanges", quietly = TRUE))
        cli::cli_abort(
            "Package {.pkg IRanges} is needed for this function to work.
            Please install it by command:
            {.code BiocManager::install('IRanges')}")
    if (!requireNamespace("psych", quietly = TRUE))
        cli::cli_abort(
            "Package {.pkg psych} is needed for this function to work.
            Please install it by command:
            {.code BiocManager::install('psych')}")
    .deprecateArgs(list(path_to_coords = "pathToCoords",
                        genes.list = "useGenes", dist = "method"))
    method <- match.arg(method)
    useDataset <- .checkUseDatasets(object, useDataset)
    if (length(useDataset) != 1)
        cli::cli_abort("Please select only one dataset")
    lad <- dataset(object, useDataset)
    if (!inherits(lad, "ligerATACDataset"))
        cli::cli_abort(
            "Specified dataset is not of `ligerATACDataset` class.
            Please try {.fn imputeKNN} with `query = '{useDataset}'`.")
    if (is.null(normData(lad)))
        cli::cli_abort("Normalized gene expression not found in specified dataset.")
    if (is.null(normPeak(lad)))
        cli::cli_abort("Normalized peak counts not found in specified dataset.")

    ### make GRanges object for peaks
    peakCounts <- normPeak(lad)
    geneCounts <- normData(lad)

    regions <- .splitPeakRegion(rownames(peakCounts))
    peaksPos <- GenomicRanges::GRanges(
        seqnames = regions$chrs,
        ranges = IRanges::IRanges(regions$chrsStart, end = regions$chrsEnd)
    )

    ### make GRanges object for genes
    geneNames <- utils::read.csv2(pathToCoords, sep = "\t",
                                  header = FALSE, stringsAsFactors = FALSE)

    geneNames <- geneNames[stats::complete.cases(geneNames), ]
    genesCoords <- GenomicRanges::GRanges(
        seqnames = geneNames$V1,
        ranges = IRanges::IRanges(as.numeric(geneNames$V2),
                                  end = as.numeric(geneNames$V3))
    )
    names(genesCoords) <- geneNames$V4

    ### Data clean-up
    # cell x genes, because psych::corr.test requires input matrix
    # with obs as rows
    geneCounts <- t(geneCounts)
    # cell x peaks
    peakCounts <- t(peakCounts)

    # find overlap peaks for each gene
    if (is.null(useGenes)) useGenes <- colnames(geneCounts)
    missingGenes <- !useGenes %in% names(genesCoords)
    if (sum(missingGenes) != 0 && isTRUE(verbose))
        cli::cli_alert_warning(
            "Ignoring {sum(missingGenes)} genes not found in given gene coordinates"
        )
    useGenes <- useGenes[!missingGenes]
    if (length(useGenes) == 0) {
        cli::cli_abort(
            "Number of genes to be tested equals 0. Please check input
            {.code useGenes} or the coordinate file."
        )
    } else {
        cli::cli_alert_info(
            "{length(useGenes)} genes to be tested against {ncol(peakCounts)} peaks"
        )
    }
    genesCoords <- genesCoords[useGenes]

    ### construct regnet
    if (isTRUE(verbose)) {
        cli::cli_alert_info("Calculating correlation for gene-peak pairs...")
        cli::cli_progress_bar("", total = length(useGenes), type = "iter")
    }

    # Result would be a sparse matrix, initialize the `i`, `p`, `x` vectors.
    ind <- numeric()
    indp <- numeric()
    values <- numeric()
    eachLen <- 0
    for (pos in seq_along(useGenes)) {
        geneUse <- useGenes[pos]
        # re-scale the window for each gene
        gene.loci <- GenomicRanges::trim(suppressWarnings(
            GenomicRanges::promoters(
                GenomicRanges::resize(genesCoords[geneUse],
                                      width = 1, fix = "start"),
                upstream = 500000,
                downstream = 500000
            )
        ))
        peaks.use <- S4Vectors::queryHits(
            GenomicRanges::findOverlaps(peaksPos, gene.loci)
        )
        if (length(peaks.use) == 0) {
            # if no peaks in window, skip this iteration
            indp <- c(indp, as.numeric(eachLen))
            next
        }

        ### compute correlation and p-adj for genes and peaks ###
        res <- suppressWarnings(psych::corr.test(
            x = geneCounts[, geneUse],
            y = as.matrix(peakCounts[, peaks.use]),
            method = method,
            adjust = "holm",
            ci = FALSE,
            use = "complete"
        ))
        # filter by p-value
        pick <- res$p < alpha
        pick[is.na(pick)] <- FALSE
        res.corr <- as.numeric(res$r[pick])
        peaks.use <- peaks.use[pick]

        eachLen <- eachLen + length(peaks.use)
        ind <- c(ind, as.numeric(peaks.use))
        indp <- c(indp, as.numeric(eachLen))
        values <- c(values, res.corr)
        if (isTRUE(verbose)) {
            cli::cli_progress_update(set = pos)
        }
    }
    # make final sparse matrix
    regnet <-  Matrix::sparseMatrix(
        i = ind, p = c(0, indp), x = values,
        dims = c(ncol(peakCounts), length(useGenes)),
        dimnames = list(colnames(peakCounts), useGenes)
    )

    return(regnet)
}

#' Export predicted gene-pair interaction
#' @description Export the predicted gene-pair interactions calculated by
#' upstream function \code{\link{linkGenesAndPeaks}} into an Interact Track file
#' which is compatible with \href{https://genome.ucsc.edu/cgi-bin/hgCustom}{UCSC
#' Genome Browser}.
#' @param corrMat A sparse matrix of correlation with peak names as rows and
#' gene names as columns.
#' @param pathToCoords Path to the gene coordinates file.
#' @param useGenes Character vector of gene names to be exported. Default
#' \code{NULL} uses all genes available in \code{corrMat}.
#' @param outputPath Path of filename where the output file will be stored. If
#' a folder, a file named \code{"Interact_Track.bed"} will be created. Default
#' current working directory.
#' @return No return value. A file located at \code{outputPath} will be created.
#' @export
#' @examples
#' \donttest{
#' bmmc <- normalize(bmmc)
#' bmmc <- selectGenes(bmmc)
#' bmmc <- scaleNotCenter(bmmc)
#' if (requireNamespace("RcppPlanc", quietly = TRUE) &&
#'     requireNamespace("GenomicRanges", quietly = TRUE) &&
#'     requireNamespace("IRanges", quietly = TRUE) &&
#'     requireNamespace("psych", quietly = TRUE)) {
#'     bmmc <- runINMF(bmmc)
#'     bmmc <- alignFactors(bmmc)
#'     bmmc <- normalizePeak(bmmc)
#'     bmmc <- imputeKNN(bmmc, reference = "atac", queries = "rna")
#'     corr <- linkGenesAndPeaks(
#'         bmmc, useDataset = "rna",
#'         pathToCoords = system.file("extdata/hg19_genes.bed", package = "rliger")
#'     )
#'     resultPath <- tempfile()
#'     exportInteractTrack(
#'         corrMat = corr,
#'         pathToCoords = system.file("extdata/hg19_genes.bed", package = "rliger"),
#'         outputPath = resultPath
#'     )
#'     head(read.table(resultPath, skip = 1))
#' }
#' }
exportInteractTrack <- function(
        corrMat,
        pathToCoords,
        useGenes = NULL,
        outputPath = getwd()
) {
    # check useGenes
    if (is.null(useGenes)) {
        useGenes <- colnames(corrMat)
    } else if (any(!useGenes %in% colnames(corrMat))) {
        cli::cli_alert_warning(
            "Removed {sum(!useGenes %in% colnames(corrMat))} genes not found in {.code corrMat}"
        )
        useGenes <- useGenes[useGenes %in% colnames(corrMat)]
    }
    # Filter useGenes by significance
    geneSel <- Matrix::colSums(corrMat[, useGenes, drop = FALSE] != 0) > 0
    if (length(useGenes) - sum(geneSel) > 0)
        cli::cli_alert_warning(
            "Totally {length(useGenes) - sum(geneSel)} selected genes do not have significant correlated peaks, out of {length(useGenes)} selected genes",
            wrap = TRUE
        )
    useGenes <- useGenes[geneSel]
    if (length(useGenes) == 0) {
        cli::cli_abort("No gene requested is either available or having
                       significant correlated peaks. ")
    }

    ### make GRanges object for genes
    genesCoords <- utils::read.csv2(
        pathToCoords, sep = "\t", header = FALSE,
        colClasses = c("character", "integer", "integer",
                       "character", "NULL", "NULL")
    )
    genesCoords <- genesCoords[stats::complete.cases(genesCoords$V4), ]
    rownames(genesCoords) <- genesCoords[, 4]
    # split peak names into chrom and coordinates
    regions <- .splitPeakRegion(rownames(corrMat))

    # check output_path
    if (dir.exists(outputPath)) {
        # If it happens to be a directory
        outputPath <- file.path(outputPath, "Interact_Track.bed")
    }
    if (!file.exists(outputPath)) file.create(outputPath)
    outputPath <- normalizePath(outputPath)

    # Start writing BED file
    trackDoc <- paste0('track type=interact name="Interaction Track" ',
                       'description="Gene-Peaks Links" ',
                       'interactDirectional=true maxHeightPixels=200:100:50 ',
                       'visibility=full')
    write(trackDoc, file = outputPath)
    for (gene in useGenes) {
        peaksSel <- corrMat[, gene] != 0
        track <- data.frame(
            chrom = regions$chrs[peaksSel],
            chromStart = regions$chrsStart[peaksSel],
            chromEnd = regions$chrsEnd[peaksSel],
            name = paste0(gene, "/", rownames(corrMat)[peaksSel]),
            score = 0,
            value = as.numeric(corrMat[peaksSel, gene]),
            exp = ".",
            color = 5,
            sourceChrom = regions$chrs[peaksSel],
            sourceStart = regions$chrsStart[peaksSel],
            sourceEnd = regions$chrsStart[peaksSel] + 1,
            sourceName = ".",
            sourceStrand = ".",
            targetChrom = genesCoords[gene, 1],
            targetStart = genesCoords[gene, 2],
            targetEnd = genesCoords[gene, 2] + 1,
            targetName = gene,
            targetStrand = "."
        )
        utils::write.table(
            track,
            file = outputPath,
            append = TRUE,
            quote = FALSE,
            sep = "\t",
            eol = "\n",
            na = "NA",
            dec = ".",
            row.names = FALSE,
            col.names = FALSE,
            qmethod = c("escape", "double"),
            fileEncoding = ""
        )
    }
    cli::cli_alert_success("Result written at: {.file {outputPath}}")
    invisible(NULL)
}

#' `r lifecycle::badge("deprecated")` Export predicted gene-pair interaction
#' @description Export the predicted gene-pair interactions calculated by
#' upstream function \code{\link{linkGenesAndPeaks}} into an Interact Track file
#' which is compatible with \href{https://genome.ucsc.edu/cgi-bin/hgCustom}{UCSC
#' Genome Browser}.
#' @param corr.mat A sparse matrix of correlation with peak names as rows and
#' gene names as columns.
#' @param path_to_coords Path to the gene coordinates file.
#' @param genes.list Character vector of gene names to be exported. Default
#' \code{NULL} uses all genes available in \code{corrMat}.
#' @param output_path Path of filename where the output file will be stored. If
#' a folder, a file named \code{"Interact_Track.bed"} will be created. Default
#' current working directory.
#' @return No return value. A file located at \code{outputPath} will be created.
#' @name makeInteractTrack-deprecated
#' @seealso \code{\link{rliger-deprecated}}, \code{\link{exportInteractTrack}}
NULL

#' @rdname rliger-deprecated
#' @section \code{makeInteractTrack}:
#' For \code{makeInteractTrack}, use \code{\link{exportInteractTrack}}.
#' @export
makeInteractTrack <- function(
        corr.mat,
        path_to_coords,
        genes.list = NULL,
        output_path = getwd()
) {
    lifecycle::deprecate_warn("1.99.0", "makeInteractTrack()",
                              "exportInteractTrack()")
    exportInteractTrack(corrMat = corr.mat, pathToCoords = path_to_coords,
                        useGenes = genes.list, outputPath = output_path)
}

.splitPeakRegion <- function(peakNames) {
    peakNames <- strsplit(peakNames, "[:-]")
    chrs <- Reduce(append, lapply(peakNames, function(peak) peak[1]))
    chrsStart <- Reduce(append, lapply(peakNames, function(peak) peak[2]))
    chrsEnd <- Reduce(append, lapply(peakNames, function(peak) peak[3]))
    list(chrs = chrs,
         chrsStart = as.numeric(chrsStart),
         chrsEnd = as.numeric(chrsEnd))
}
