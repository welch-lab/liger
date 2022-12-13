#' @importFrom Matrix colSums rowSums t
NULL

#' Normalize raw datasets to column sums
#'
#' This function normalizes data to account for total gene expression across a
#' cell.
#'
#' @param object \code{liger} object.
#' @param chunk size of chunks in hdf5 file. (default 1000)
#' @param format.type string of HDF5 format (10X CellRanger by default).
#' @param remove.missing Whether to remove cells not expressing any measured
#' genes, and genes not expressed in any cells (if take.gene.union = TRUE,
#' removes only genes not expressed in any dataset) (default TRUE).
#' @param verbose Print progress bar/messages (TRUE by default)
#' @return \code{liger} object with norm.data slot set.
#' @export
#' @examples
#' # Demonstration using matrices with randomly generated numbers
#' Y <- matrix(runif(5000,0,2), 10,500)
#' Z <- matrix(runif(5000,0,2), 10,500)
#' ligerex <- createLiger(list(y_set = Y, z_set = Z))
#' ligerex <- normalize(ligerex)
normalizeOld <- function(object,
                      chunk = 1000,
                      format.type = "10X",
                      remove.missing = TRUE,
                      verbose = TRUE) {
    if (class(object@raw.data[[1]])[1] == "H5File") {
        hdf5_files = names(object@raw.data)
        nUMI = c()
        nGene = c()
        for (i in 1:length(hdf5_files))
        {
            if (verbose) {
                message(hdf5_files[i])
            }
            chunk_size = chunk
            #fname = hdf5_files[[i]]
            num_entries = object@h5file.info[[i]][["data"]]$dims
            num_cells = object@h5file.info[[i]][["barcodes"]]$dims
            num_genes = object@h5file.info[[i]][["genes"]]$dims


            prev_end_col = 1
            prev_end_data = 1
            prev_end_ind = 0
            gene_sum_sq = rep(0, num_genes)
            gene_means = rep(0, num_genes)
            #file.h5$close_all()

            safe_h5_create(
                object = object,
                idx = i,
                dataset_name = "/norm.data",
                dims = num_entries,
                mode = hdf5r::h5types$double,
                chunk_size = chunk_size
            )
            safe_h5_create(
                object = object,
                idx = i,
                dataset_name = "/cell_sums",
                dims = num_cells,
                mode = hdf5r::h5types$int,
                chunk_size = chunk_size
            )

            #file.h5 = H5File$new(fname, mode="r+")
            num_chunks = ceiling(num_cells / chunk_size)
            if (verbose) {
                pb = txtProgressBar(0, num_chunks, style = 3)
            }
            ind = 0
            while (prev_end_col < num_cells)
            {
                ind = ind + 1
                if (num_cells - prev_end_col < chunk_size)
                {
                    chunk_size = num_cells - prev_end_col + 1
                }
                start_inds = object@h5file.info[[i]][["indptr"]][prev_end_col:(prev_end_col +
                                                                                   chunk_size)]
                row_inds = object@h5file.info[[i]][["indices"]][(prev_end_ind +
                                                                     1):(tail(start_inds, 1))]
                counts = object@h5file.info[[i]][["data"]][(prev_end_ind +
                                                                1):(tail(start_inds, 1))]
                raw.data = sparseMatrix(
                    i = row_inds[1:length(counts)] + 1,
                    p = start_inds[1:(chunk_size + 1)] - prev_end_ind,
                    x = counts,
                    dims = c(num_genes, chunk_size)
                )
                nUMI = c(nUMI, colSums(raw.data))
                nGene = c(nGene, colSums(raw.data > 0))
                norm.data = Matrix.column_norm(raw.data)
                object@raw.data[[i]][["norm.data"]][(prev_end_ind + 1):(tail(start_inds, 1))] = norm.data@x
                object@raw.data[[i]][["cell_sums"]][prev_end_col:(prev_end_col +
                                                                      chunk_size - 1)] = Matrix::colSums(raw.data)
                #h5write(norm.data,file=fname,name="/norm.data",index=list(prev_end_ind:tail(start_inds, 1)))
                #h5write(colSums(raw.data),file=fname,name="/cell_sums",index=list(prev_end_col:(prev_end_col+chunk_size)))
                prev_end_col = prev_end_col + chunk_size
                prev_end_data = prev_end_data + length(norm.data@x)
                prev_end_ind = tail(start_inds, 1)

                # calculate row sum and sum of squares using normalized data
                row_sums = Matrix::rowSums(norm.data)
                gene_sum_sq = gene_sum_sq + rowSums(norm.data * norm.data)
                gene_means = gene_means + row_sums
                if (verbose) {
                    setTxtProgressBar(pb, ind)
                }
            }
            if (verbose) {
                setTxtProgressBar(pb, num_chunks)
                cat("\n")
            }
            gene_means = gene_means / num_cells
            safe_h5_create(
                object = object,
                idx = i,
                dataset_name = "gene_means",
                dims = num_genes,
                mode = hdf5r::h5types$double
            )
            safe_h5_create(
                object = object,
                idx = i,
                dataset_name = "gene_sum_sq",
                dims = num_genes,
                mode = hdf5r::h5types$double
            )
            object@raw.data[[i]][["gene_means"]][1:length(gene_means)] = gene_means
            object@raw.data[[i]][["gene_sum_sq"]][1:length(gene_sum_sq)] = gene_sum_sq
            object@norm.data[[i]] = object@raw.data[[i]][["norm.data"]]
            rm(row_sums)
            rm(raw.data)
        }
        object@cell.data$nUMI = nUMI
        object@cell.data$nGene = nGene

        for (i in 1:length(object@raw.data)) {
            if (!object@raw.data[[i]]$exists("cell.data")) {
                cell.data.i = object@cell.data[object@cell.data$dataset == names(object@raw.data)[i],]
                cell.data.i$barcode = rownames(cell.data.i)
                object@raw.data[[i]][["cell.data"]] = cell.data.i
            }
        }

        names(object@norm.data) = names(object@raw.data)
    } else {
        if (remove.missing) {
            object <-
                removeMissingObs(object, slot.use = "raw.data", use.cols = TRUE)
        }
        if (class(object@raw.data[[1]])[1] == "dgTMatrix" |
            class(object@raw.data[[1]])[1] == "dgCMatrix") {
            object@norm.data <- lapply(object@raw.data, Matrix.column_norm)
        } else {
            object@norm.data <- lapply(object@raw.data, function(x) {
                sweep(x, 2, colSums(x), "/")
            })
        }
    }
    return(object)
}

#' Calculate variance of gene expression across cells in an online fashion
#'
#' This function calculates the variance of gene expression values across cells
#' for hdf5 files.
#'
#' @param object \code{liger} object. The input raw.data should be a list of
#' hdf5 files. Should call normalize and selectGenes before calling.
#' @param chunk size of chunks in hdf5 file. (default 1000)
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return \code{liger} object with scale.data slot set.
calcGeneVars = function(object,
                        chunk = 1000,
                        verbose = TRUE)
{
    hdf5_files = names(object@raw.data)
    for (i in 1:length(hdf5_files)) {
        if (verbose) {
            message(hdf5_files[i])
        }
        chunk_size = chunk
        num_cells = object@h5file.info[[i]][["barcodes"]]$dims
        num_genes = object@h5file.info[[i]][["genes"]]$dims
        num_entries = object@h5file.info[[i]][["data"]]$dims

        prev_end_col = 1
        prev_end_data = 1
        prev_end_ind = 0
        gene_vars = rep(0, num_genes)
        gene_means = object@raw.data[[i]][["gene_means"]][]
        gene_num_pos = rep(0, num_genes)

        num_chunks = ceiling(num_cells / chunk_size)
        if (verbose) {
            pb = txtProgressBar(0, num_chunks, style = 3)
        }
        ind = 0
        while (prev_end_col < num_cells) {
            ind = ind + 1
            if (num_cells - prev_end_col < chunk_size) {
                chunk_size = num_cells - prev_end_col + 1
            }
            start_inds = object@h5file.info[[i]][["indptr"]][prev_end_col:(prev_end_col +
                                                                               chunk_size)]
            row_inds = object@h5file.info[[i]][["indices"]][(prev_end_ind +
                                                                 1):(tail(start_inds, 1))]
            counts = object@norm.data[[i]][(prev_end_ind + 1):(tail(start_inds, 1))]
            norm.data = sparseMatrix(
                i = row_inds[1:length(counts)] + 1,
                p = start_inds[1:(chunk_size + 1)] - prev_end_ind,
                x = counts,
                dims = c(num_genes, chunk_size)
            )

            num_read = length(counts)
            prev_end_col = prev_end_col + chunk_size
            prev_end_data = prev_end_data + num_read
            prev_end_ind = tail(start_inds, 1)
            gene_vars = gene_vars + sumSquaredDeviations(norm.data, gene_means)
            if (verbose) {
                setTxtProgressBar(pb, ind)
            }
        }
        if (verbose) {
            setTxtProgressBar(pb, num_chunks)
            cat("\n")
        }
        gene_vars = gene_vars / (num_cells - 1)
        safe_h5_create(
            object = object,
            idx = i,
            dataset_name = "/gene_vars",
            dims = num_genes,
            mode = hdf5r::h5types$double
        )
        object@raw.data[[i]][["gene_vars"]][1:num_genes] = gene_vars
    }
    return(object)
}

#' Select a subset of informative genes
#'
#' This function identifies highly variable genes from each dataset and combines
#' these gene sets (either by union or intersection) for use in downstream
#' analysis. Assuming that gene expression approximately follows a Poisson
#' distribution, this function identifies genes with gene expression variance
#' above a given variance threshold (relative to mean gene expression). It also
#' provides a log plot of gene variance vs gene expression (with a line
#' indicating expected expression across genes and cells). Selected genes are
#' plotted in green.
#'
#' @param object \code{liger} object. Should have already called normalize.
#' @param var.thresh Variance threshold. Main threshold used to identify
#' variable genes. Genes with expression variance greater than threshold
#' (relative to mean) are selected. (higher threshold -> fewer selected genes).
#' Accepts single value or vector with separate var.thresh for each dataset.
#' (default 0.1)
#' @param alpha.thresh Alpha threshold. Controls upper bound for expected mean
#' gene expression (lower threshold -> higher upper bound). (default 0.99)
#' @param num.genes Number of genes to find for each dataset. Optimises the
#' value of var.thresh for each dataset to get this number of genes. Accepts
#' single value or vector with same length as number of datasets (optional,
#' default=NULL).
#' @param tol Tolerance to use for optimization if num.genes values passed in
#' (default 0.0001).
#' @param datasets.use List of datasets to include for discovery of highly
#' variable genes. (default 1:length(object@raw.data))
#' @param combine How to combine variable genes across experiments. Either
#' "union" or "intersection". (default "union")
#' @param capitalize Capitalize gene names to match homologous genes (ie. across
#' species)
#'   (default FALSE)
#' @param do.plot Display log plot of gene variance vs. gene expression for each
#' dataset. Selected genes are plotted in green. (default FALSE)
#' @param cex.use Point size for plot.
#' @param chunk size of chunks in hdf5 file. (default 1000)
#' @param unshared.features Whether to consider unshared features
#' @param unshared.datasets A list of the datasets to consider unshared features
#' for, i.e. list(2), to use the second dataset
#' @param unshared.thresh A list of threshold values to apply to each unshared
#' dataset. If only one value is provided, it will apply to all unshared
#' datasets. If a list is provided, it must match the length of the unshared
#' datasets submitted.
#' @return \code{liger} object with var.genes slot set.
#'
#' @importFrom stats optimize
#' @importFrom graphics abline plot points title
#' @importFrom stats qnorm
#'
#' @export
#' @examples
#' \dontrun{
#' # Given datasets Y and Z
#' ligerex <- createLiger(list(y_set = Y, z_set = Z))
#' ligerex <- normalize(ligerex)
#' # use default selectGenes settings (var.thresh = 0.1)
#' ligerex <- selectGenes(ligerex)
#' # select a smaller subset of genes
#' ligerex <- selectGenes(ligerex, var.thresh = 0.3)
#' }

selectGenes <-
    function(object,
             var.thresh = 0.1,
             alpha.thresh = 0.99,
             num.genes = NULL,
             tol = 0.0001,
             datasets.use = 1:length(object@raw.data),
             combine = "union",
             capitalize = FALSE,
             do.plot = FALSE,
             cex.use = 0.3,
             chunk = 1000,
             unshared = F,
             unshared.datasets = NULL,
             unshared.thresh = NULL)
    {
        if (class(object@raw.data[[1]])[1] == "H5File") {
            if (!object@raw.data[[1]]$exists("gene_vars")) {
                object = calcGeneVars(object, chunk)
            }
            hdf5_files = names(object@raw.data)
            if (length(var.thresh) == 1) {
                var.thresh <- rep(var.thresh, length(hdf5_files))
            }
            genes.use <- c()
            for (i in datasets.use) {
                if (object@h5file.info[[i]][["format.type"]] == "AnnData") {
                    genes = object@h5file.info[[i]][["genes"]][]$index
                } else {
                    genes = object@h5file.info[[i]][["genes"]][]
                }

                if (capitalize) {
                    genes = toupper(genes)
                }
                trx_per_cell = object@raw.data[[i]][["cell_sums"]][]
                gene_expr_mean = object@raw.data[[i]][["gene_means"]][]
                gene_expr_var = object@raw.data[[i]][["gene_vars"]][]

                names(gene_expr_mean) <-
                    names(gene_expr_var) <- genes # assign gene names
                nolan_constant <- mean((1 / trx_per_cell))
                alphathresh.corrected <- alpha.thresh / length(genes)
                genemeanupper <-
                    gene_expr_mean + qnorm(1 - alphathresh.corrected / 2) *
                    sqrt(gene_expr_mean * nolan_constant / length(trx_per_cell))
                genes.new <-
                    names(gene_expr_var)[which(
                        gene_expr_var / nolan_constant >
                            genemeanupper &
                            log10(gene_expr_var) > log10(gene_expr_mean) +
                            (log10(nolan_constant) + var.thresh[i])
                    )]
                if (do.plot) {
                    plot(
                        log10(gene_expr_mean),
                        log10(gene_expr_var),
                        cex = cex.use,
                        xlab = "Gene Expression Mean (log10)",
                        ylab = "Gene Expression Variance (log10)"
                    )
                    points(
                        log10(gene_expr_mean[genes.new]),
                        log10(gene_expr_var[genes.new]),
                        cex = cex.use,
                        col = "green"
                    )
                    abline(log10(nolan_constant), 1, col = "purple")
                    legend(
                        "bottomright",
                        paste0("Selected genes: ",
                               length(genes.new)),
                        pch = 20,
                        col = "green"
                    )
                    title(main = hdf5_files[i])
                }
                if (combine == "union") {
                    genes.use <- union(genes.use, genes.new)
                }
                if (combine == "intersection") {
                    if (length(genes.use) == 0) {
                        genes.use <- genes.new
                    }
                    genes.use <- intersect(genes.use, genes.new)
                }
            }

            for (i in 1:length(hdf5_files)) {
                if (object@h5file.info[[i]][["format.type"]] == "AnnData") {
                    genes = object@h5file.info[[i]][["genes"]][]$index
                } else {
                    genes = object@h5file.info[[i]][["genes"]][]
                }
                genes.use <- genes.use[genes.use %in% genes]
            }

            if (length(genes.use) == 0) {
                warning(
                    "No genes were selected; lower var.thresh values or choose 'union' for combine parameter",
                    immediate. = TRUE
                )
            }
            object@var.genes = genes.use
        } else {
            # Expand if only single var.thresh passed
            if (length(var.thresh) == 1) {
                var.thresh <- rep(var.thresh, length(object@raw.data))
            }
            if (length(num.genes) == 1) {
                num.genes <- rep(num.genes, length(object@raw.data))
            }
            if (!identical(intersect(datasets.use, 1:length(object@raw.data)), datasets.use)) {
                datasets.use = intersect(datasets.use, 1:length(object@raw.data))
            }
            genes.use <- c()
            for (i in datasets.use) {
                if (capitalize) {
                    rownames(object@raw.data[[i]]) <-
                        toupper(rownames(object@raw.data[[i]]))
                    rownames(object@norm.data[[i]]) <-
                        toupper(rownames(object@norm.data[[i]]))
                }
                trx_per_cell <- colSums(object@raw.data[[i]])
                # Each gene's mean expression level (across all cells)
                gene_expr_mean <- rowMeansFast(object@norm.data[[i]])
                # Each gene's expression variance (across all cells)
                gene_expr_var <-
                    rowVarsFast(object@norm.data[[i]], gene_expr_mean)
                names(gene_expr_mean) <-
                    names(gene_expr_var) <- rownames(object@norm.data[[i]])
                nolan_constant <- mean((1 / trx_per_cell))
                alphathresh.corrected <-
                    alpha.thresh / nrow(object@raw.data[[i]])
                genemeanupper <-
                    gene_expr_mean + qnorm(1 - alphathresh.corrected / 2) *
                    sqrt(gene_expr_mean * nolan_constant / ncol(object@raw.data[[i]]))
                basegenelower <- log10(gene_expr_mean * nolan_constant)

                num_varGenes <- function(x, num.genes.des) {
                    # This function returns the difference between the desired number of genes and
                    # the number actually obtained when thresholded on x
                    y <-
                        length(
                            which(
                                gene_expr_var / nolan_constant > genemeanupper &
                                    log10(gene_expr_var) > basegenelower + x
                            )
                        )
                    return(abs(num.genes.des - y))
                }

                if (!is.null(num.genes)) {
                    # Optimize to find value of x which gives the desired number of genes for this dataset
                    # if very small number of genes requested, var.thresh may need to exceed 1
                    optimized <-
                        optimize(num_varGenes,
                                 c(0, 1.5),
                                 tol = tol,
                                 num.genes.des = num.genes[i])
                    var.thresh[i] <- optimized$minimum
                    if (optimized$objective > 1) {
                        warning(
                            paste0(
                                "Returned number of genes for dataset ",
                                i,
                                " differs from requested by ",
                                optimized$objective,
                                ". Lower tol or alpha.thresh for better results."
                            )
                        )
                    }
                }

                genes.new <-
                    names(gene_expr_var)[which(
                        gene_expr_var / nolan_constant > genemeanupper &
                            log10(gene_expr_var) > basegenelower + var.thresh[i]
                    )]

                if (do.plot) {
                    graphics::plot(
                        log10(gene_expr_mean),
                        log10(gene_expr_var),
                        cex = cex.use,
                        xlab = 'Gene Expression Mean (log10)',
                        ylab = 'Gene Expression Variance (log10)'
                    )

                    graphics::points(
                        log10(gene_expr_mean[genes.new]),
                        log10(gene_expr_var[genes.new]),
                        cex = cex.use,
                        col = "green"
                    )
                    graphics::abline(log10(nolan_constant), 1, col = "purple")

                    legend(
                        "bottomright",
                        paste0("Selected genes: ", length(genes.new)),
                        pch = 20,
                        col = "green"
                    )
                    graphics::title(main = names(object@raw.data)[i])
                }
                if (combine == "union") {
                    genes.use <- union(genes.use, genes.new)
                }
                if (combine == "intersection") {
                    if (length(genes.use) == 0) {
                        genes.use <- genes.new
                    }
                    genes.use <- intersect(genes.use, genes.new)
                }
            }

            for (i in 1:length(object@raw.data)) {
                genes.use <-
                    genes.use[genes.use %in% rownames(object@raw.data[[i]])]
            }

            if (length(genes.use) == 0) {
                warning(
                    "No genes were selected; lower var.thresh values or choose 'union' for combine parameter",
                    immediate. = TRUE
                )
            }
            object@var.genes <- genes.use
        }
        # Only for unshared Features
        if (unshared == T) {
            ind.thresh = c()
            # If only one threshold is provided, apply to all unshared datasets
            if (length(unshared.thresh == 1)) {
                ind.thresh = rep(unshared.thresh, length(object@raw.data))
            }    else{
                # If thresholds are provided for every dataset, use the respective threshold for each datatset
                if (length(unshared.thresh) != length(unshared.datasets)) {
                    warning(
                        "The number of thresholds does not match the number of datasets; Please provide either a single threshold value or a value for each unshared dataset.",
                        immediate. = T
                    )
                }
                names(unshared.thresh) = unshared.datasets
                for (i in unshared.datasets) {
                    ind.thresh[[i]] = unshared.thresh$i
                }
            }
            unshared.feats <- c()

            for (i in 1:length(object@raw.data)) {
                unshared.feats[i] <- list(NULL)
            }

            #construct a list of shared features
            shared_names = rownames(object@raw.data[[1]])
            for (matrix in 2:length(object@raw.data)) {
                shared_names = subset(shared_names,
                                      shared_names %in% rownames(object@raw.data[[i]]))
            }

            for (i in unshared.datasets) {
                unshared.use <- c()
                #Provides normalized subset of unshared features
                normalized_unshared = object@norm.data[[i]][!rownames(object@norm.data[[i]]) %in% shared_names, ]
                #Selects top variable features
                genes.unshared <- c()
                trx_per_cell <- colSums(object@raw.data[[i]])
                # Each gene's mean expression level (across all cells)
                gene_expr_mean <- rowMeansFast(normalized_unshared)
                # Each gene's expression variance (across all cells)
                gene_expr_var <-
                    rowVarsFast(normalized_unshared, gene_expr_mean)
                names(gene_expr_mean) <-
                    names(gene_expr_var) <- rownames(normalized_unshared)
                nolan_constant <- mean((1 / trx_per_cell))
                alphathresh.corrected <-
                    alpha.thresh / nrow(object@raw.data[[i]])
                genemeanupper <-
                    gene_expr_mean + qnorm(1 - alphathresh.corrected / 2) *
                    sqrt(gene_expr_mean * nolan_constant / ncol(object@raw.data[[i]]))
                basegenelower <- log10(gene_expr_mean * nolan_constant)
                genes.unshared <-
                    names(gene_expr_var)[which(
                        gene_expr_var / nolan_constant > genemeanupper &
                            log10(gene_expr_var) > basegenelower + ind.thresh[[i]]
                    )]
                if (length(genes.unshared) == 0) {
                    warning(
                        'Dataset ',
                        i ,
                        ' does not contain any unshared features. Please remove this dataset from the unshared.datasets list and rerun the function',
                        immediate. = TRUE
                    )
                }
                if (length(genes.unshared != 0)) {
                    unshared.feats[[i]] <- c(genes.unshared)
                }
            }
            names(unshared.feats) <- names(object@raw.data)
            object@var.unshared.features <- unshared.feats
            for (i in unshared.datasets) {
                print(
                    paste0(
                        "Selected ",
                        length(unshared.feats[[i]]),
                        " unshared features from ",
                        names(unshared.feats)[i],
                        " Dataset"
                    )
                )
            }
        }
        return(object)
    }

#' Scale genes by root-mean-square across cells
#'
#' This function scales normalized gene expression data after variable genes have been selected.
#' Note that the data is not mean-centered before scaling because expression values must remain
#' positive (NMF only accepts positive values). It also removes cells which do not have any
#' expression across the genes selected, by default.
#'
#' @param object \code{liger} object. Should call normalize and selectGenes before calling.
#' @param remove.missing Whether to remove cells from scale.data with no gene expression
#'   (default TRUE).
#' @param chunk size of chunks in hdf5 file. (default 1000)
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return \code{liger} object with scale.data slot set.
#'
#' @export
#' @examples
#' \dontrun{
#' # Given datasets Y and Z
#' ligerex <- createLiger(list(y_set = Y, z_set = Z))
#' ligerex <- normalize(ligerex)
#' # use default selectGenes settings (var.thresh = 0.1)
#' ligerex <- selectGenes(ligerex)
#' ligerex <- scaleNotCenter(ligerex)
#' }

scaleNotCenter <-
    function(object,
             remove.missing = TRUE,
             chunk = 1000,
             verbose = TRUE) {
        if (class(object@raw.data[[1]])[1] == "H5File") {
            hdf5_files = names(object@raw.data)
            vargenes = object@var.genes
            for (i in 1:length(hdf5_files)) {
                if (verbose) {
                    message(hdf5_files[i])
                }
                chunk_size = chunk

                if (object@h5file.info[[i]][["format.type"]] == "AnnData") {
                    genes = object@raw.data[[i]][["raw.var"]][]$index
                } else {
                    genes = object@h5file.info[[i]][["genes"]][]
                }
                num_cells = object@h5file.info[[i]][["barcodes"]]$dims
                num_genes = length(genes)
                num_entries = object@h5file.info[[i]][["data"]]$dims

                prev_end_col = 1
                prev_end_data = 1
                prev_end_ind = 0
                gene_vars = rep(0, num_genes)
                gene_means = object@raw.data[[i]][["gene_means"]][1:num_genes]
                gene_sum_sq = object@raw.data[[i]][["gene_sum_sq"]][1:num_genes]

                gene_inds = which(genes %in% vargenes)
                gene_root_mean_sum_sq = sqrt(gene_sum_sq / (num_cells - 1))
                safe_h5_create(
                    object = object,
                    idx = i,
                    dataset_name = "scale.data",
                    dims = c(length(vargenes), num_cells),
                    mode = h5types$double,
                    chunk_size = c(length(vargenes), chunk_size)
                )
                num_chunks = ceiling(num_cells / chunk_size)
                if (verbose) {
                    pb = txtProgressBar(0, num_chunks, style = 3)
                }
                ind = 0
                while (prev_end_col < num_cells) {
                    ind = ind + 1
                    if (num_cells - prev_end_col < chunk_size) {
                        chunk_size = num_cells - prev_end_col + 1
                    }
                    start_inds = object@h5file.info[[i]][["indptr"]][prev_end_col:(prev_end_col +
                                                                                       chunk_size)]
                    row_inds = object@h5file.info[[i]][["indices"]][(prev_end_ind +
                                                                         1):(tail(start_inds, 1))]
                    counts = object@norm.data[[i]][(prev_end_ind + 1):(tail(start_inds, 1))]
                    scaled = sparseMatrix(
                        i = row_inds[1:length(counts)] + 1,
                        p = start_inds[1:(chunk_size + 1)] - prev_end_ind,
                        x = counts,
                        dims = c(num_genes, chunk_size)
                    )
                    scaled = scaled[gene_inds,]
                    scaled = as.matrix(scaled)
                    root_mean_sum_sq = gene_root_mean_sum_sq[gene_inds]
                    scaled = sweep(scaled, 1, root_mean_sum_sq, "/")
                    rownames(scaled) = genes[gene_inds]
                    scaled = scaled[vargenes,]
                    scaled[is.na(scaled)] = 0
                    scaled[scaled == Inf] = 0
                    object@raw.data[[i]][["scale.data"]][1:length(vargenes), prev_end_col:(prev_end_col +
                                                                                               chunk_size - 1)] = scaled
                    num_read = length(counts)
                    prev_end_col = prev_end_col + chunk_size
                    prev_end_data = prev_end_data + num_read
                    prev_end_ind = tail(start_inds, 1)
                    if (verbose) {
                        setTxtProgressBar(pb, ind)
                    }
                }
                object@scale.data[[i]] = object@raw.data[[i]][["scale.data"]]
                if (verbose) {
                    setTxtProgressBar(pb, num_chunks)
                    cat("\n")
                }
            }
            names(object@scale.data) <- names(object@raw.data)
        } else {
            object@scale.data <-
                lapply(1:length(object@norm.data), function(i) {
                    scaleNotCenterFast(t(object@norm.data[[i]][object@var.genes,]))
                })
            # TODO: Preserve sparseness later on (convert inside optimizeALS)
            object@scale.data <- lapply(object@scale.data, function(x) {
                as.matrix(x)
            })

            names(object@scale.data) <- names(object@norm.data)
            for (i in 1:length(object@scale.data)) {
                object@scale.data[[i]][is.na(object@scale.data[[i]])] <- 0
                rownames(object@scale.data[[i]]) <-
                    colnames(object@raw.data[[i]])
                colnames(object@scale.data[[i]]) <- object@var.genes
            }
            # may want to remove such cells before scaling -- should not matter for large datasets?
        }
        #Scale unshared features
        if (length(object@var.unshared.features) != 0) {
            for (i in 1:length(object@raw.data)) {
                if (!is.null(object@var.unshared.features[[i]])) {
                    if (class(object@raw.data[[i]])[1] == "dgTMatrix" |
                        class(object@raw.data[[i]])[1] == "dgCMatrix") {
                        object@scale.unshared.data[[i]] <-
                            scaleNotCenterFast(t(object@norm.data[[i]][object@var.unshared.features[[i]],]))
                        object@scale.unshared.data[[i]] <-
                            as.matrix(object@scale.unshared.data[[i]])
                    } else {
                        object@scale.unshared.data[[i]] <-
                            scale(t(object@norm.data[[i]][object@var.unshared.features[[i]],]),
                                  center = F,
                                  scale = T)
                    }
                    #names(object@scale.unshared.data) <- names(object@norm.data)
                    object@scale.unshared.data[[i]][is.na(object@scale.unshared.data[[i]])] <-
                        0
                    rownames(object@scale.unshared.data[[i]]) <-
                        colnames(object@raw.data[[i]])
                    colnames(object@scale.unshared.data[[i]]) <-
                        object@var.unshared.features[[i]]
                    #Remove cells that were deemed missing for the shared features
                    object@scale.unshared.data[[i]] <-
                        t(object@scale.unshared.data[[i]][rownames(object@scale.data[[i]]), ])
                } else{
                    object@scale.unshared.data[i] <- NA
                }
            }
            names(object@scale.unshared.data) <- names(object@norm.data)
        }
        return(object)
    }

#' Remove cells/genes with no expression across any genes/cells
#'
#' Removes cells/genes from chosen slot with no expression in any genes or cells respectively.
#'
#' @param object \code{liger} object (scale.data or norm.data must be set).
#' @param slot.use The data slot to filter (takes "raw.data" and "scale.data") (default "raw.data").
#' @param use.cols Treat each column as a cell (default TRUE).
#' @param verbose Print messages (TRUE by default)
#'
#' @return \code{liger} object with modified raw.data (or chosen slot) (dataset names preserved).
#'
#' @export
#' @examples
#' \dontrun{
#' # liger object: ligerex
#' ligerex <- removeMissingObs(ligerex)
#' }

removeMissingObs <-
    function(object,
             slot.use = "raw.data",
             use.cols = TRUE,
             verbose = TRUE) {
        filter.data <- slot(object, slot.use)
        removed <-
            ifelse(((
                slot.use %in% c("raw.data", "norm.data")
            ) & (use.cols == TRUE)) |
                ((slot.use == "scale.data") &
                     (use.cols == FALSE)) ,
            yes = "cells", no = "genes")
        expressed <-
            ifelse(removed == "cells", yes = " any genes", no = "")
        filter.data <- lapply(seq_along(filter.data), function(x) {
            if (use.cols) {
                missing <- which(colSums(filter.data[[x]]) == 0)
            } else {
                missing <- which(rowSums(filter.data[[x]]) == 0)
            }
            if (length(missing) > 0) {
                if (verbose) {
                    message(
                        "Removing ",
                        length(missing),
                        " ",
                        removed,
                        " not expressing",
                        expressed,
                        " in ",
                        names(object@raw.data)[x],
                        "."
                    )
                }
                if (use.cols) {
                    if (length(missing) < 25) {
                        if (verbose) {
                            message(writeLines(colnames(filter.data[[x]])[missing]))
                        }
                    }
                    subset <- filter.data[[x]][, -missing]
                } else {
                    if (length(missing) < 25) {
                        if (verbose) {
                            message(writeLines(rownames(filter.data[[x]])[missing]))
                        }
                    }
                    subset <- filter.data[[x]][-missing,]
                }
            } else {
                subset <- filter.data[[x]]
            }
            subset
        })
        names(filter.data) <- names(object@raw.data)
        slot(object, slot.use) <- filter.data
        return(object)
    }

#helper function for readSubset
#Samples cell barcodes from specified datasets
#balance=NULL (default) means that max_cells are sampled from among all cells.
#balance="cluster" samples up to max_cells from each cluster in each dataset
#balance="dataset" samples up to max_cells from each dataset
#datasets.use uses only the specified datasets for sampling. Default is NULL (all datasets)
#rand.seed for reproducibility (default 1).
#verbose for printing messages
#Returns: vector of cell barcodes
downsample <-
    function(object,
             balance = NULL,
             max_cells = 1000,
             datasets.use = NULL,
             seed = 1,
             verbose = TRUE)
    {
        set.seed(seed)
        if (is.null(datasets.use))
        {
            datasets.use = names(object@H)
            if (verbose) {
                message(datasets.use)
            }
        }
        inds = c()
        inds_ds = list()
        if (is.null(balance))
        {
            for (ds in 1:length(datasets.use))
            {
                inds = c(inds, rownames(object@H[[ds]]))
            }
            num_to_samp = min(max_cells, length(inds))
            inds = sample(inds, num_to_samp)
            for (ds in 1:length(datasets.use))
            {
                inds_ds[[ds]] = intersect(inds, rownames(object@H[[ds]]))
            }
        }
        else if (balance == "dataset")
        {
            for (ds in 1:length(datasets.use))
            {
                num_to_samp = min(max_cells, nrow(object@H[[ds]]))
                inds_ds[[ds]] = rownames(object@H[[ds]])[sample(1:nrow(object@H[[ds]]), num_to_samp)]
            }
        }
        else
            #balance clusters
        {
            if (nrow(object@cell.data) == 0)
            {
                dataset <- unlist(lapply(seq_along(object@H), function(i) {
                    rep(names(object@H)[i], nrow(object@H[[i]]))
                }), use.names = FALSE)
                object@cell.data <- data.frame(dataset)
                rownames(object@cell.data) <- unlist(lapply(object@H,
                                                            function(x) {
                                                                rownames(x)
                                                            }), use.names = FALSE)
            }
            for (ds in 1:length(datasets.use))
            {
                for (i in levels(object@clusters))
                {
                    inds_to_samp = names(object@clusters)[object@clusters == i &
                                                              object@cell.data[["dataset"]] == ds]
                    num_to_samp = min(max_cells, length(inds_to_samp))
                    inds_ds[[ds]] = sample(inds_to_samp, num_to_samp)
                }
            }
        }
        return(inds_ds)
    }

#' Sample data for plotting
#'
#' This function samples raw/normalized/scaled data from on-disk HDF5 files for plotting.
#' This function assumes that the cell barcodes are unique across all datasets.
#'
#' @param object \code{liger} object. Should call normalize and selectGenes before calling.
#' @param slot.use Type of data for sampling (raw.data, norm.data(default), scale.data).
#' @param balance Type of sampling. NULL means that max_cells are sampled from among all cells;
#'                balance="dataset" samples up to max_cells from each dataset;
#'                balance="cluster" samples up to max_cells from each cluster.
#' @param chunk is the max number of cells at a time to read from disk (default 1000).
#' @param max.cells Total number of cell to sample (default 5000).
#' @param rand.seed  (default 1).
#' @param datasets.use uses only the specified datasets for sampling. Default is NULL (all datasets)
#' @param genes.use samples from only the specified genes. Default is NULL (all genes)
#' @param rand.seed for reproducibility (default 1).
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return \code{liger} object with sample.data slot set.
#'
#' @export
#' @examples
#' \dontrun{
#' # Only for online liger object (based on HDF5 files)
#' # Example: sample a total amount of 5000 cells from norm.data for downstream analysis
#' ligerex <- readSubset(ligerex, slot.use = "norm.data", max.cells = 5000)
#' }

readSubset <- function(object,
                       slot.use = "norm.data",
                       balance = NULL,
                       max.cells = 1000,
                       chunk = 1000,
                       datasets.use = NULL,
                       genes.use = NULL,
                       rand.seed = 1,
                       verbose = TRUE) {
    if (class(object@raw.data[[1]])[1] == "H5File") {
        if (verbose) {
            message("Start sampling")
        }
        if (is.null(datasets.use))
        {
            datasets.use = names(object@H)
        }
        cell_inds = downsample(
            object,
            balance = balance,
            max_cells = max.cells,
            datasets.use = datasets.use,
            seed = rand.seed,
            verbose = verbose
        )

        hdf5_files = names(object@raw.data)
        #vargenes = object@var.genes

        # find the intersect of genes from each input datasets
        genes = c()
        if (slot.use != "scale.data") {
            for (i in 1:length(hdf5_files)) {
                if (object@h5file.info[[i]][["format.type"]] == "AnnData") {
                    genes_i = object@h5file.info[[i]][["genes"]][]$index
                } else {
                    genes_i = object@h5file.info[[i]][["genes"]][]
                }
                if (i == 1)
                    genes = genes_i
                else
                    genes = intersect(genes, genes_i)
            }
        } else {
            genes = object@var.genes
        }

        if (is.null(genes.use))
        {
            genes.use = genes
        }

        for (i in 1:length(hdf5_files)) {
            if (verbose) {
                message(hdf5_files[i])
            }
            if (slot.use == "scale.data") {
                data.subset = c()
            } else {
                data.subset = Matrix(
                    nrow = length(genes.use),
                    ncol = 0,
                    sparse = TRUE
                )
            }
            chunk_size = chunk
            if (object@h5file.info[[i]][["format.type"]] == "AnnData") {
                barcodes = object@h5file.info[[i]][["barcodes"]][]$cell
                genes = object@h5file.info[[i]][["genes"]][]$index
            } else {
                barcodes = object@h5file.info[[i]][["barcodes"]][]
                genes = object@h5file.info[[i]][["genes"]][]
            }
            num_cells = length(barcodes)
            num_genes = length(genes)

            prev_end_col = 1
            prev_end_data = 1
            prev_end_ind = 0


            #gene_inds = which(genes %in% vargenes)

            num_chunks = ceiling(num_cells / chunk_size)
            if (verbose) {
                pb = txtProgressBar(0, num_chunks, style = 3)
            }
            ind = 0

            while (prev_end_col < num_cells) {
                ind = ind + 1
                if (num_cells - prev_end_col < chunk_size) {
                    chunk_size = num_cells - prev_end_col + 1
                }
                if (slot.use != "scale.data") {
                    start_inds = object@h5file.info[[i]][["indptr"]][prev_end_col:(prev_end_col +
                                                                                       chunk_size)]
                    row_inds = object@h5file.info[[i]][["indices"]][(prev_end_ind +
                                                                         1):(tail(start_inds, 1))]
                    if (slot.use == "raw.data")
                    {
                        counts = object@h5file.info[[i]][["data"]][(prev_end_ind + 1):(tail(start_inds, 1))]
                    }
                    if (slot.use == "norm.data")
                    {
                        counts = object@norm.data[[i]][(prev_end_ind + 1):(tail(start_inds, 1))]
                    }
                    one_chunk = sparseMatrix(
                        i = row_inds[1:length(counts)] + 1,
                        p = start_inds[1:(chunk_size + 1)] - prev_end_ind,
                        x = counts,
                        dims = c(num_genes, chunk_size)
                    )
                    rownames(one_chunk) = genes
                    colnames(one_chunk) = barcodes[(prev_end_col):(prev_end_col +
                                                                       chunk_size - 1)]
                    use_these = intersect(colnames(one_chunk), cell_inds[[i]])
                    one_chunk = one_chunk[genes.use, use_these]
                    data.subset = cbind(data.subset, one_chunk)

                    num_read = length(counts)
                    prev_end_col = prev_end_col + chunk_size
                    prev_end_data = prev_end_data + num_read
                    prev_end_ind = tail(start_inds, 1)
                    setTxtProgressBar(pb, ind)
                } else {
                    one_chunk = object@scale.data[[i]][, prev_end_col:(prev_end_col + chunk_size - 1)]
                    rownames(one_chunk) = object@var.genes
                    colnames(one_chunk) = barcodes[(prev_end_col):(prev_end_col +
                                                                       chunk_size - 1)]
                    use_these = intersect(colnames(one_chunk), cell_inds[[i]])
                    one_chunk = one_chunk[genes.use, use_these]
                    data.subset = cbind(data.subset, one_chunk)

                    prev_end_col = prev_end_col + chunk_size
                    if (verbose) {
                        setTxtProgressBar(pb, ind)
                    }
                }
                if (class(object@raw.data[[i]])[1] == "H5File") {
                    object@sample.data[[i]] = data.subset
                } else if (class(object@raw.data[[i]])[1] != "H5File" &
                           slot.use == "scale.data") {
                    object@sample.data[[i]] = t(data.subset)
                }
                object@h5file.info[[i]][["sample.data.type"]] = slot.use
            }
            if (verbose) {
                setTxtProgressBar(pb, num_chunks)
                cat("\n")
            }
        }
    } else {
        if (verbose) {
            message("Start sampling")
        }
        if (is.null(datasets.use))
        {
            datasets.use = names(object@H)
        }
        cell_inds = downsample(
            object,
            balance = balance,
            max_cells = max.cells,
            datasets.use = datasets.use,
            verbose = verbose
        )

        files = names(object@raw.data)
        # find the intersect of genes from each input datasets
        genes = c()
        for (i in 1:length(files)) {
            genes_i = rownames(object@raw.data[[i]])
            if (i == 1)
                genes = genes_i
            else
                genes = intersect(genes, genes_i)
        }
        if (is.null(genes.use))
        {
            genes.use = genes
        }
        if (verbose) {
            pb = txtProgressBar(0, length(files), style = 3)
        }
        for (i in 1:length(files)) {
            if (slot.use == "raw.data")
            {
                data.subset_i = object@raw.data[[i]][genes.use, cell_inds[[i]]]
            }
            if (slot.use == "norm.data")
            {
                data.subset_i = object@norm.data[[i]][genes.use, cell_inds[[i]]]
            }
            if (slot.use == "scale.data")
            {
                data.subset_i = t(object@scale.data[[i]][cell_inds[[i]], genes.use])
            }
            if (verbose) {
                setTxtProgressBar(pb, i)
            }
        }
        if (verbose) {
            cat("\n")
        }
        object@sample.data[[i]] = data.subset_i
        object@h5file.info[[i]][["sample.data.type"]] = slot.use
    }
    names(object@sample.data) = names(object@raw.data)
    return(object)
}
