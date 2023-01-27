#' Louvain algorithm for community detection
#'
#' @description
#' After quantile normalization, users can additionally run the Louvain algorithm
#' for community detection, which is widely used in single-cell analysis and excels at merging
#' small clusters into broad cell classes.
#'
#' @param object \code{liger} object. Should run quantile_norm before calling.
#' @param k The maximum number of nearest neighbours to compute. (default 20)
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want
#' to obtain a larger (smaller) number of communities. (default 1.0)
#' @param prune Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
#' prune everything). (default 1/15)
#' @param eps The error bound of the nearest neighbor search. (default 0.1)
#' @param nRandomStarts Number of random starts. (default 10)
#' @param nIterations Maximal number of iterations per random start. (default 100)
#' @param random.seed Seed of the random number generator. (default 1)
#' @param verbose Print messages (TRUE by default)
#' @param dims.use Indices of factors to use for Louvain clustering (default 1:ncol(H[[1]])).
#'
#' @return \code{liger} object with refined 'clusters' slot set.
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' ligerex <- louvainCluster(ligerex, resulotion = 0.3)
#' }
#'

louvainClusterOld <-
    function(object,
             resolution = 1.0,
             k = 20,
             prune = 1 / 15,
             eps = 0.1,
             nRandomStarts = 10,
             nIterations = 100,
             random.seed = 1,
             verbose = TRUE,
             dims.use = NULL) {
        output_path <- paste0('edge_', sub('\\s', '_', Sys.time()), '.txt')
        output_path = sub(":", "_", output_path)
        output_path = sub(":", "_", output_path)

        if (is.null(dims.use)) {
            use_these_factors <- 1:ncol(object@H[[1]])
        } else {
            use_these_factors <- dims.use
        }

        if (dim(object@H.norm)[1] == 0) {
            if (verbose) {
                message("Louvain Clustering on unnormalized cell factor loadings.")
            }
            knn <-
                RANN::nn2(Reduce(rbind, object@H)[, use_these_factors],
                          k = k,
                          eps = eps)
        } else {
            if (verbose) {
                message("Louvain Clustering on quantile normalized cell factor loadings.")
            }
            knn <-
                RANN::nn2(object@H.norm[, use_these_factors], k = k, eps = eps)
        }
        snn <- ComputeSNN(knn$nn.idx, prune = prune)
        WriteEdgeFile(snn, output_path, display_progress = FALSE)
        clusts <- RunModularityClusteringCpp(
            snn,
            modularityFunction = 1,
            resolution = resolution,
            nRandomStarts = nRandomStarts,
            nIterations = nIterations,
            algorithm = 1,
            randomSeed = random.seed,
            printOutput = FALSE,
            edgefilename = output_path
        )
        names(clusts) = rownames(object@cell.data)
        rownames(snn) = rownames(object@cell.data)
        colnames(snn) = rownames(object@cell.data)
        clusts <-
            GroupSingletons(ids = clusts,
                            SNN = snn,
                            verbose = FALSE)
        object@clusters = as.factor(clusts)
        unlink(output_path)
        return(object)
    }

# Group single cells that make up their own cluster in with the cluster they are
# most connected to. (Adopted from Seurat v3)
#
# @param ids Named vector of cluster ids
# @param SNN SNN graph used in clustering
# @param group.singletons Group singletons into nearest cluster (TRUE by default). If FALSE, assign all singletons to
# @param verbose Print message
# a "singleton" group
#
# @return Returns updated cluster assignment with all singletons merged with most connected cluster
#
GroupSingletons <-
    function(ids,
             SNN,
             group.singletons = TRUE,
             verbose = FALSE) {
        # identify singletons
        singletons <- c()
        singletons <- names(x = which(x = table(ids) == 1))
        singletons <- intersect(x = unique(x = ids), singletons)
        if (!group.singletons) {
            ids[which(ids %in% singletons)] <- "singleton"
            return(ids)
        }
        # calculate connectivity of singletons to other clusters, add singleton
        # to cluster it is most connected to
        cluster_names <- as.character(x = unique(x = ids))
        cluster_names <- setdiff(x = cluster_names, y = singletons)
        connectivity <-
            vector(mode = "numeric", length = length(x = cluster_names))
        names(x = connectivity) <- cluster_names
        new.ids <- ids
        for (i in singletons) {
            i.cells <- names(which(ids == i))
            for (j in cluster_names) {
                j.cells <- names(which(ids == j))
                subSNN <- SNN[i.cells, j.cells]
                set.seed(1) # to match previous behavior, random seed being set in WhichCells
                if (is.object(x = subSNN)) {
                    connectivity[j] <-
                        sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
                } else {
                    connectivity[j] <- mean(x = subSNN)
                }
            }
            m <- max(connectivity, na.rm = TRUE)
            mi <- which(x = connectivity == m, arr.ind = TRUE)
            closest_cluster <-
                sample(x = names(x = connectivity[mi]), 1)
            ids[i.cells] <- closest_cluster
        }
        if (length(x = singletons) > 0 && verbose) {
            message(paste(
                length(x = singletons),
                "singletons identified.",
                length(x = unique(x = ids)),
                "final clusters."
            ))
        }
        return(ids)
    }


#' Impute the query cell expression matrix
#'
#' Impute query features from a reference dataset using KNN.
#'
#' @param object \code{liger} object.
#' @param knn_k The maximum number of nearest neighbors to search. (default 20)
#' @param reference Dataset containing values to impute into query dataset(s).
#' @param queries Dataset to be augmented by imputation. If not specified, will pass in all datasets.
#' @param weight Whether to use KNN distances as weight matrix (default FALSE).
#' @param norm Whether normalize the imputed data with default parameters (default TRUE).
#' @param scale Whether scale but not center the imputed data with default parameters (default TRUE).
#' @param verbose Print progress bar/messages (TRUE by default)
#'
#' @return \code{liger} object with raw data in raw.data slot replaced by imputed data (genes by cells)
#'
#' @importFrom FNN get.knnx
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' # impute every dataset other than the reference dataset
#' ligerex <- imputeKNN(ligerex, reference = "y_set", weight = FALSE)
#' # impute only z_set dataset
#' ligerex <- imputeKNN(ligerex, reference = "y_set", queries = list("z_set"), knn_k = 50)
#' }

imputeKNN <-
    function(object,
             reference,
             queries,
             knn_k = 20,
             weight = TRUE,
             norm = TRUE,
             scale = FALSE,
             verbose = TRUE) {
        if (verbose) {
            cat(
                "NOTE: This function will discard the raw data previously stored in the liger object and",
                "replace the raw.data slot with the imputed data.\n\n"
            )
        }

        if (length(reference) > 1) {
            stop("Can only have ONE reference dataset")
        }
        if (missing(queries)) {
            # all datasets
            queries <- names(object@raw.data)
            queries <- as.list(queries[!queries %in% reference])
            if (verbose) {
                cat(
                    "Imputing ALL the datasets except the reference dataset\n",
                    "Reference dataset:\n",
                    paste("  ", reference, "\n"),
                    "Query datasets:\n",
                    paste("  ", as.character(queries), "\n")
                )
            }
        }
        else {
            # only given query datasets
            queries <- as.list(queries)
            if (reference %in% queries) {
                stop("Reference dataset CANNOT be inclued in the query datasets")
            }
            else {
                if (verbose) {
                    cat(
                        "Imputing given query datasets\n",
                        "Reference dataset:\n",
                        paste("  ", reference, "\n"),
                        "Query datasets:\n",
                        paste("  ", as.character(queries), "\n")
                    )
                }
            }
        }

        reference_cells <-
            colnames(object@raw.data[[reference]]) # cells by genes
        for (query in queries) {
            query_cells <- colnames(object@raw.data[[query]])

            # creating a (reference cell numbers X query cell numbers) weights matrix for knn weights and unit weights
            nn.k <-
                get.knnx(object@H.norm[reference_cells, ],
                         object@H.norm[query_cells, ],
                         k = knn_k,
                         algorithm = "CR")
            weights <-
                Matrix(
                    0,
                    nrow = ncol(object@raw.data[[reference]]),
                    ncol = nrow(nn.k$nn.index),
                    sparse = TRUE
                )
            if (weight == TRUE) {
                # for weighted situation
                # find nearest neighbors for query cell in normed ref datasets
                for (n in 1:nrow(nn.k$nn.index)) {
                    # record ref-query cell-cell distances
                    weights[nn.k$nn.index[n, ], n] <-
                        exp(-nn.k$nn.dist[n, ]) / sum(exp(-nn.k$nn.dist[n, ]))
                }
            }
            else{
                # for unweighted situation
                for (n in 1:nrow(nn.k$nn.index)) {
                    weights[nn.k$nn.index[n, ], n] <- 1 / knn_k # simply count the mean
                }
            }

            # (genes by ref cell num) multiply by the weight matrix (ref cell num by query cell num)
            imputed_vals <- object@raw.data[[reference]] %*% weights
            # assigning dimnames
            colnames(imputed_vals) <- query_cells
            rownames(imputed_vals) <-
                rownames(object@raw.data[[reference]])

            # formatiing the matrix
            if (class(object@raw.data[[reference]])[1] == "dgTMatrix" |
                class(object@raw.data[[reference]])[1] == "dgCMatrix") {
                imputed_vals <- as(imputed_vals, "dgCMatrix")
            } else {
                imputed_vals <- as.matrix(imputed_vals)
            }

            object@raw.data[[query]] <- imputed_vals
        }

        if (norm) {
            if (verbose) {
                cat('\nNormalizing data...\n')
            }
            object <- rliger::normalize(object)
        }
        if (scale) {
            if (verbose) {
                cat('Scaling (but not centering) data...')
            }
            object <- rliger::scaleNotCenter(object)
        }

        return(object)
    }

#' Perform Wilcoxon rank-sum test
#'
#' Perform Wilcoxon rank-sum tests on specified dataset using given method.
#'
#' @param object \code{liger} object.
#' @param data.use This selects which dataset(s) to use. (default 'all')
#' @param compare.method This indicates the metric of the test. Either 'clusters' or 'datasets'.
#'
#' @return A 10-columns data.frame with test results.
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object based on in-memory datasets), factorization complete
#' wilcox.results <- runWilcoxon(ligerex, compare.method = "cluster")
#' wilcox.results <- runWilcoxon(ligerex, compare.method = "datastes", data.use = c(1, 2))
#' # HDF5 input
#' # ligerex (liger object based on datasets in HDF5 format), factorization complete
#' # Need to sample cells before implementing Wilcoxon test
#' ligerex <- readSubset(ligerex, slot.use = "norm.data", max.cells = 1000)
#' de_genes <- runWilcoxon(ligerex, compare.method = "clusters")
#' }

runWilcoxon <- function(object, data.use = "all", compare.method) {
    # check parameter inputs
    if (missing(compare.method)) {
        stop("Parameter *compare.method* cannot be empty!")
    }
    if (compare.method != "datasets" &
        compare.method != "clusters") {
        stop("Parameter *compare.method* should be either *clusters* or *datasets*.")
    }
    if (compare.method == "datasets") {
        if (length(names(object@norm.data)) < 2) {
            stop("Should have at least TWO inputs to compare between datasets")
        }
        if (!missing(data.use) & length(data.use) < 2) {
            stop("Should have at least TWO inputs to compare between datasets")
        }
    }

    if (class(object@raw.data[[1]])[1] == "H5File") {
        if (is.null(object@h5file.info[[1]][["sample.data.type"]])) {
            message("Need to sample data before Wilcoxon test for HDF5 input.")
        } else {
            message("Running Wilcoxon test on ", object@h5file.info[[1]][["sample.data.type"]])
        }
    }

    ### create feature x sample matrix
    if (data.use[1] == "all" |
        length(data.use) > 1) {
        # at least two datasets
        if (data.use[1] == "all") {
            message("Performing Wilcoxon test on ALL datasets: ",
                    toString(names(object@norm.data)))
            sample.list <- attributes(object@norm.data)$names
        }
        else {
            message("Performing Wilcoxon test on GIVEN datasets: ",
                    toString(data.use))
            sample.list <- data.use
        }
        # get all shared genes of every datasets
        genes <-
            Reduce(intersect, lapply(sample.list, function(sample) {
                if (class(object@norm.data[[sample]])[[1]] == "dgCMatrix")
                {
                    return(object@norm.data[[sample]]@Dimnames[[1]])
                }
                else
                {
                    return(rownames(object@sample.data[[sample]]))
                }
            }))
        if (class(object@norm.data[[sample.list[1]]])[[1]] == "dgCMatrix") {
            feature_matrix <-
                Reduce(cbind, lapply(sample.list, function(sample) {
                    object@norm.data[[sample]][genes, ]
                })) # get feature matrix, shared genes as rows and all barcodes as columns
        } else {
            feature_matrix <- Reduce(cbind, object@sample.data[sample.list])
        }
        # get labels of clusters and datasets
        cell_source <-
            object@cell.data[["dataset"]] # from which dataset
        names(cell_source) <- names(object@clusters)
        cell_source <-
            cell_source[colnames(feature_matrix), drop = TRUE]
        clusters <-
            object@clusters[colnames(feature_matrix), drop = TRUE] # from which cluster
    } else {
        # for one dataset only
        message("Performing Wilcoxon test on GIVEN dataset: ", data.use)
        if (class(object@norm.data[[data.use]])[[1]] == "dgCMatrix") {
            feature_matrix <- object@norm.data[[data.use]]
            clusters <-
                object@clusters[object@norm.data[[data.use]]@Dimnames[[2]], drop = TRUE] # from which cluster
        } else {
            feature_matrix <- object@sample.data[[data.use]]
            clusters <-
                object@clusters[colnames(object@sample.data[[data.use]]), drop = TRUE] # from which cluster
        }
    }

    ### perform wilcoxon test
    if (compare.method == "clusters") {
        # compare between clusters across datasets
        len <- nrow(feature_matrix)
        if (len > 100000) {
            message("Calculating Large-scale Input...")
            results <-
                Reduce(rbind, lapply(suppressWarnings(split(
                    seq(len), seq(len / 100000)
                )), function(index) {
                    wilcoxauc(log(feature_matrix[index, ] + 1e-10), clusters)
                }))
        } else {
            results <- wilcoxauc(log(feature_matrix + 1e-10), clusters)
        }
    }

    if (compare.method == "datasets") {
        # compare between datasets within each cluster
        results <-
            Reduce(rbind, lapply(levels(clusters), function(cluster) {
                sub_barcodes <-
                    names(clusters[clusters == cluster]) # every barcode within this cluster
                sub_label <-
                    paste0(cluster, "-", cell_source[sub_barcodes]) # data source for each cell
                sub_matrix <- feature_matrix[, sub_barcodes]
                if (length(unique(cell_source[sub_barcodes])) == 1) {
                    # if cluster has only 1 data source
                    message("Note: Skip Cluster ",
                            cluster,
                            " since it has only ONE data source.")
                    return()
                }
                return(wilcoxauc(log(sub_matrix + 1e-10), sub_label))
            }))
    }
    return(results)
}

#' Linking genes to putative regulatory elements
#'
#' Evaluate the relationships between pairs of genes and peaks based on specified distance metric.
#'
#' @param gene_counts A gene expression matrix (genes by cells) of normalized counts.
#' This matrix has to share the same column names (cell barcodes) as the matrix passed to peak_counts
#' @param peak_counts A peak-level matrix (peaks by cells) of normalized accessibility values, such as the one resulting from imputeKNN.
#' This matrix must share the same column names (cell barcodes) as the matrix passed to gene_counts.
#' @param genes.list A list of the genes symbols to be tested. If not specified,
#' this function will use all the gene symbols from the matrix passed to gmat by default.
#' @param alpha Significance threshold for correlation p-value. Peak-gene correlations with p-values below
#' this threshold are considered significant. The default is 0.05.
#' @param dist This indicates the type of correlation to calculate -- one of “spearman” (default), "pearson", or "kendall".
#' @param path_to_coords Path to the gene coordinates file.
#' @param verbose Print messages (TRUE by default)
#'
#' @return a sparse matrix with peak names as rows and gene symbols as columns, with each element indicating the
#' correlation between peak i and gene j (or 0 if the gene and peak are not significantly linked).
#'
#' @importFrom utils read.csv2
#' @export
#' @examples
#' \dontrun{
#' # some gene counts matrix: gmat.small
#' # some peak counts matrix: pmat.small
#' regnet <- linkGenesAndPeaks(gmat.small, pmat.small, dist = "spearman",
#' alpha = 0.05, path_to_coords = 'some_path')
#' }

linkGenesAndPeaks <-
    function(gene_counts,
             peak_counts,
             genes.list = NULL,
             dist = "spearman",
             alpha = 0.05,
             path_to_coords,
             verbose = TRUE) {
        ## check dependency
        if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
            stop(
                "Package \"GenomicRanges\" needed for this function to work. Please install it by command:\n",
                "BiocManager::install('GenomicRanges')",
                call. = FALSE
            )
        }

        if (!requireNamespace("IRanges", quietly = TRUE)) {
            stop(
                "Package \"IRanges\" needed for this function to work. Please install it by command:\n",
                "BiocManager::install('IRanges')",
                call. = FALSE
            )
        }

        ### make Granges object for peaks
        peak.names <- strsplit(rownames(peak_counts), "[:-]")
        chrs <- Reduce(append, lapply(peak.names, function(peak) {
            peak[1]
        }))
        chrs.start <-
            Reduce(append, lapply(peak.names, function(peak) {
                peak[2]
            }))
        chrs.end <-
            Reduce(append, lapply(peak.names, function(peak) {
                peak[3]
            }))
        peaks.pos <- GenomicRanges::GRanges(
            seqnames = chrs,
            ranges = IRanges::IRanges(as.numeric(chrs.start), end = as.numeric(chrs.end))
        )

        ### make Granges object for genes
        gene.names <-
            read.csv2(
                path_to_coords,
                sep = "\t",
                header = FALSE,
                stringsAsFactors = FALSE
            )
        gene.names <- gene.names[complete.cases(gene.names), ]
        genes.coords <- GenomicRanges::GRanges(
            seqnames = gene.names$V1,
            ranges = IRanges::IRanges(as.numeric(gene.names$V2), end = as.numeric(gene.names$V3))
        )
        names(genes.coords) <- gene.names$V4

        ### construct regnet
        gene_counts <- t(gene_counts) # cell x genes
        peak_counts <- t(peak_counts) # cell x genes

        # find overlap peaks for each gene
        if (missing(genes.list)) {
            genes.list <- colnames(gene_counts)
        }
        missing_genes <- !genes.list %in% names(genes.coords)
        if (sum(missing_genes) != 0 && verbose) {
            message("Removing ",
                    sum(missing_genes),
                    " genes not found in given gene coordinates...")
        }
        genes.list <- genes.list[!missing_genes]
        genes.coords <- genes.coords[genes.list]

        if (verbose) {
            message("Calculating correlation for gene-peak pairs...")
        }
        each.len <- 0
        # assign('each.len', 0, envir = globalenv())

        elements <- lapply(seq(length(genes.list)), function(pos) {
            gene.use <- genes.list[pos]
            # re-scale the window for each gene
            gene.loci <-
                GenomicRanges::trim(suppressWarnings(
                    GenomicRanges::promoters(
                        GenomicRanges::resize(genes.coords[gene.use],
                                              width = 1, fix = "start"),
                        upstream = 500000,
                        downstream = 500000
                    )
                ))
            peaks.use <-
                S4Vectors::queryHits(GenomicRanges::findOverlaps(peaks.pos, gene.loci))
            if ((x <-
                 length(peaks.use)) == 0L) {
                # if no peaks in window, skip this iteration
                return(list(NULL, as.numeric(each.len), NULL))
            }
            ### compute correlation and p-adj for genes and peaks ###
            res <- suppressWarnings(
                psych::corr.test(
                    x = gene_counts[, gene.use],
                    y = as.matrix(peak_counts[, peaks.use]),
                    method = dist,
                    adjust = "holm",
                    ci = FALSE,
                    use = "complete"
                )
            )
            pick <- res[["p"]] < alpha # filter by p-value
            pick[is.na(pick)] <- FALSE

            if (sum(pick) == 0) {
                # if no peaks are important, skip this iteration
                return(list(NULL, as.numeric(each.len), NULL))
            }
            else {
                res.corr <- as.numeric(res[["r"]][pick])
                peaks.use <- peaks.use[pick]
            }
            # each.len <<- each.len + length(peaks.use)
            assign('each.len', each.len + length(peaks.use), envir = parent.frame(2))
            return(list(
                as.numeric(peaks.use),
                as.numeric(each.len),
                res.corr
            ))
        })

        i_index <- Reduce(append, lapply(elements, function(ele) {
            ele[[1]]
        }))
        p_index <-
            c(0, Reduce(append, lapply(elements, function(ele) {
                ele[[2]]
            })))
        value_list <-
            Reduce(append, lapply(elements, function(ele) {
                ele[[3]]
            }))

        # make final sparse matrix
        regnet <- sparseMatrix(
            i = i_index,
            p = p_index,
            x = value_list,
            dims = c(ncol(peak_counts), length(genes.list)),
            dimnames = list(colnames(peak_counts), genes.list)
        )

        return(regnet)
    }


#' Export predicted gene-pair interaction
#'
#' Export the predicted gene-pair interactions calculated by upstream function 'linkGenesAndPeaks'
#'  into an Interact Track file which is compatible with UCSC Genome Browser.
#'
#' @param corr.mat A sparse matrix with peak names as rows and gene symbols as columns.
#' @param genes.list A list of the genes symbols to be tested. If not specified,
#' this function will use all the gene symbols from the matrix passed to gmat by default.
#' @param output_path Path in which the output file will be stored.
#' @param path_to_coords Path to the gene coordinates file.
#'
#' @return An Interact Track file stored in the specified path.
#'
#' @importFrom stats complete.cases
#' @importFrom utils write.table
#'
#' @export
#' @examples
#' \dontrun{
#' # some gene-peak correlation matrix: regent
#' makeInteractTrack(regnet, path_to_coords = 'some_path_to_gene_coordinates/hg19_genes.bed')
#' }

makeInteractTrack <-
    function(corr.mat,
             genes.list,
             output_path,
             path_to_coords) {
        # get genomic coordinates
        if (missing(path_to_coords)) {
            stop("Parameter 'path_to_coords' cannot be empty.")
        }

        ### make Granges object for genes
        genes.coords <- read.csv2(
            path_to_coords,
            sep = "\t",
            header = FALSE,
            colClasses =
                c(
                    "character",
                    "integer",
                    "integer",
                    "character",
                    "NULL",
                    "NULL"
                )
        )
        genes.coords <-
            genes.coords[complete.cases(genes.coords$V4), ]
        rownames(genes.coords) <- genes.coords[, 4]

        # split peak names into chrom and coordinates
        peak.names <- strsplit(rownames(corr.mat), "[:-]")
        chrs <- Reduce(append, lapply(peak.names, function(peak) {
            peak[1]
        }))
        chrs.start <-
            as.numeric(Reduce(append, lapply(peak.names, function(peak) {
                peak[2]
            })))
        chrs.end <-
            as.numeric(Reduce(append, lapply(peak.names, function(peak) {
                peak[3]
            })))

        # check genes.list
        if (missing(genes.list)) {
            genes.list <- colnames(corr.mat)
        }

        # check output_path
        if (missing(output_path)) {
            output_path <- getwd()
        }

        output_path <- paste0(output_path, "/Interact_Track.bed")
        track.doc <-
            paste0(
                'track type=interact name="Interaction Track" description="Gene-Peaks Links"',
                ' interactDirectional=true maxHeightPixels=200:100:50 visibility=full'
            )
        write(track.doc, file = output_path)

        genes_not_existed <- 0
        filtered_genes <- 0

        for (gene in genes.list) {
            if (!gene %in% colnames(corr.mat)) {
                # if gene not in the corr.mat
                genes_not_existed <- genes_not_existed + 1
                next
            }
            peaks.sel <- which(corr.mat[, gene] != 0)
            if (sum(peaks.sel) == 0) {
                filtered_genes <- filtered_genes + 1
                next
            }

            track <- data.frame(
                chrom = chrs[peaks.sel],
                chromStart = chrs.start[peaks.sel],
                chromEnd = chrs.end[peaks.sel],
                name = paste0(gene, "/", rownames(corr.mat)[peaks.sel]),
                score = 0,
                value = as.numeric(corr.mat[peaks.sel, gene]),
                exp = ".",
                color = 5,
                sourceChrom = chrs[peaks.sel],
                sourceStart = chrs.start[peaks.sel],
                sourceEnd = chrs.start[peaks.sel] + 1,
                sourceName = ".",
                sourceStrand = ".",
                targetChrom = genes.coords[gene, 1],
                targetStart = genes.coords[gene, 2],
                targetEnd = genes.coords[gene, 2] + 1,
                targetName = gene,
                targetStrand = "."
            )
            write.table(
                track,
                file = output_path,
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

        message("A total of ",
                genes_not_existed,
                " genes do not exist in input matrix.")
        message("A total of ",
                filtered_genes,
                " genes do not have significant correlated peaks.")
        message("The Interaction Track is stored in Path: ", output_path)
    }


#' Analyze biological interpretations of metagene
#'
#' Identify the biological pathways (gene sets from Reactome) that each metagene (factor) might belongs to.
#'
#' @param object \code{liger} object.
#' @param gene_sets A list of the Reactome gene sets names to be tested. If not specified,
#' this function will use all the gene sets from the Reactome by default
#' @param mat_w This indicates whether to use the shared factor loadings 'W' (default TRUE)
#' @param mat_v This indicates which V matrix to be added to the analysis. It can be a numeric number or a list
#' of the numerics.
#' @param custom_gene_sets A named list of character vectors of entrez gene ids. If not specified,
#' this function will use all the gene symbols from the input matrix by default
#'
#' @return A list of matrices with GSEA analysis for each factor
#'
#' @importFrom methods .hasSlot
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' wilcox.results <- runGSEA(ligerex)
#' wilcox.results <- runGSEA(ligerex, mat_v = c(1, 2))
#' }
#'
runGSEA <-
    function(object,
             gene_sets = c(),
             mat_w = TRUE,
             mat_v = 0,
             custom_gene_sets = c()) {
        if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
            stop(
                "Package \"org.Hs.eg.db\" needed for this function to work. Please install it by command:\n",
                "BiocManager::install('org.Hs.eg.db')",
                call. = FALSE
            )
        }

        if (!requireNamespace("reactome.db", quietly = TRUE)) {
            stop(
                "Package \"reactome.db\" needed for this function to work. Please install it by command:\n",
                "BiocManager::install('reactome.db')",
                call. = FALSE
            )
        }

        if (!requireNamespace("fgsea", quietly = TRUE)) {
            stop(
                "Package \"fgsea\" needed for this function to work. Please install it by command:\n",
                "BiocManager::install('fgsea')",
                call. = FALSE
            )
        }

        if (length(mat_v) > length(object@V)) {
            stop("The gene loading input is invalid.", call. = FALSE)
        }

        if (!.hasSlot(object, "W") | !.hasSlot(object, "V")) {
            stop("There is no W or V matrix. Please do iNMF first.", call. = FALSE)
        }

        if (mat_w) {
            gene_loadings <- object@W
            if (mat_v) {
                gene_loadings <-
                    gene_loadings + Reduce("+", lapply(mat_v, function(v) {
                        object@V[[v]]
                    }))
            }
        } else {
            gene_loadings <- Reduce("+", lapply(mat_v, function(v) {
                object@V[[v]]
            }))
        }

        gene_ranks <-
            t(apply(gene_loadings, MARGIN = 1, function(x) {
                rank(x)
            }))

        colnames(gene_ranks) <-
            sapply(colnames(gene_ranks), toupper)
        gene_id <-
            as.character(
                AnnotationDbi::mapIds(
                    org.Hs.eg.db::org.Hs.eg.db,
                    colnames(gene_ranks),
                    "ENTREZID",
                    "SYMBOL"
                )
            )
        colnames(gene_ranks) <- gene_id
        gene_ranks <- gene_ranks[, !is.na(colnames(gene_ranks))]
        if (inherits((custom_gene_sets)[1], "tbl_df")) {
            pathways <-
                split(custom_gene_sets,
                      x = custom_gene_sets$entrez_gene,
                      f = custom_gene_sets$gs_name)
            pathways <- lapply(pathways, function(x) {
                as.character(x)
            })
        } else if (length(custom_gene_sets)) {
            pathways <- custom_gene_sets
        } else {
            pathways <- fgsea::reactomePathways(colnames(gene_ranks))
            if (length(gene_sets)) {
                pathways <- pathways[intersect(gene_sets, names(pathways))]
            }
        }
        # gsea <- list()
        gsea <- apply(gene_ranks, MARGIN = 1, function(x) {
            fgsea::fgsea(
                pathways,
                x,
                minSize = 15,
                maxSize = 500,
                nperm = 10000
            )
        })
        gsea <- lapply(gsea, function(x) {
            as.matrix(x[order(x$pval), ])
        })
        return(gsea)
    }

#' Perform t-SNE dimensionality reduction
#'
#' Runs t-SNE on the normalized cell factors (or raw cell factors) to generate a 2D embedding for
#' visualization. Has option to run on subset of factors. Note that running multiple times will
#' reset tsne.coords values.
#'
#' In order to run fftRtsne (recommended for large datasets), you must first install FIt-SNE as
#' detailed \href{https://github.com/KlugerLab/FIt-SNE}{here}. Include the path to the cloned
#' FIt-SNE directory as the fitsne.path parameter, though this is only necessary for the first call
#' to runTSNE. For more detailed FIt-SNE installation instructions, see the liger repo README.
#'
#' @param object \code{liger} object. Should run quantile_norm before calling with defaults.
#' @param use.raw Whether to use un-aligned cell factor loadings (H matrices) (default FALSE).
#' @param dims.use Factors to use for computing tSNE embedding (default 1:ncol(H.norm)).
#' @param use.pca Whether to perform initial PCA step for Rtsne (default FALSE).
#' @param perplexity Parameter to pass to Rtsne (expected number of neighbors) (default 30).
#' @param theta Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE
#'   (default 0.5).
#' @param method Supports two methods for estimating tSNE values: Rtsne (Barnes-Hut implementation
#'   of t-SNE) and fftRtsne (FFT-accelerated Interpolation-based t-SNE) (using Kluger Lab
#'   implementation). (default Rtsne)
#' @param fitsne.path Path to the cloned FIt-SNE directory (ie. '/path/to/dir/FIt-SNE') (required
#'   for using fftRtsne -- only first time runTSNE is called) (default NULL).
#' @param rand.seed Random seed for reproducibility (default 42).
#'
#' @return \code{liger} object with tsne.coords slot set.
#'
#' @importFrom Rtsne Rtsne
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' # generate H.norm by quantile normalizig factor loadings
#' ligerex <- quantile_norm(ligerex)
#' # get tsne.coords for normalized data
#' ligerex <- runTSNE(ligerex)
#' # get tsne.coords for raw factor loadings
#' ligerex <- runTSNE(ligerex, use.raw = TRUE)
#' }

runTSNE <- function(object,
                    use.raw = FALSE,
                    dims.use = 1:ncol(object@H.norm),
                    use.pca = FALSE,
                    perplexity = 30,
                    theta = 0.5,
                    method = "Rtsne",
                    fitsne.path = NULL,
                    rand.seed = 42) {
    if (use.raw) {
        data.use <- do.call(rbind, object@H)
        if (identical(dims.use, 1:0)) {
            dims.use <- 1:ncol(data.use)
        }
    } else {
        data.use <- object@H.norm
    }
    if (method == "Rtsne") {
        set.seed(rand.seed)
        object@tsne.coords <- Rtsne(
            data.use[, dims.use],
            pca = use.pca,
            check_duplicates = FALSE,
            theta = theta,
            perplexity = perplexity
        )$Y
    } else if (method == "fftRtsne") {
        # if (!exists('fftRtsne')) {
        #   if (is.null(fitsne.path)) {
        #     stop('Please pass in path to FIt-SNE directory as fitsne.path.')
        #   }
        # source(paste0(fitsne.path, '/fast_tsne.R'), chdir = TRUE)
        # }
        # object@tsne.coords <- fftRtsne(data.use[, dims.use], rand_seed = rand.seed,
        #                                theta = theta, perplexity = perplexity)
        object@tsne.coords <- fftRtsne(
            X = data.use[, dims.use],
            rand_seed = rand.seed,
            fast_tsne_path = fitsne.path,
            theta = theta,
            perplexity = perplexity
        )
    } else {
        stop("Invalid method: Please choose Rtsne or fftRtsne")
    }
    rownames(object@tsne.coords) <- rownames(data.use)
    return(object)
}

#' Perform UMAP dimensionality reduction
#'
#' @description
#' Run UMAP on the normalized cell factors (or raw cell factors) to generate a 2D embedding for
#' visualization (or general dimensionality reduction). Has option to run on subset of factors.
#' Note that running multiple times will overwrite tsne.coords values. It is generally
#' recommended to use this method for dimensionality reduction with extremely large datasets.
#'
#' Note that this method requires that the package uwot is installed. It does not depend
#' on reticulate or python umap-learn.
#'
#' @param object \code{liger} object. Should run quantile_norm before calling with defaults.
#' @param use.raw Whether to use un-aligned cell factor loadings (H matrices) (default FALSE).
#' @param dims.use Factors to use for computing tSNE embedding (default 1:ncol(H.norm)).
#' @param k Number of dimensions to reduce to (default 2).
#' @param distance Mtric used to measure distance in the input space. A wide variety of metrics are
#'   already coded, and a user defined function can be passed as long as it has been JITd by numba.
#'   (default "euclidean", alternatives: "cosine", "manhattan", "hamming")
#' @param n_neighbors Number of neighboring points used in local approximations of manifold
#'   structure. Larger values will result in more global structure being preserved at the loss of
#'   detailed local structure. In general this parameter should often be in the range 5 to 50, with
#'   a choice of 10 to 15 being a sensible default. (default 10)
#' @param min_dist Controls how tightly the embedding is allowed compress points together. Larger
#'   values ensure embedded points are more evenly distributed, while smaller values allow the
#'   algorithm to optimise more accurately with regard to local structure. Sensible values are in
#'   the range 0.001 to 0.5, with 0.1 being a reasonable default. (default 0.1)
#' @param rand.seed Random seed for reproducibility (default 42).
#'
#' @return \code{liger} object with tsne.coords slot set.
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' # generate H.norm by quantile normalizig factor loadings
#' ligerex <- quantileAlignSNF(ligerex)
#' # get tsne.coords for normalized data
#' ligerex <- runUMAP(ligerex)
#' # get tsne.coords for raw factor loadings
#' ligerex <- runUMAP(ligerex, use.raw = TRUE)
#' }
runUMAP <- function(object,
                    use.raw = FALSE,
                    dims.use = 1:ncol(object@H.norm),
                    k = 2,
                    distance = "euclidean",
                    n_neighbors = 10,
                    min_dist = 0.1,
                    rand.seed = 42) {
    set.seed(rand.seed)
    if (use.raw) {
        raw.data <- do.call(rbind, object@H)
        # if H.norm not set yet
        if (identical(dims.use, 1:0)) {
            dims.use <- 1:ncol(raw.data)
        }
        object@tsne.coords <- uwot::umap(
            raw.data[, dims.use],
            n_components = as.integer(k),
            metric = distance,
            n_neighbors = as.integer(n_neighbors),
            min_dist = min_dist
        )
        rownames(object@tsne.coords) <- rownames(raw.data)
    } else {
        object@tsne.coords <- uwot::umap(
            object@H.norm[, dims.use],
            n_components = as.integer(k),
            metric = distance,
            n_neighbors = as.integer(n_neighbors),
            min_dist = min_dist
        )
        rownames(object@tsne.coords) <- rownames(object@H.norm)
    }
    return(object)
}

#' Find shared and dataset-specific markers
#'
#' Applies various filters to genes on the shared (W) and dataset-specific (V) components of the
#' factorization, before selecting those which load most significantly on each factor (in a shared
#' or dataset-specific way).
#'
#' @param object \code{liger} object. Should call optimizeALS before calling.
#' @param dataset1 Name of first dataset (default first dataset by order)
#' @param dataset2 Name of second dataset (default second dataset by order)
#' @param factor.share.thresh Use only factors with a dataset specificity less than or equalt to
#'   threshold (default 10).
#' @param dataset.specificity Pre-calculated dataset specificity if available. Will calculate if not
#'   available.
#' @param log.fc.thresh Lower log-fold change threshold for differential expression in markers
#'   (default 1).
#' @param pval.thresh Upper p-value threshold for Wilcoxon rank test for gene expression
#'   (default 0.05).
#' @param num.genes Max number of genes to report for each dataset (default 30).
#' @param print.genes Print ordered markers passing logfc, umi and frac thresholds (default FALSE).
#' @param verbose Print messages (TRUE by default)
#'
#' @return List of shared and specific factors. First three elements are dataframes of dataset1-
#'   specific, shared, and dataset2-specific markers. Last two elements are tables indicating the
#'   number of factors in which marker appears.
#'
#' @importFrom stats wilcox.test
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete input
#' markers <- getFactorMarkers(ligerex, num.genes = 10)
#' # look at shared markers
#' head(markers[[2]])
#' }

getFactorMarkers <- function(object, dataset1 = NULL, dataset2 = NULL, factor.share.thresh = 10,
                             dataset.specificity = NULL, log.fc.thresh = 1, pval.thresh = 0.05,
                             num.genes = 30, print.genes = FALSE, verbose = TRUE) {
    if (is.null(dataset1) | is.null(dataset2)) {
        dataset1 <- names(object@H)[1]
        dataset2 <- names(object@H)[2]
    }
    if (is.null(num.genes)) {
        num.genes <- length(object@var.genes)
    }
    if (is.null(dataset.specificity)) {
        dataset.specificity <- calcDatasetSpecificity(object, dataset1 = dataset1,
                                                      dataset2 = dataset2, do.plot = FALSE)
    }
    factors.use <- which(abs(dataset.specificity[[3]]) <= factor.share.thresh)

    if (length(factors.use) < 2 && verbose) {
        message(
            "Warning: only ", length(factors.use),
            " factors passed the dataset specificity threshold."
        )
    }

    Hs_scaled <- lapply(object@H, function(x) {
        scale(x, scale = TRUE, center = TRUE)
    })
    labels <- list()
    for (i in 1:length(Hs_scaled)) {
        if (class(object@raw.data[[1]])[1] == "H5File") {
            if (object@h5file.info[[1]][["sample.data.type"]] != "norm.data") {
                stop("norm.data should be sampled for obtaining factor markers.")
            }
            labels[[i]] <- factors.use[as.factor(apply(Hs_scaled[[i]][colnames(object@sample.data[[i]]), factors.use], 1, which.max))]
        } else {
            labels[[i]] <- factors.use[as.factor(apply(Hs_scaled[[i]][, factors.use], 1, which.max))]
        }
    }
    names(labels) <- names(object@H)

    V1_matrices <- list()
    V2_matrices <- list()
    W_matrices <- list()
    for (j in 1:length(factors.use)) {
        i <- factors.use[j]

        W <- t(object@W)
        V1 <- t(object@V[[dataset1]])
        V2 <- t(object@V[[dataset2]])
        rownames(W) <- rownames(V1) <- rownames(V2) <- object@var.genes
        # if not max factor for any cell in either dataset
        if (sum(labels[[dataset1]] == i) <= 1 | sum(labels[[dataset2]] == i) <= 1) {
            message("Warning: factor", i, "did not appear as max in any cell in either dataset")
            next
        }
        # filter genes by gene_count and cell_frac thresholds
        if (class(object@raw.data[[1]])[1] == "H5File") {
            if (object@h5file.info[[1]][["sample.data.type"]] != "norm.data") {
                stop("Sampled norm.data are required for this analysis.")
            }
            expr_mat = Reduce(cbind, object@sample.data[c(dataset1,dataset2)])[object@var.genes, c(labels[[dataset1]] == i, labels[[dataset2]] == i)]
            cell_label = rep(c(dataset1, dataset2), c(sum(labels[[dataset1]] == i), sum(labels[[dataset2]] == i)))
            wilcoxon_result = wilcoxauc(log(expr_mat + 1e-10), cell_label)

        } else {
            expr_mat = cbind(object@norm.data[[dataset1]][object@var.genes, labels[[dataset1]] == i],
                             object@norm.data[[dataset2]][object@var.genes, labels[[dataset2]] == i])
            cell_label = rep(c(dataset1, dataset2), c(sum(labels[[dataset1]] == i), sum(labels[[dataset2]] == i)))
            wilcoxon_result = wilcoxauc(log(expr_mat + 1e-10), cell_label)
        }
        log2fc = wilcoxon_result[wilcoxon_result$group == dataset1, ]$logFC
        names(log2fc) = wilcoxon_result[wilcoxon_result$group == dataset1, ]$feature
        filtered_genes_V1 = wilcoxon_result[wilcoxon_result$logFC > log.fc.thresh & wilcoxon_result$pval < pval.thresh, ]$feature
        filtered_genes_V2 = wilcoxon_result[-wilcoxon_result$logFC > log.fc.thresh & wilcoxon_result$pval < pval.thresh, ]$feature

        W <- pmin(W + V1, W + V2)
        V1 <- V1[filtered_genes_V1, , drop = FALSE]
        V2 <- V2[filtered_genes_V2, , drop = FALSE]

        if (length(filtered_genes_V1) == 0) {
            top_genes_V1 <- character(0)
        } else {
            top_genes_V1 <- row.names(V1)[order(V1[, i], decreasing = TRUE)[1:num.genes]]
            top_genes_V1 <- top_genes_V1[!is.na(top_genes_V1)]
            top_genes_V1 <- top_genes_V1[which(V1[top_genes_V1, i] > 0)]
        }
        if (length(filtered_genes_V2) == 0) {
            top_genes_V2 <- character(0)
        } else {
            top_genes_V2 <- row.names(V2)[order(V2[, i], decreasing = TRUE)[1:num.genes]]
            top_genes_V2 <- top_genes_V2[!is.na(top_genes_V2)]
            top_genes_V2 <- top_genes_V2[which(V2[top_genes_V2, i] > 0)]
        }
        top_genes_W <- row.names(W)[order(W[, i], decreasing = TRUE)[1:num.genes]]
        top_genes_W <- top_genes_W[!is.na(top_genes_W)]
        top_genes_W <- top_genes_W[which(W[top_genes_W, i] > 0)]

        if (print.genes && verbose) {
            message("Factor ", i)
            message('Dataset 1')
            message(top_genes_V1)
            message('Shared')
            message(top_genes_W)
            message('Dataset 2')
            message(top_genes_V2)
        }

        pvals <- list() # order is V1, V2, W
        top_genes <- list(top_genes_V1, top_genes_V2, top_genes_W)
        for (k in 1:length(top_genes)) {
            pvals[[k]] <- wilcoxon_result[wilcoxon_result$feature %in% top_genes[[k]] & wilcoxon_result$group == dataset1, ]$pval
        }
        # bind values in matrices
        V1_matrices[[j]] <- Reduce(cbind, list(
            rep(i, length(top_genes_V1)), top_genes_V1,
            log2fc[top_genes_V1], pvals[[1]]
        ))
        V2_matrices[[j]] <- Reduce(cbind, list(
            rep(i, length(top_genes_V2)), top_genes_V2,
            log2fc[top_genes_V2], pvals[[2]]
        ))
        W_matrices[[j]] <- Reduce(cbind, list(
            rep(i, length(top_genes_W)), top_genes_W,
            log2fc[top_genes_W], pvals[[3]]
        ))
    }
    V1_genes <- data.frame(Reduce(rbind, V1_matrices), stringsAsFactors = FALSE)
    V2_genes <- data.frame(Reduce(rbind, V2_matrices), stringsAsFactors = FALSE)
    W_genes <- data.frame(Reduce(rbind, W_matrices), stringsAsFactors = FALSE)
    df_cols <- c("factor_num", "gene", "log2fc", "p_value")
    output_list <- list(V1_genes, W_genes, V2_genes)
    output_list <- lapply(seq_along(output_list), function(x) {
        df <- output_list[[x]]
        colnames(df) <- df_cols
        df <- transform(df,
                        factor_num = as.numeric(df$'factor_num'), gene = as.character(df$'gene'),
                        log2fc = as.numeric(df$'log2fc'), p_value = as.numeric(df$'p_value')
        )
        # Cutoff only applies to dataset-specific dfs
        if (x != 2) {
            df[which(df$p_value < pval.thresh), ]
        } else {
            df
        }
    })
    names(output_list) <- c(dataset1, "shared", dataset2)
    output_list[["num_factors_V1"]] <- table(output_list[[dataset1]]$'gene')
    output_list[["num_factors_V2"]] <- table(output_list[[dataset2]]$'gene')
    return(output_list)
}
