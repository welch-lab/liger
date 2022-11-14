#' @importFrom Matrix colSums t
NULL

#' Calculate a dataset-specificity score for each factor
#'
#' This score represents the relative magnitude of the dataset-specific components of each factor's
#' gene loadings compared to the shared components for two datasets. First, for each dataset we
#' calculate the norm of the sum of each factor's shared loadings (W) and dataset-specific loadings
#' (V). We then determine the ratio of these two values and subtract from 1... TODO: finish
#' description.
#'
#' @param object \code{liger} object. Should run optimizeALS before calling.
#' @param dataset1 Name of first dataset (by default takes first two datasets for dataset1 and 2)
#' @param dataset2 Name of second dataset
#' @param do.plot Display barplot of dataset specificity scores (by factor) (default TRUE).
#'
#' @return List containing three elements. First two elements are the norm of each metagene factor
#' for each dataset. Last element is the vector of dataset specificity scores.
#'
#' @importFrom graphics barplot
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' # generate H.norm by quantile normalizig factor loadings
#' ligerex <- quantile_norm(ligerex)
#' dataset_spec <- calcDatasetSpecificity(ligerex, do.plot = F)
#' }

calcDatasetSpecificity <- function(object,
                                   dataset1 = NULL,
                                   dataset2 = NULL,
                                   do.plot = TRUE) {
    if (is.null(dataset1) | is.null(dataset2)) {
        dataset1 <- names(object@H)[1]
        dataset2 <- names(object@H)[2]
    }
    k <- ncol(object@H[[1]])
    pct1 <- rep(0, k)
    pct2 <- rep(0, k)
    for (i in 1:k) {
        pct1[i] <-
            norm(as.matrix(object@V[[dataset1]][i,] + object@W[i,]), "F")
        # norm(object@H[[1]][,i] %*% t(object@W[i,] + object@V[[1]][i,]),"F")
        pct2[i] <-
            norm(as.matrix(object@V[[dataset2]][i,] + object@W[i,]), "F")
        # norm(object@H[[2]][,i] %*% t(object@W[i,] + object@V[[2]][i,]),"F")
    }
    # pct1 = pct1/sum(pct1)
    # pct2 = pct2/sum(pct2)
    if (do.plot) {
        graphics::barplot(
            100 * (1 - (pct1 / pct2)),
            xlab = "Factor",
            ylab = "Percent Specificity",
            main = "Dataset Specificity of Factors",
            names.arg = 1:k,
            cex.names = 0.75,
            mgp = c(2, 0.5, 0)
        ) # or possibly abs(pct1-pct2)
    }
    return(list(pct1, pct2, 100 * (1 - (pct1 / pct2))))
}

#' Calculate agreement metric
#'
#' @description
#' This metric quantifies how much the factorization and alignment distorts the geometry of the
#' original datasets. The greater the agreement, the less distortion of geometry there is. This is
#' calculated by performing dimensionality reduction on the original and quantile aligned (or just
#' factorized) datasets, and measuring similarity between the k nearest neighbors for each cell in
#' original and aligned datasets. The Jaccard index is used to quantify similarity, and is the final
#' metric averages across all cells.
#'
#' Note that for most datasets, the greater the chosen k, the greater the agreement in general.
#' There are several options for dimensionality reduction, with the default being 'NMF' as it is
#' expected to be most similar to iNMF. Although agreement can theoretically approach 1, in practice
#' it is usually no higher than 0.2-0.3 (particularly for non-deterministic approaches like NMF).
#'
#' @param object \code{liger} object. Should call quantile_norm before calling.
#' @param dr.method Dimensionality reduction method to use for assessing pre-alignment geometry
#'   (either "PCA", "NMF", or "ICA"). (default "NMF")
#' @param ndims Number of dimensions to use in dimensionality reduction (recommended to use the
#'   same as number of factors) (default 40).
#' @param k Number of nearest neighbors to use in calculating Jaccard index (default 15).
#' @param use.aligned Whether to use quantile aligned or unaligned cell factor loadings (default
#'   TRUE).
#' @param rand.seed Random seed for reproducibility (default 42).
#' @param by.dataset Return agreement calculated for each dataset (default FALSE).
#'
#' @return Agreement metric (or vector of agreement per dataset).
#'
#' @importFrom FNN get.knn
#' @importFrom ica icafast
#' @importFrom irlba prcomp_irlba
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object based on in-memory datasets), factorization complete
#' # generate H.norm by quantile normalizig factor loadings
#' ligerex <- quantile_norm(ligerex)
#' agreement <- calcAgreement(ligerex, dr.method = "NMF")
#' # ligerex (liger object based on datasets in HDF5 format), factorization complete
#' ligerex <- quantile_norm(ligerex)
#' ligerex <- readSubset(ligerex, slot.use = "scale.data", max.cells = 5000)
#' agreement <- calcAgreement(ligerex, dr.method = "NMF")
#' }

calcAgreement <- function(object,
                          dr.method = "NMF",
                          ndims = 40,
                          k = 15,
                          use.aligned = TRUE,
                          rand.seed = 42,
                          by.dataset = FALSE) {
    # if (!requireNamespace("NNLM", quietly = TRUE) & dr.method == "NMF") {
    #   stop("Package \"NNLM\" needed for this function to perform NMF. Please install it.",
    #        call. = FALSE
    #   )
    # }
    if (class(object@raw.data[[1]])[1] == "H5File") {
        if (object@h5file.info[[1]][["sample.data.type"]] != "scale.data") {
            stop(
                "HDF5-based Liger object requires sampled scale.data for calculating agreement."
            )
        }
    }

    message("Reducing dimensionality using ", dr.method)
    set.seed(rand.seed)
    dr <- list()
    if (dr.method == "NMF") {
        if (class(object@raw.data[[1]])[1] == "H5File") {
            dr <- lapply(object@sample.data, function(x) {
                nmf_hals(t(x), k = ndims)[[1]]
            })
        } else {
            dr <- lapply(object@scale.data, function(x) {
                nmf_hals(x, k = ndims)[[1]]
            })
        }
    }
    else if (dr.method == "ICA") {
        if (class(object@raw.data[[1]])[1] == "H5File") {
            dr <- list()
            for (i in 1:length(object@H)) {
                dr[[i]] = icafast(t(object@sample.data[[i]]), nc = ndims)$S
            }

        } else {
            dr <- lapply(object@scale.data, function(x) {
                icafast(x, nc = ndims)$S
            })
        }
    } else {
        if (class(object@raw.data[[1]])[1] == "H5File") {
            dr <- list()
            for (i in 1:length(object@H)) {
                dr[[i]] = suppressWarnings(
                    prcomp_irlba(
                        object@sample.data[[i]],
                        n = ndims,
                        scale. = (colSums(t(
                            object@sample.data[[i]]
                        )) > 0),
                        center = FALSE
                    )$rotation
                )
                rownames(dr[[i]]) = colnames(object@sample.data[[i]])
            }

        } else {
            dr <- lapply(object@scale.data, function(x) {
                suppressWarnings(prcomp_irlba(
                    t(x),
                    n = ndims,
                    scale. = (colSums(x) > 0),
                    center = FALSE
                )$rotation)
            })
            for (i in 1:length(dr)) {
                rownames(dr[[i]]) <- rownames(object@scale.data[[i]])
            }
        }
    }
    ns <- sapply(dr, nrow)
    n <- sum(ns)
    jaccard_inds <- c()
    distorts <- c()

    for (i in 1:length(dr)) {
        jaccard_inds_i <- c()
        if (use.aligned) {
            original <- object@H.norm[rownames(dr[[i]]),]
        } else {
            original <- object@H[[i]]
        }
        fnn.1 <- get.knn(dr[[i]], k = k)
        fnn.2 <- get.knn(original, k = k)
        jaccard_inds_i <-
            c(jaccard_inds_i, sapply(1:ns[i], function(i) {
                intersect <- intersect(fnn.1$nn.index[i,], fnn.2$nn.index[i,])
                union <- union(fnn.1$nn.index[i,], fnn.2$nn.index[i,])
                length(intersect) / length(union)
            }))
        jaccard_inds_i <- jaccard_inds_i[is.finite(jaccard_inds_i)]
        jaccard_inds <- c(jaccard_inds, jaccard_inds_i)

        distorts <- c(distorts, mean(jaccard_inds_i))
    }
    if (by.dataset) {
        return(distorts)
    }
    return(mean(jaccard_inds))
}

#' Calculate alignment metric
#'
#' This metric quantifies how well-aligned two or more datasets are. Alignment is defined as in the
#' documentation for Seurat. We randomly downsample all datasets to have as many cells as the
#' smallest one. We construct a nearest-neighbor graph and calculate for each cell how many of its
#' neighbors are from the same dataset. We average across all cells and compare to the expected
#' value for perfectly mixed datasets, and scale the value from 0 to 1. Note that in practice,
#' alignment can be greater than 1 occasionally.
#'
#' @param object \code{liger} object. Should call quantile_norm before calling.
#' @param k Number of nearest neighbors to use in calculating alignment. By default, this will be
#'   floor(0.01 * total number of cells), with a lower bound of 10 in all cases except where the
#'   total number of sampled cells is less than 10.
#' @param rand.seed Random seed for reproducibility (default 1).
#' @param cells.use Vector of cells across all datasets to use in calculating alignment
#' @param cells.comp Vector of cells across all datasets to compare to cells.use when calculating
#'   alignment (instead of dataset designations). These can be from the same dataset as cells.use.
#'   (default NULL)
#' @param clusters.use Names of clusters to use in calculating alignment (default NULL).
#' @param by.cell Return alignment calculated individually for each cell (default FALSE).
#' @param by.dataset Return alignment calculated for each dataset (default FALSE).
#'
#' @return Alignment metric.
#'
#' @importFrom FNN get.knn
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object ), factorization complete
#' ligerex <- quantile_norm(ligerex)
#' alignment <- calcAlignment(ligerex)
#' }

calcAlignment <- function(object,
                          k = NULL,
                          rand.seed = 1,
                          cells.use = NULL,
                          cells.comp = NULL,
                          clusters.use = NULL,
                          by.cell = FALSE,
                          by.dataset = FALSE) {
    if (is.null(cells.use)) {
        cells.use <- rownames(object@H.norm)
    }
    if (!is.null(clusters.use)) {
        cells.use <-
            names(object@clusters)[which(object@clusters %in% clusters.use)]
    }
    if (!is.null(cells.comp)) {
        nmf_factors <- object@H.norm[c(cells.use, cells.comp),]
        num_cells <- length(c(cells.use, cells.comp))
        func_H <- list(cells1 = nmf_factors[cells.use,],
                       cells2 = nmf_factors[cells.comp,])
        message('Using designated sets cells.use and cells.comp as subsets to compare')
    } else {
        nmf_factors <- object@H.norm[cells.use,]
        num_cells <- length(cells.use)
        func_H <- lapply(seq_along(object@H), function(x) {
            cells.overlap <- intersect(cells.use, rownames(object@H[[x]]))
            if (length(cells.overlap) > 0) {
                object@H[[x]][cells.overlap,]
            } else {
                warning(paste0(
                    "Selected subset eliminates dataset ",
                    names(object@H)[x]
                ),
                immediate. = TRUE)
                return(NULL)
            }
        })
        func_H <- func_H[!sapply(func_H, is.null)]
    }
    num_factors <- ncol(object@H.norm)
    N <- length(func_H)
    if (N == 1) {
        warning("Alignment null for single dataset", immediate. = TRUE)
    }
    set.seed(rand.seed)
    min_cells <- min(sapply(func_H, function(x) {
        nrow(x)
    }))
    sampled_cells <- unlist(lapply(1:N, function(x) {
        sample(rownames(func_H[[x]]), min_cells)
    }))
    max_k <- length(sampled_cells) - 1
    if (is.null(k)) {
        k <- min(max(floor(0.01 * num_cells), 10), max_k)
    } else if (k > max_k) {
        stop(paste0("Please select k <=", max_k))
    }
    knn_graph <-
        get.knn(nmf_factors[sampled_cells, 1:num_factors], k)
    # Generate new "datasets" for desired cell groups
    if (!is.null(cells.comp)) {
        dataset <- unlist(sapply(1:N, function(x) {
            rep(paste0('group', x), nrow(func_H[[x]]))
        }))
    } else {
        dataset <- unlist(sapply(1:N, function(x) {
            rep(names(object@H)[x], nrow(func_H[[x]]))
        }))
    }
    names(dataset) <- rownames(nmf_factors)
    dataset <- dataset[sampled_cells]

    num_sampled <- N * min_cells
    num_same_dataset <- rep(k, num_sampled)

    alignment_per_cell <- c()
    for (i in 1:num_sampled) {
        inds <- knn_graph$nn.index[i,]
        num_same_dataset[i] <- sum(dataset[inds] == dataset[i])
        alignment_per_cell[i] <-
            1 - (num_same_dataset[i] - (k / N)) / (k - k / N)
    }
    if (by.dataset) {
        alignments <- c()
        for (i in 1:N) {
            start <- 1 + (i - 1) * min_cells
            end <- i * min_cells
            alignment <- mean(alignment_per_cell[start:end])
            alignments <- c(alignments, alignment)
        }
        return(alignments)
    } else if (by.cell) {
        names(alignment_per_cell) <- sampled_cells
        return(alignment_per_cell)
    }
    return(mean(alignment_per_cell))
}

#' Calculate alignment for each cluster
#'
#' Returns alignment for each cluster in analysiss (see documentation for calcAlignment).
#'
#' @param object \code{liger} object. Should call quantileAlignSNF before calling.
#' @param rand.seed Random seed for reproducibility (default 1).
#' @param k Number of nearest neighbors in calculating alignment (see calcAlignment for default).
#'   Can pass in single value or vector with same length as number of clusters.
#' @param by.dataset Return alignment calculated for each dataset in cluster (default FALSE).
#'
#' @return Vector of alignment statistics (with names of clusters).
#'
#' @importFrom FNN get.knn
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' ligerex <- quantile_norm(ligerex)
#' # get alignment for each cluster
#' alignment_per_cluster <- calcAlignmentPerCluster(ligerex)
#' }

calcAlignmentPerCluster <- function(object,
                                    rand.seed = 1,
                                    k = NULL,
                                    by.dataset = FALSE) {
    clusters <- levels(object@clusters)
    if (typeof(k) == "double") {
        if (length(k) == 1) {
            k <- rep(k, length(clusters))
        } else if (length(k) != length(clusters)) {
            stop("Length of k does not match length of clusters")
        }
    }
    align_metrics <- sapply(seq_along(clusters), function(x) {
        calcAlignment(
            object,
            k = k[x],
            rand.seed = rand.seed,
            clusters.use = clusters[x],
            by.dataset = by.dataset
        )
    })
    if (by.dataset) {
        colnames(align_metrics) <- levels(object@clusters)
        rownames(align_metrics) <- names(object@H)
    } else {
        names(align_metrics) <- levels(object@clusters)
    }
    return(align_metrics)
}

#' Calculate adjusted Rand index
#'
#' Computes adjusted Rand index for \code{liger} clustering and external clustering.
#' The Rand index ranges from 0 to 1, with 0 indicating no agreement between clusterings and 1
#' indicating perfect agreement.
#'
#' @param object \code{liger} object. Should run quantileAlignSNF before calling.
#' @param clusters.compare Clustering with which to compare (named vector).
#' @param verbose Print messages (TRUE by default)
#'
#' @return Adjusted Rand index value.
#'
#' @importFrom mclust adjustedRandIndex
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' ligerex <- quantile_norm(ligerex)
#' # toy clusters
#' cluster1 <- sample(c('type1', 'type2', 'type3'), ncol(ligerex@raw.data[[1]]), replace = TRUE)
#' names(cluster1) <- colnames(ligerex@raw.data[[1]])
#' cluster2 <- sample(c('type4', 'type5', 'type6'), ncol(ligerex@raw.data[[2]]), replace = TRUE)
#' names(cluster2) <- colnames(ligerex@raw.data[[2]])
#' # get ARI for first clustering
#' ari1 <- calcARI(ligerex, cluster1)
#' # get ARI for second clustering
#' ari2 <- calcARI(ligerex, cluster2)
#' }

calcARI <- function(object, clusters.compare, verbose = TRUE) {
    if (length(clusters.compare) < length(object@clusters) && verbose) {
        message("Calculating ARI for subset of all cells")
    }
    return(adjustedRandIndex(object@clusters[names(clusters.compare)],
                             clusters.compare))
}

#' Calculate purity
#'
#' Calculates purity for \code{liger} clustering and external clustering (true clusters/classes).
#' Purity can sometimes be a more useful metric when the clustering to be tested contains more
#' subgroups or clusters than the true clusters (or classes). Purity also ranges from 0 to 1,
#' with a score of 1 representing a pure, or accurate, clustering.
#'
#' @param object \code{liger} object. Should run quantileAlignSNF before calling.
#' @param classes.compare Clustering with which to compare (named vector).
#' @param verbose Print messages (TRUE by default)
#'
#' @return Purity value.
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' ligerex <- quantile_norm(ligerex)
#' # toy clusters
#' cluster1 <- sample(c('type1', 'type2', 'type3'), ncol(ligerex@raw.data[[1]]), replace = TRUE)
#' names(cluster1) <- colnames(ligerex@raw.data[[1]])
#' cluster2 <- sample(c('type4', 'type5', 'type6'), ncol(ligerex@raw.data[[2]]), replace = TRUE)
#' names(cluster2) <- colnames(ligerex@raw.data[[2]])
#' # get ARI for first clustering
#' ari1 <- calcPurity(ligerex, cluster1)
#' # get ARI for second clustering
#' ari2 <- calcPurity(ligerex, cluster2)
#' }

calcPurity <- function(object, classes.compare, verbose = TRUE) {
    if (length(classes.compare) < length(object@clusters) && verbose) {
        print("Calculating purity for subset of full cells")
    }
    clusters <- object@clusters[names(classes.compare)]
    purity <-
        sum(apply(table(classes.compare, clusters), 2, max)) / length(clusters)

    return(purity)
}

#' Calculate proportion mitochondrial contribution
#'
#' Calculates proportion of mitochondrial contribution based on raw or normalized data.
#'
#' @param object \code{liger} object.
#' @param use.norm Whether to use cell normalized data in calculating contribution (default FALSE).
#'
#' @return Named vector containing proportion of mitochondrial contribution for each cell.
#'
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' ligerex@cell.data[["percent_mito"]] <- getProportionMito(ligerex)
#' }

getProportionMito <- function(object, use.norm = FALSE) {
    all.genes <- Reduce(union, lapply(object@raw.data, rownames))
    mito.genes <-
        grep(pattern = "^mt-",
             x = all.genes,
             value = TRUE)
    data.use <- object@raw.data
    if (use.norm) {
        data.use <- object@norm.data
    }
    percent_mito <- unlist(lapply(unname(data.use), function(x) {
        colSums(x[mito.genes,]) / colSums(x)
    }), use.names = TRUE)

    return(percent_mito)
}
