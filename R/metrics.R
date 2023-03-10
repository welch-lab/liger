#' [Deprecated] Calculate agreement metric
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
#' @name calcAgreement-deprecated
#' @seealso \code{\link{rliger-deprecated}}
NULL

#' @rdname rliger-deprecated
#' @export
#' @section \code{calcAgreement}:
#' \code{calcAgreement} is deprecated and is now defunct.
calcAgreement <- function(object,
                          dr.method = "NMF",
                          ndims = 40,
                          k = 15,
                          use.aligned = TRUE,
                          rand.seed = 42,
                          by.dataset = FALSE) {
    lifecycle::deprecate_stop("1.99.0", "calcAgreement()")
    # if (!requireNamespace("NNLM", quietly = TRUE) & dr.method == "NMF") {
    #   stop("Package \"NNLM\" needed for this function to perform NMF. Please install it.",
    #        call. = FALSE
    #   )
    # }
    if (class(object@raw.data[[1]])[1] == "H5File") {
        if (object@h5fileInfo[[1]][["sample.data.type"]] != "scale.data") {
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
        if (!requireNamespace("ica", quietly = TRUE))
            stop("Package \"ica\" needed for this function to work. ",
                 "Please install it by command:\n",
                 "BiocManager::install('ica')",
                 call. = FALSE)
        if (class(object@raw.data[[1]])[1] == "H5File") {
            dr <- list()
            for (i in 1:length(object@H)) {
                dr[[i]] = ica::icafast(t(object@sample.data[[i]]), nc = ndims)$S
            }

        } else {
            dr <- lapply(object@scale.data, function(x) {
                ica::icafast(x, nc = ndims)$S
            })
        }
    } else {
        if (!requireNamespace("irlba", quietly = TRUE))
            stop("Package \"irlba\" needed for this function to work. ",
                 "Please install it by command:\n",
                 "BiocManager::install('irlba')",
                 call. = FALSE)
        if (class(object@raw.data[[1]])[1] == "H5File") {
            dr <- list()
            for (i in 1:length(object@H)) {
                dr[[i]] = suppressWarnings(
                    irlba::prcomp_irlba(
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
                suppressWarnings(irlba::prcomp_irlba(
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
        fnn.1 <- FNN::get.knn(dr[[i]], k = k)
        fnn.2 <- FNN::get.knn(original, k = k)
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
#' @export
#' @examples
#' \dontrun{
#' # ligerex (liger object ), factorization complete
#' ligerex <- quantile_norm(ligerex)
#' alignment <- calcAlignment(ligerex)
#' }
calcAlignment <- function(
        object,
        k = NULL,
        rand.seed = 1,
        cells.use = NULL,
        cells.comp = NULL,
        clusters.use = NULL,
        by.cell = FALSE,
        by.dataset = FALSE
) {
    if (is.null(cells.use)) {
        cells.use <- colnames(object)
    }
    if (!is.null(clusters.use)) {
        cells.use <- names(object@clusters)[which(object@clusters %in% clusters.use)]
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
        FNN::get.knn(nmf_factors[sampled_cells, 1:num_factors], k)
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

#' [Deprecated] Calculate alignment for each cluster
#'
#' Returns alignment for each cluster in analysiss (see documentation for calcAlignment).
#'
#' @param object \code{liger} object. Should run
#' \code{\link{quantileNorm}} before calling.
#' @param rand.seed Random seed for reproducibility (default 1).
#' @param k Number of nearest neighbors in calculating alignment (see calcAlignment for default).
#'   Can pass in single value or vector with same length as number of clusters.
#' @param by.dataset Return alignment calculated for each dataset in cluster (default FALSE).
#'
#' @return Vector of alignment statistics (with names of clusters).
#'
#' @examples
#' \dontrun{
#' # ligerex (liger object), factorization complete
#' ligerex <- quantile_norm(ligerex)
#' # get alignment for each cluster
#' alignment_per_cluster <- calcAlignmentPerCluster(ligerex)
#' }
#' @name calcAlignmentPerCluster-deprecated
#' @seealso \code{\link{rliger-deprecated}}
NULL

#' @rdname rliger-deprecated
#' @export
#' @section \code{calcAlignmentPerCluster}:
#' \code{calcAlignmentPerCluster} is deprecated and is now defunct.
calcAlignmentPerCluster <- function(object,
                                    rand.seed = 1,
                                    k = NULL,
                                    by.dataset = FALSE) {
    lifecycle::deprecate_stop("1.99.0", "calcAlignmentPerCluster()")
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

#' [Deprecated] Calculate adjusted Rand index
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
#' @name calcARI-deprecated
#' @seealso \code{\link{rliger-deprecated}}
NULL

#' @rdname rliger-deprecated
#' @export
#' @section \code{calcARI}:
#' \code{calcARI} is deprecated and is now defunct.
calcARI <- function(object, clusters.compare, verbose = TRUE) {
    lifecycle::deprecate_stop("1.99.0", "calcARI()")
    if (!requireNamespace("mclust", quietly = TRUE))
        stop("Package \"mclust\" needed for this function to work. ",
             "Please install it by command:\n",
             "BiocManager::install('mclust')",
             call. = FALSE)
    if (length(clusters.compare) < length(object@clusters) && verbose) {
        message("Calculating ARI for subset of all cells")
    }
    return(mclust::adjustedRandIndex(object@clusters[names(clusters.compare)],
                             clusters.compare))
}

#' [Deprecated] Calculate purity
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
#' @name calcPurity-deprecated
#' @seealso \code{\link{rliger-deprecated}}
NULL

#' @rdname rliger-deprecated
#' @export
#' @section \code{calcPurity}:
#' \code{calcPurity} is deprecated and is now defunct.
calcPurity <- function(object, classes.compare, verbose = TRUE) {
    lifecycle::deprecate_stop("1.99.0", "calcPurity()")
    if (length(classes.compare) < length(object@clusters) && verbose) {
        print("Calculating purity for subset of full cells")
    }
    clusters <- object@clusters[names(classes.compare)]
    purity <-
        sum(apply(table(classes.compare, clusters), 2, max)) / length(clusters)

    return(purity)
}

#' Calculate loadings for each factor
#' @description Calculates the contribution of each factor of W,V, and U to the
#' reconstruction.
#' @param object A \linkS4class{liger} object. Should call
#' \code{\link{quantileNorm}} before calling.
#' @return A dataframe, such that each column represents the contribution of a
#' specific matrix (W, V_1, V_2, etc. )
#' @export
calcNormLoadings = function(object) {
    .checkValidFactorResult(object, names(object))
    # c x k
    H_norm = getMatrix(object, "H.norm")
    # k x g
    W_norm = t(getMatrix(object, "W"))
    # k x g
    V_norm = lapply(getMatrix(object, "V"), t)
    # k x g
    U_norm = lapply(getMatrix(object, "U"), function(x) {
        if (!is.null(x)) t(x)
        else NULL
    })
    considerU <- !any(sapply(U_norm, is.null))
    ##### Calculation of Contribution #########################
    w_loadings = list()
    u_loadings = rep(list(list()), length(object))
    v_loadings = rep(list(list()), length(object))
    k <- object@uns$factorization$k
    for (i in seq(k)) {
        hi = as.matrix(H_norm[, i])
        ####### Calculate W
        wi = t(as.matrix(W_norm[i, ]))
        hw = hi %*% wi
        forb_hw = norm(hw, type = "F") / dim(W_norm)[[2]]
        w_loadings = append(w_loadings, forb_hw)

        ###### Calculate V
        for (j in seq_along(object)) {
            temp_v = t(as.matrix(V_norm[[j]][i, ]))
            hv_temp = hi %*% temp_v
            forb_hv = norm(hv_temp, type = "F") / dim(V_norm[[j]])[[2]]
            v_loadings[[j]] = append(v_loadings[[j]], forb_hv)
        }
        if (considerU) {
            ###### Calculate U
            for (j in seq_along(object)) {
                if (length(U_norm[[j]]) != 0) {
                    temp_u = t(as.matrix(U_norm[[j]][i, ]))
                    hu_temp = hi %*% temp_u
                    forb_hu = norm(hu_temp, type = "F") / dim(U_norm[[j]])[[2]]
                    u_loadings[[j]] = append(u_loadings[[j]], forb_hu)
                }
            }
        }
    }

    ################# Format the return object
    w_loadings = unlist(w_loadings)
    factors = seq(k)
    results = data.frame(factors, w_loadings)

    # For all V
    for (j in seq_along(object)) {
        results = cbind(results, unlist(v_loadings[[j]]))
        colnames(results)[[2 + j]] = paste0("V_", j, "_loadings")
    }
    if (considerU) {
        # For all U
        for (j in seq_along(object)) {
            name_di = dim(results)[[2]]
            if (length(U_norm[[j]]) != 0) {
                results = cbind(results, unlist(u_loadings[[j]]))
                colnames(results)[[name_di + 1]] = paste0("U_", j, "_loadings")
            }
        }
    }

    return(results)
}

# Hierarchical alternating least squares for regular NMF
nmf_hals <- function(A, k, max_iters = 500, thresh = 1e-4, reps = 20, W0 = NULL, H0 = NULL) {
    m <- nrow(A)
    n <- ncol(A)
    if (is.null(W0)) {
        W0 <- matrix(abs(stats::runif(m * k, min = 0, max = 2)), m, k)
    }
    if (is.null(H0)) {
        H0 <- matrix(abs(stats::runif(n * k, min = 0, max = 2)), n, k)
    }
    W <- W0
    rownames(W) <- rownames(A)
    H <- H0
    rownames(H) <- colnames(A)

    # alpha = frobenius_prod(A %*% H, W)/frobenius_prod(t(W)%*%W,t(H)%*%H)
    # W = alpha*W

    for (i in 1:k)
    {
        W[, i] <- W[, i] / sum(W[, i]^2)
    }

    delta <- 1
    iters <- 0
    # pb <- txtProgressBar(min=0,max=max_iters,style=3)
    iter_times <- rep(0, length(max_iters))
    objs <- rep(0, length(max_iters))
    obj <- norm(A - W %*% t(H), "F")^2
    # print(obj)
    while (delta > thresh & iters < max_iters) {
        start_time <- Sys.time()
        obj0 <- obj
        HtH <- t(H) %*% H
        AH <- A %*% H
        for (i in 1:k)
        {
            W[, i] <- nonneg(W[, i] + (AH[, i] - (W %*% HtH[, i])) / (HtH[i, i]))
            W[, i] <- W[, i] / sqrt(sum(W[, i]^2))
        }

        AtW <- t(A) %*% W
        WtW <- t(W) %*% W
        for (i in 1:k)
        {
            H[, i] <- nonneg(H[, i] + (AtW[, i] - (H %*% WtW[, i])) / (WtW[i, i]))
        }

        obj <- norm(A - W %*% t(H), "F")^2
        # print(obj)
        delta <- abs(obj0 - obj) / (mean(obj0, obj))
        iters <- iters + 1
        end_time <- Sys.time()
        iter_times[iters] <- end_time - start_time
        objs[iters] <- obj
        # setTxtProgressBar(pb,iters)
    }
    .log("Converged in", end_time - start_time, "seconds,", iters, "iterations. Objective:", obj, "\n")

    return(list(W, H, cumsum(iter_times), objs))
}
