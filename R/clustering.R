#' Louvain algorithm for community detection
#' @description After quantile normalization, users can additionally run the
#' Louvain algorithm for community detection, which is widely used in
#' single-cell analysis and excels at merging small clusters into broad cell
#' classes.
#'
#' While using quantile normalized factor loadings (result from
#' \code{\link{quantileNorm}}) is recommended, this function looks for
#' unnormalized factor loadings (result from \code{\link{optimizeALS}} or
#' \code{\link{online_iNMF}}) when the former is not available.
#' @param object \linkS4class{liger} object. Should run
#' \code{\link{quantileNorm}} before calling.
#' @param nNeighbors Integer, the maximum number of nearest neighbors to
#' compute. Default \code{20}.
#' @param resolution Numeric, value of the resolution parameter, a larger value
#' results in a larger number of communities with smaller sizes. Default
#' \code{1.0}.
#' @param prune Numeric. Sets the cutoff for acceptable Jaccard index when
#' computing the neighborhood overlap for the SNN construction. Any edges with
#' values less than or equal to this will be set to 0 and removed from the SNN
#' graph. Essentially sets the stringency of pruning. \code{0} for no pruning,
#' while \code{1} prunes everything. Default \code{1/15}.
#' @param eps Numeric, the error bound of the nearest neighbor search. Default
#' \code{0.1}.
#' @param nRandomStarts Integer number of random starts. Default \code{10}.
#' @param nIterations Integer, maximal number of iterations per random start.
#' Default \code{100}.
#' @param useDims Indices of factors to use for Louvain clustering. Default
#' \code{NULL} uses all available factors.
#' @param groupSingletons Whether to group single cells that make up their own
#' cluster in with the cluster they are most connected to. Default \code{TRUE},
#' if \code{FALSE}, assign all singletons to a \code{"singleton"} group.
#' @param clusterName Name of the variable that will store the clustering result
#' in \code{cellMeta} slot of \code{object}. Default \code{"louvain_cluster"}.
#' @param seed Seed of the random number generator. Default \code{1}.
#' @param verbose Logical. Whether to show information of the progress.
#' Default \code{TRUE}.
#' @param k,dims.use,random.seed \bold{Deprecated}. See Usage section for
#' replacement.
#' @return \code{object} with refined cluster assignment updated in
#' \code{louvain_cluster} variable in \code{cellMeta} slot. Can be fetched
#' with \code{object$louvain_cluster}
#' @export
#' @examples
#' data("pbmcPlot", package = "rliger")
#' pbmc <- louvainCluster(pbmcPlot)
louvainCluster <- function(
        object,
        resolution = 1.0,
        nNeighbors = 20,
        prune = 1 / 15,
        eps = 0.1,
        nRandomStarts = 10,
        nIterations = 100,
        useDims = NULL,
        groupSingletons = TRUE,
        clusterName = "louvain_cluster",
        seed = 1,
        verbose = TRUE,
        # Deprecated coding style
        k = nNeighbors,
        dims.use = NULL,
        random.seed = 1
) {
    .deprecateArgs(list(k = "nNeighbors", dims.use = "useDims",
                        random.seed = "seed"))
    object <- recordCommand(object, dependencies = "RANN")
    H.norm <- getMatrix(object, "H.norm")
    if (is.null(H.norm)) {
        type <- " unnormalized "
        H.norm <- Reduce(cbind, getMatrix(object, "H"))
    } else type <- " quantile normalized "
    if (is.null(H.norm))
        stop("No factor loading ('H.norm' or 'H') found in `object`.")
    # Not transposing when cbind'ing becausing `t(NULL)` causes error
    if (type == " unnormalized ") H.norm <- t(H.norm)

    edgeOutPath <- paste0('edge_', sub('\\s', '_', Sys.time()), '.txt')
    edgeOutPath <- gsub("-", "", edgeOutPath)
    edgeOutPath <- gsub(":", "", edgeOutPath)

    if (!is.null(useDims)) H.norm <- H.norm[, useDims]

    if (isTRUE(verbose))
        .log("Louvain clustering on", type, "cell factor loadings...")
    knn <- RANN::nn2(H.norm, k = nNeighbors, eps = eps)
    snn <- ComputeSNN(knn$nn.idx, prune = prune)
    WriteEdgeFile(snn, edgeOutPath, display_progress = FALSE)
    clusts <- RunModularityClusteringCpp(
        snn,
        modularityFunction = 1,
        resolution = resolution,
        nRandomStarts = nRandomStarts,
        nIterations = nIterations,
        algorithm = 1,
        randomSeed = seed,
        printOutput = FALSE,
        edgefilename = edgeOutPath
    )
    names(clusts) <- colnames(object)
    rownames(snn) <- colnames(object)
    colnames(snn) <- colnames(object)
    # clusts must not be a factor at this point
    clusts <- groupSingletons(ids = clusts, SNN = snn,
                              groupSingletons = groupSingletons,
                              verbose = verbose)
    cellMeta(object, columns = clusterName, check = FALSE) <- factor(clusts)
    unlink(edgeOutPath)
    return(object)
}

# Group single cells that make up their own cluster in with the cluster they are
# most connected to. (Adopted from Seurat v3)
#
# @param ids Named vector of cluster ids
# @param SNN SNN graph used in clustering
# @param groupSingletons Group singletons into nearest cluster (TRUE by
# default). If FALSE, assign all singletons to a "singleton" group.
# @param verbose Print message
# @return Returns updated cluster assignment (vector) with all singletons merged
# with most connected cluster
groupSingletons <- function(
        ids,
        SNN,
        groupSingletons = TRUE,
        verbose = FALSE
) {
    # identify singletons
    singletons <- names(which(table(ids) == 1))
    singletons <- intersect(unique(ids), singletons)
    if (!isTRUE(groupSingletons)) {
        ids[ids %in% singletons] <- "singleton"
        return(ids)
    }
    # calculate connectivity of singletons to other clusters, add singleton
    # to cluster it is most connected to
    clusterNames <- as.character(unique(ids))
    clusterNames <- setdiff(clusterNames, singletons)
    connectivity <- vector(mode = "numeric", length = length(clusterNames))
    names(connectivity) <- clusterNames
    for (i in singletons) {
        i.cells <- names(which(ids == i))
        for (j in clusterNames) {
            j.cells <- names(which(ids == j))
            subSNN <- SNN[i.cells, j.cells]
            # to match previous behavior, random seed being set in WhichCells
            set.seed(1)
            if (is.object(subSNN))
                connectivity[j] <- sum(subSNN) / (nrow(subSNN) * ncol(subSNN))
            else
                connectivity[j] <- mean(subSNN)
        }
        m <- max(connectivity, na.rm = TRUE)
        mi <- which(connectivity == m, arr.ind = TRUE)
        closest_cluster <- sample(names(connectivity[mi]), 1)
        ids[i.cells] <- closest_cluster
    }
    if (length(singletons) > 0 && isTRUE(verbose))
        .log(length(singletons), " singletons identified. ",
             length(unique(ids)), " final clusters.")
    return(ids)
}

mapAnnotation <- function(lig, original, mapping, newName) {
    if (length(original) == 1) original <- lig[[original]]
    else if (length(original) != ncol(lig))
        stop("Wrong length from `original`.")

    newAnn <- rep(NA, ncol(lig))
    for (origClust in names(mapping)) {
        newAnn[original == origClust] <- mapping[[origClust]]
    }
    newAnn <- factor(newAnn)
    cellMeta(lig, columns = newName, check = FALSE) <- newAnn
    return(lig)
}
