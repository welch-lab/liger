#' @title Deprecated functions in package \pkg{rliger}.
#' @description The functions listed below are deprecated and will be defunct in
#'   the near future. When possible, alternative functions with similar
#'   functionality or a replacement are also mentioned. Help pages for
#'   deprecated functions are available at \code{help("<function>-deprecated")}.
#' @name rliger-deprecated
#' @keywords internal
NULL


# These are deprecated functions likely to be removed in future versions.
# Documentation for these functions is incomplete.

#' Quantile align (normalize) factor loadings
#'
#' This is a deprecated function. Calling 'quantileNorm' instead.
#'
#' This process builds a shared factor neighborhood graph to jointly cluster cells, then quantile
#' normalizes corresponding clusters.
#'
#' The first step, building the shared factor neighborhood graph, is performed in SNF(), and
#' produces a graph representation where edge weights between cells (across all datasets)
#' correspond to their similarity in the shared factor neighborhood space. An important parameter
#' here is knn_k, the number of neighbors used to build the shared factor space (see SNF()). Afterwards,
#' modularity-based community detection is performed on this graph (Louvain clustering) in order
#' to identify shared clusters across datasets. The method was first developed by Waltman and van Eck
#' (2013) and source code is available at http://www.ludowaltman.nl/slm/. The most important parameter
#' here is resolution, which corresponds to the number of communities detected.
#'
#' Next we perform quantile alignment for each dataset, factor, and cluster (by
#' stretching/compressing datasets' quantiles to better match those of the reference dataset). These
#' aligned factor loadings are combined into a single matrix and returned as H.norm.
#'
#' @param object \code{liger} object. Should run optimizeALS before calling.
#' @param knn_k Number of nearest neighbors for within-dataset knn graph (default 20).
#' @param k2 Horizon parameter for shared nearest factor graph. Distances to all but the k2 nearest
#'   neighbors are set to 0 (cuts down on memory usage for very large graphs). (default 500)
#' @param prune.thresh Minimum allowed edge weight. Any edges below this are removed (given weight
#'  0) (default 0.2)
#' @param ref_dataset Name of dataset to use as a "reference" for normalization. By default,
#'   the dataset with the largest number of cells is used.
#' @param min_cells Minimum number of cells to consider a cluster shared across datasets (default 2)
#' @param quantiles Number of quantiles to use for quantile normalization (default 50).
#' @param nstart Number of times to perform Louvain community detection with different random
#'   starts (default 10).
#' @param resolution Controls the number of communities detected. Higher resolution -> more
#'   communities. (default 1)
#' @param dims.use Indices of factors to use for shared nearest factor determination (default
#'   1:ncol(H[[1]])).
#' @param dist.use Distance metric to use in calculating nearest neighbors (default "CR").
#' @param center Centers the data when scaling factors (useful for less sparse modalities like
#'   methylation data). (default FALSE)
#' @param small.clust.thresh Extracts small clusters loading highly on single factor with fewer
#'   cells than this before regular alignment (default 0 -- no small cluster extraction).
#' @param id.number Number to use for identifying edge file (when running in parallel)
#'   (generates random value by default).
#' @param print.mod Print modularity output from clustering algorithm (default FALSE).
#' @param print.align.summary Print summary of clusters which did not align normally (default FALSE).
#'
#' @return \code{liger} object with H.norm and cluster slots set.
#' @export
#'
#' @examples
#' \dontrun{
#' # liger object, factorization complete
#' ligerex
#' # do basic quantile alignment
#' ligerex <- quantileAlignSNF(ligerex)
#' # higher resolution for more clusters (note that SNF is conserved)
#' ligerex <- quantileAlignSNF(ligerex, resolution = 1.2)
#' # change knn_k for more fine-grained local clustering
#' ligerex <- quantileAlignSNF(ligerex, knn_k = 15, resolution = 1.2)
#' }
#'
quantileAlignSNF <- function(
                             object,
                             knn_k = 20,
                             k2 = 500,
                             prune.thresh = 0.2,
                             ref_dataset = NULL,
                             min_cells = 20,
                             quantiles = 50,
                             nstart = 10,
                             resolution = 1,
                             dims.use = 1:ncol(x = object@H[[1]]),
                             dist.use = "CR",
                             center = FALSE,
                             small.clust.thresh = 0,
                             id.number = NULL,
                             print.mod = FALSE,
                             print.align.summary = FALSE) {
  .Deprecated(
    msg = paste(
      "This is a deprecated function. Calling 'quantileNorm' instead.",
      "Note that not all parameters can be passed to 'quantileNorm'.",
      "It's suggested to run 'louvainCluster' subsequently as well."
    )
  )
  quantileNorm(object, quantiles = quantiles, reference = ref_dataset,
               minCells = min_cells, nNeighbors = knn_k, useDims = dims.use,
               center = center, maxSample = 1000, eps = 0.9, refineKNN = TRUE)
}





# Variance for sparse matrices
sparse.var = function(x) {
  rms <- rowMeans(x)
  rowSums((x - rms) ^ 2) / (dim(x)[2] - 1)
}

# Transposition for sparse matrices
sparse.transpose = function(x) {
  h = summary(x)
  Matrix::sparseMatrix(i = h[, 2], j = h[, 1], x = h[, 3])
}

# After running modularity clustering, assign singleton communities to the mode of the cluster
# assignments of the within-dataset neighbors
assign.singletons <-
  function(object,
           idents,
           k.use = 15,
           center = FALSE) {
    if (!inherits(x = object, what = 'liger')) {
      stop("'asign.singletons' expects an object of class 'liger'")
    }
    return(
      assign.singletons.list(
        object = object@H,
        idents = idents,
        k.use = k.use,
        center = center
      )
    )
  }

assign.singletons.list <-
  function(object,
           idents,
           k.use = 15,
           center = FALSE) {
    if (!is.list(x = object) ||
        !all(sapply(X = object, FUN = is.matrix))) {
      stop("'assign.singletons.list' expects a list of matrices")
    }
    message('Assigning singletons')
    singleton.clusters <-
      names(x = table(idents))[which(x = table(idents) == 1)]
    singleton.cells <-
      names(x = idents)[which(x = idents %in% singleton.clusters)]
    if (length(x = singleton.cells) > 0) {
      singleton.list <- lapply(
        X = object,
        FUN = function(H) {
          return(intersect(x = rownames(x = H), y = singleton.cells))
        }
      )
      out <- unlist(x = lapply(
        X = 1:length(x = object),
        FUN = function(d) {
          H = object[[d]]
          H = t(x = apply(
            X = H,
            MARGIN = 1,
            FUN = ifelse(
              test = center,
              yes = scale,
              no = scaleL2norm
            )
          ))
          cells.use <- singleton.list[[d]]
          if (length(x = cells.use) > 0) {
            knn.idx <-  RANN::nn2(
              data = H,
              query = matrix(H[cells.use,], ncol = ncol(x = H)),
              k = k.use,
            )$nn.idx
            o <- sapply(
              X = 1:length(x = cells.use),
              FUN = function(i) {
                ids <-  idents[knn.idx[i, ]]
                return(getMode(x = ids[which(x = !(ids %in% singleton.clusters))]))
              }
            )
            return(o)
          } else {
            return(NULL)
          }
        }
      ))
      names(x = out) <- unlist(x = singleton.list)
      idx <- names(x = idents) %in% names(x = out)
      idents[idx] <- out
      idents <- droplevels(idents)
    }
    return(idents)
  }

# Run modularity based clustering on edge file
SLMCluster <-
  function(edge,
           prune.thresh = 0.2,
           nstart = 100,
           iter.max = 10,
           algorithm = 1,
           R = 1,
           modularity = 1,
           ModularityJarFile = "",
           random.seed = 1,
           id.number = NULL,
           print.mod = F) {
    # Prepare data for modularity based clustering
    edge = edge[which(edge[, 3] > prune.thresh), ]

    message("making edge file.")
    edge_file <- paste0("edge", id.number, fileext = ".txt")
    # Make sure no scientific notation written to edge file

    # restore default settings when the current function exits
    init_option <- options()
    on.exit(options(init_option))

    saveScipen = options(scipen = 10000)[[1]]
    utils::write.table(
      x = edge,
      file = edge_file,
      sep = "\t",
      row.names = F,
      col.names = F
    )
    options(scipen = saveScipen)
    output_file <- tempfile("louvain.out", fileext = ".txt")
    if (is.null(random.seed)) {
      # NULL random.seed disallowed for this program.
      random.seed = 0
    }
    liger.dir <- "."#system.file(package = "rliger")
    ModularityJarFile <-
      paste0(liger.dir, "/java/ModularityOptimizer.jar")
    command <-
      paste(
        "java -jar",
        ModularityJarFile,
        edge_file,
        output_file,
        modularity,
        R,
        algorithm,
        nstart,
        iter.max,
        random.seed,
        as.numeric(print.mod),
        sep = " "
      )
    message("Starting SLM")
    ret = system(command, wait = TRUE)
    if (ret != 0) {
      stop(paste0(ret, " exit status from ", command))
    }
    unlink(edge_file)
    ident.use <-
      factor(utils::read.table(
        file = output_file,
        header = FALSE,
        sep = "\t"
      )[, 1])

    return(ident.use)
  }


# Scale to unit vector (scale by l2-norm)
scaleL2norm <- function(x) {
  return(x / sqrt(sum(x ^ 2)))
}

# get mode of identities
getMode <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x = x[!is.na(x)]
  }
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# Helper function for plotFactors
fplot <- function(tsne,
                 NMFfactor,
                 title,
                 cols.use = grDevices::heat.colors(10),
                 pt.size = 0.7,
                 pch.use = 20) {
  data.cut = as.numeric(as.factor(cut(
    as.numeric(NMFfactor), breaks = length(cols.use)
  )))
  data.col = rev(cols.use)[data.cut]
  plot(
    tsne[, 1],
    tsne[, 2],
    col = data.col,
    cex = pt.size,
    pch = pch.use,
    main = title
  )

}
