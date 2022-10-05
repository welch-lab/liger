#' @importFrom Matrix colSums rowSums rowMeans t sparseMatrix
#' @importFrom grDevices heat.colors
#' @importFrom methods as is
#' @importFrom rlang .data
NULL

# Utility functions for iiger methods. Some published, some not.

# Helper function for plotFactors
fplot = function(tsne,NMFfactor,title,cols.use=heat.colors(10),pt.size=0.7,pch.use=20) {
  data.cut=as.numeric(as.factor(cut(as.numeric(NMFfactor),breaks=length(cols.use))))
  data.col=rev(cols.use)[data.cut]
  plot(tsne[,1],tsne[,2],col=data.col,cex=pt.size,pch=pch.use,main=title)
  
}

# Binds list of matrices row-wise (vertical stack)
rbindlist = function(mat_list)
{
  return(do.call(rbind,mat_list))
}

# helper function for calculating KL divergence from uniform distribution
# (related to Shannon entropy) for factorization
#' @importFrom methods is
kl_divergence_uniform = function(object, Hs=NULL)
{
  if (is.null(Hs)) {Hs = object@H}
  n_cells = sum(sapply(Hs, nrow))
  n_factors = ncol(Hs[[1]])
  dataset_list = list()
  for (i in 1:length(Hs)) {
    scaled = scale(Hs[[i]], center=FALSE, scale=TRUE)
    
    inflated = t(apply(scaled, 1, function(x) {
      replace(x, x == 0, 1e-20)
    }))
    inflated = inflated/rowSums(inflated)
    divs = apply(inflated, 1, function(x) {log2(n_factors) + sum(log2(x) * x)})
    dataset_list[[i]] = divs
  }
  return(dataset_list)
}

# Function takes in a list of DGEs, with gene rownames and cell colnames, 
# and merges them into a single DGE.
# Also adds library.names to cell.names if expected to be overlap (common with 10X barcodes)
MergeSparseDataAll <- function(datalist, library.names = NULL) {
  
  # Use summary to convert the sparse matrices into three-column indexes where i are the
  # row numbers, j are the column numbers, and x are the nonzero entries
  col_offset <- 0
  allGenes <- unique(unlist(lapply(datalist, rownames)))
  allCells <- c()
  for (i in 1:length(datalist)) {
    curr <- datalist[[i]]
    curr_s <- summary(curr)
    
    # Now, alter the indexes so that the two 3-column matrices can be properly merged.
    # First, make the current and full column numbers non-overlapping.
    curr_s[, 2] <- curr_s[, 2] + col_offset
    
    # Update full cell names
    if (!is.null(library.names)) {
      cellnames <- paste0(library.names[i], "_", colnames(curr))
    } else {
      cellnames <- colnames(curr)
    }
    allCells <- c(allCells, cellnames)
    
    # Next, change the row (gene) indexes so that they index on the union of the gene sets,
    # so that proper merging can occur.
    idx <- match(rownames(curr), allGenes)
    newgenescurr <- idx[curr_s[, 1]]
    curr_s[, 1] <- newgenescurr
    
    # Now bind the altered 3-column matrices together, and convert into a single sparse matrix.
    if (!exists("full_mat")) {
      full_mat <- curr_s
    } else {
      full_mat <- rbind(full_mat, curr_s)
    }
    col_offset <- length(allCells)
  }
  M <- sparseMatrix(
    i = full_mat[, 1],
    j = full_mat[, 2],
    x = full_mat[, 3],
    dims = c(
      length(allGenes),
      length(allCells)
    ),
    dimnames = list(
      allGenes,
      allCells
    )
  )
  return(M)
}

# Perform fast and memory-efficient normalization operation on sparse matrix data.
# param A Sparse matrix DGE.

Matrix.column_norm <- function(A){
  if (class(A)[1] == "dgTMatrix") {
    temp = summary(A)
    A = sparseMatrix(i=temp[,1],j=temp[,2],x=temp[,3])
  }
  A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
  return(A)
}

# Variance for sparse matrices
sparse.var = function(x){
  rms <- rowMeans(x)
  rowSums((x-rms)^2)/(dim(x)[2]-1)
}

# Transposition for sparse matrices
sparse.transpose = function(x){
  h = summary(x)
  sparseMatrix(i = h[,2],j=h[,1],x=h[,3])
  
}

# After running modularity clustering, assign singleton communities to the mode of the cluster
# assignments of the within-dataset neighbors
assign.singletons <- function(object, idents, k.use = 15, center = FALSE) {
  if (!inherits(x = object, what = 'liger')) {
    stop("'asign.singletons' expects an object of class 'liger'")
  }
  return(assign.singletons.list(
    object = object@H,
    idents = idents,
    k.use = k.use,
    center = center
  ))
}

#' @importFrom FNN get.knnx
#
assign.singletons.list <- function(object, idents, k.use = 15, center = FALSE) {
  if (!is.list(x = object) || !all(sapply(X = object, FUN = is.matrix))) {
    stop("'assign.singletons.list' expects a list of matrices")
  }
  message('Assigning singletons')
  singleton.clusters <- names(x = table(idents))[which(x = table(idents) == 1)]
  singleton.cells <- names(x = idents)[which(x = idents %in% singleton.clusters)]
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
        H = t(x = apply(X = H, MARGIN = 1, FUN = ifelse(test = center, yes = scale, no = scaleL2norm)))
        cells.use <- singleton.list[[d]]
        if (length(x = cells.use) > 0) {
          knn.idx <-  get.knnx(
            data = H,
            query = matrix(H[cells.use, ], ncol = ncol(x = H)),
            k = k.use,
            algorithm = "CR"
          )$nn.index
          o <- sapply(
            X = 1:length(x = cells.use),
            FUN = function(i) {
              ids <-  idents[knn.idx[i,]]
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
SLMCluster<-function(edge,prune.thresh=0.2,nstart=100,iter.max=10,algorithm=1,R=1, 
                     modularity=1, ModularityJarFile="",random.seed=1, 
                     id.number=NULL, print.mod=F) {
  
  # Prepare data for modularity based clustering
  edge = edge[which(edge[,3]>prune.thresh),]
  
  message("making edge file.")
  edge_file <- paste0("edge", id.number, fileext=".txt")
  # Make sure no scientific notation written to edge file

  # restore default settings when the current function exits
  init_option <- options() 
  on.exit(options(init_option))

  saveScipen=options(scipen=10000)[[1]]
  write.table(x=edge,file=edge_file,sep="\t",row.names=F,col.names=F)
  options(scipen=saveScipen)
  output_file <- tempfile("louvain.out", fileext=".txt")
  if (is.null(random.seed)) {
    # NULL random.seed disallowed for this program.
    random.seed = 0
  }
  liger.dir <- system.file(package = "rliger")
  ModularityJarFile <- paste0(liger.dir, "/java/ModularityOptimizer.jar")
  command <- paste("java -jar", ModularityJarFile, edge_file, output_file, modularity, R, 
                   algorithm, nstart,
                   iter.max, random.seed, as.numeric(print.mod), sep = " ")
  message ("Starting SLM")
  ret = system(command, wait = TRUE)
  if (ret != 0) {
    stop(paste0(ret, " exit status from ", command))
  }
  unlink(edge_file)
  ident.use <- factor(read.table(file = output_file, header = FALSE, sep = "\t")[, 1])
  
  return(ident.use)
}

# Scale to unit vector (scale by l2-norm)
scaleL2norm <- function(x) {
  return(x / sqrt(sum(x^2)))
}

# get mode of identities 
getMode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# utility function for seuratToLiger function
# Compares colnames in reference matrix1 and adds back any missing 
# column names to matrix.subset as rows 
# Set transpose = TRUE if rownames of matrix1 should be referenced
addMissingCells <- function(matrix1, matrix.subset, transpose = F) {
  if (transpose) {
    matrix1 <- t(matrix1)
  }
  if (ncol(matrix1) != nrow(matrix.subset)) {
    extra <- matrix(NA, nrow = ncol(matrix1) - nrow(matrix.subset),
                    ncol = ncol(matrix.subset))
    colnames(extra) <- colnames(matrix.subset)
    rownames(extra) <- setdiff(colnames(matrix1), rownames(matrix.subset))
    matrix.subset <- rbind(matrix.subset, extra)
    # get original order
    matrix.subset <- matrix.subset[colnames(matrix1), ]
  }
  return(matrix.subset)
}

#' Get gene expression values from list of expression matrices.
#'
#' @description
#' Returns single vector of gene values across all datasets in list provided. Data can be in raw, 
#' normalized or scaled form. If matrices are in cell x gene format, set use.cols = TRUE. 
#'
#' @param list List of gene x cell (or cell x gene) matrices
#' @param gene Gene for which to return values (if gene is not found in appropriate dimnames will
#'   return vector of NA).
#' @param use.cols Whether to query columns for desired gene (set to TRUE if matrices are cell x 
#'   gene) (default FALSE).
#' @param methylation.indices Indices of datasets with methylation data (never log2scaled) 
#'   (default NULL).
#' @param log2scale Whether to log2+1 scale (with multiplicative factor) values (default FALSE).
#' @param scale.factor Scale factor to use with log2 scaling (default 10000).
#'
#' @return Plots to console (1-2 pages per factor)
#' @export
#' @examples
#' \dontrun{
#' # liger object with factorization complete
#' # ligerex
#' gene_values <- getGeneValues(ligerex@raw.data, 'MALAT1')
#' }

getGeneValues <- function(list, gene, use.cols = FALSE, methylation.indices = NULL, log2scale = FALSE,
                          scale.factor = 10000) {
  gene_vals <- unlist(lapply(seq_along(list), function(i) {
    mtx <- unname(list)[[i]]
    if (use.cols) {
      mtx <- t(mtx)
    } 
    if (gene %in% rownames(mtx)) {
      gene_vals_int <- mtx[gene, ]
    } else {
      gene_vals_int <- rep(0, ncol(mtx))
      names(gene_vals_int) <- colnames(mtx)
    }
    if (log2scale & !(i %in% methylation.indices)) {
      gene_vals_int <- log2(scale.factor * gene_vals_int + 1)
    } 
    return(gene_vals_int)
  }),
  use.names = TRUE)
  return(gene_vals)
}

# helper function for refining clusers by KNN
# related to function quantile_norm
refine_clusts_knn = function(H,clusts,k,eps=0.1)
{
  for (i in 1:length(H))
  {
    clusts_H = clusts[rownames(H[[i]])]
    H_knn = RANN::nn2(H[[i]],eps=eps,k=k,searchtype="standard")
    #for (j in 1:length(clusts_H))
    #{
    #  clusts_H[j] = names(which.max(table(clusts_H[H_knn$nn.idx[j,]])))
    #}
    #clusts[rownames(H[[i]])] = clusts_H
    new_clusts = cluster_vote(H_knn$nn.idx,clusts_H)
    clusts[rownames(H[[i]])] = new_clusts
  }
  return(clusts)
}


################################## For fast Wilcoxon test ################################
# helper function for wilcoxon tests on general variables like matrix and dgCMatrix
# related to function runWilcoxon
#' @importFrom stats p.adjust
wilcoxauc <- function(X, y, groups_use=NULL, verbose=TRUE) {
  ## Check and possibly correct input values
  if (is(X, 'dgeMatrix')) X <- as.matrix(X)
  if (is(X, 'data.frame')) X <- as.matrix(X)
  # if (is(X, 'DataFrame')) X <- as.matrix(X)
  # if (is(X, 'data.table')) X <- as.matrix(X)
  if (is(X, 'dgTMatrix')) X <- as(X, 'dgCMatrix')
  if (is(X, 'TsparseMatrix')) X <- as(X, 'dgCMatrix')
  if (ncol(X) != length(y)) stop("number of columns of X does not
                                 match length of y")
  if (!is.null(groups_use)) {
    idx_use <- which(y %in% intersect(groups_use, y))
    y <- y[idx_use]
    X <- X[, idx_use]
  }
  
  
  y <- factor(y)
  idx_use <- which(!is.na(y))
  if (length(idx_use) < length(y)) {
    y <- y[idx_use]
    X <- X[, idx_use]
    if (verbose) 
      message('Removing NA values from labels')        
  }
  
  group.size <- as.numeric(table(y))
  if (length(group.size[group.size > 0]) < 2) {
    stop('Must have at least 2 groups defined.')
  }
  
  #     features_use <- which(apply(!is.na(X), 1, all))
  #     if (verbose & length(features_use) < nrow(X)) {
  #         message('Removing features with NA values')
  #     }
  #     X <- X[features_use, ]
  if (is.null(row.names(X))) {
    row.names(X) <- paste0('Feature', seq_len(nrow(X)))
  }
  
  ## Compute primary statistics
  group.size <- as.numeric(table(y))
  n1n2 <- group.size * (ncol(X) - group.size)
  if (is(X, 'dgCMatrix')) {
    rank_res <- rank_matrix(Matrix::t(X))        
  } else {
    rank_res <- rank_matrix(X)
  }
  
  ustat <- compute_ustat(rank_res$X_ranked, y, n1n2, group.size) 
  auc <- t(ustat / n1n2)
  pvals <- compute_pval(ustat, rank_res$ties, ncol(X), n1n2) 
  fdr <- apply(pvals, 2, function(x) p.adjust(x, 'BH'))
  
  ### Auxiliary Statistics (AvgExpr, PctIn, LFC, etc)
  group_sums <- sumGroups(X, y, 1)
  group_nnz <- nnzeroGroups(X, y, 1)
  group_pct <- sweep(group_nnz, 1, as.numeric(table(y)), "/") %>% t()
  group_pct_out <- -group_nnz %>% 
    sweep(2, colSums(group_nnz) , "+") %>% 
    sweep(1, as.numeric(length(y) - table(y)), "/") %>% t()
  group_means <- sweep(group_sums, 1, as.numeric(table(y)), "/") %>% t()
  cs <- colSums(group_sums)
  gs <- as.numeric(table(y))
  lfc <- Reduce(cbind, lapply(seq_len(length(levels(y))), function(g) {
    group_means[, g] - ((cs - group_sums[g, ]) / (length(y) - gs[g]))
  }))
  
  res_list <- list(auc = auc, 
                   pval = pvals,
                   padj = fdr, 
                   pct_in = 100 * group_pct, 
                   pct_out = 100 * group_pct_out,
                   avgExpr = group_means, 
                   statistic = t(ustat),
                   logFC = lfc)
  return(tidy_results(res_list, row.names(X), levels(y)))
}


tidy_results <- function(wide_res, features, groups) {
  res <- Reduce(cbind, lapply(wide_res, as.numeric)) %>% data.frame() 
  colnames(res) <- names(wide_res)
  res$feature <- rep(features, times = length(groups))
  res$group <- rep(groups, each = length(features))
  res %>% dplyr::select(
    .data$feature, 
    .data$group, 
    .data$avgExpr, 
    .data$logFC, 
    .data$statistic, 
    .data$auc, 
    .data$pval, 
    .data$padj, 
    .data$pct_in, 
    .data$pct_out
  )
}


compute_ustat <- function(Xr, cols, n1n2, group.size) {
  grs <- sumGroups(Xr, cols)
  
  if (is(Xr, 'dgCMatrix')) {
    gnz <- (group.size - nnzeroGroups(Xr, cols))
    zero.ranks <- (nrow(Xr) - diff(Xr@p) + 1) / 2
    ustat <- t((t(gnz) * zero.ranks)) + grs - group.size *
      (group.size + 1 ) / 2        
  } else {
    ustat <- grs - group.size * (group.size + 1 ) / 2
  }
  return(ustat)
}

#' @importFrom stats pnorm
compute_pval <- function(ustat, ties, N, n1n2) {
  z <- ustat - .5 * n1n2
  z <- z - sign(z) * .5
  .x1 <- N ^ 3 - N
  .x2 <- 1 / (12 * (N^2 - N))
  rhs <- lapply(ties, function(tvals) {
    (.x1 - sum(tvals ^ 3 - tvals)) * .x2
  }) %>% unlist
  usigma <- sqrt(matrix(n1n2, ncol = 1) %*% matrix(rhs, nrow = 1))
  z <- t(z / usigma)
  
  pvals <- matrix(2 * pnorm(-abs(as.numeric(z))), ncol = ncol(z))
  return(pvals)
}


#' rank_matrix
#' 
#' Utility function to rank columns of matrix
#' 
#' @param X feature by observation matrix. 
#' 
#' @return List with 2 items


rank_matrix <- function(X) {
  UseMethod('rank_matrix')
}

##' @rdname rank_matrix
rank_matrix.dgCMatrix <- function(X) {
  Xr <- Matrix(X, sparse = TRUE)
  ties <- cpp_rank_matrix_dgc(Xr@x, Xr@p, nrow(Xr), ncol(Xr))
  return(list(X_ranked = Xr, ties = ties))
}

##' @rdname rank_matrix
rank_matrix.matrix <- function(X) {
  cpp_rank_matrix_dense(X)
}

#' sumGroups
#' 
#' Utility function to sum over group labels
#' 
#' @param X matrix
#' @param y group labels
#' @param MARGIN whether observations are rows (=2) or columns (=1)
#' 
#' @return Matrix of groups by features

sumGroups <- function(X, y, MARGIN=2) {
  if (MARGIN == 2 & nrow(X) != length(y)) {
    stop('wrong dims')
  } else if (MARGIN == 1 & ncol(X) != length(y)) {
    stop('wrong dims') 
  }
  UseMethod('sumGroups')
}

##' @rdname sumGroups
sumGroups.dgCMatrix <- function(X, y, MARGIN=2) {
  if (MARGIN == 1) {
    cpp_sumGroups_dgc_T(X@x, X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1,
                        length(unique(y)))        
  } else {
    cpp_sumGroups_dgc(X@x, X@p, X@i, ncol(X), as.integer(y) - 1,
                      length(unique(y)))
  }
}

##' @rdname sumGroups
sumGroups.matrix <- function(X, y, MARGIN=2) {
  if (MARGIN == 1) {
    cpp_sumGroups_dense_T(X, as.integer(y) - 1, length(unique(y)))        
  } else {
    cpp_sumGroups_dense(X, as.integer(y) - 1, length(unique(y)))
  }
}



#' nnzeroGroups
#' 
#' Utility function to compute number of zeros-per-feature within group
#' 
#' @param X matrix
#' @param y group labels
#' @param MARGIN whether observations are rows (=2) or columns (=1)
#' 
#' @return Matrix of groups by features

nnzeroGroups <- function(X, y, MARGIN=2) {
  if (MARGIN == 2 & nrow(X) != length(y)) {
    stop('wrong dims')
  } else if (MARGIN == 1 & ncol(X) != length(y)) {
    stop('wrong dims')        
  }
  UseMethod('nnzeroGroups')
}

##' @rdname nnzeroGroups
nnzeroGroups.dgCMatrix <- function(X, y, MARGIN=2) {
  if (MARGIN == 1) {
    cpp_nnzeroGroups_dgc_T(X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1,
                           length(unique(y)))        
  } else {
    cpp_nnzeroGroups_dgc(X@p, X@i, ncol(X), as.integer(y) - 1,
                         length(unique(y)))
  }
}

##' @rdname nnzeroGroups
nnzeroGroups.matrix <- function(X, y, MARGIN=2) {
  if (MARGIN == 1) {        
    cpp_nnzeroGroups_dense_T(X, as.integer(y) - 1, length(unique(y)))        
  } else {
    cpp_nnzeroGroups_dense(X, as.integer(y) - 1, length(unique(y)))
  }
}


# helper function of nmf_hals
nonneg <- function(x, eps = 1e-16) {
  x[x < eps] <- eps
  return(x)
}

frobenius_prod <- function(X, Y) {
  sum(X * Y)
}

# Hierarchical alternating least squares for regular NMF
nmf_hals <- function(A, k, max_iters = 500, thresh = 1e-4, reps = 20, W0 = NULL, H0 = NULL) {
  m <- nrow(A)
  n <- ncol(A)
  if (is.null(W0)) {
    W0 <- matrix(abs(runif(m * k, min = 0, max = 2)), m, k)
  }
  if (is.null(H0)) {
    H0 <- matrix(abs(runif(n * k, min = 0, max = 2)), n, k)
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
  cat("\nConverged in", end_time - start_time, "seconds,", iters, "iterations. Objective:", obj, "\n")
  # boxplot(iter_times)

  return(list(W, H, cumsum(iter_times), objs))
}


# FIt-SNE helper function for calling fast_tsne from R
#
# Codes from Seurat (https://github.com/satijalab/seurat)
#
# Originally Based on Kluger Lab FIt-SNE v1.1.0 code on https://github.com/KlugerLab/FIt-SNE/blob/master/fast_tsne.R
# commit d2cf403 on Feb 8, 2019
#
#' @importFrom utils file_test
#
fftRtsne <- function(X,
                     dims = 2,
                     perplexity = 30,
                     theta = 0.5,
                     check_duplicates = TRUE,
                     max_iter = 1000,
                     fft_not_bh = TRUE,
                     ann_not_vptree = TRUE,
                     stop_early_exag_iter = 250,
                     exaggeration_factor = 12.0,
                     no_momentum_during_exag = FALSE,
                     start_late_exag_iter = -1.0,
                     late_exag_coeff = 1.0,
                     mom_switch_iter = 250,
                     momentum = 0.5,
                     final_momentum = 0.8,
                     learning_rate = 200,
                     n_trees = 50,
                     search_k = -1,
                     rand_seed = -1,
                     nterms = 3,
                     intervals_per_integer = 1,
                     min_num_intervals = 50,
                     K = -1,
                     sigma = -30,
                     initialization = NULL,
                     data_path = NULL,
                     result_path = NULL,
                     load_affinities = NULL,
                     fast_tsne_path = NULL,
                     nthreads = getOption("mc.cores", default = 1),
                     perplexity_list = NULL,
                     get_costs = FALSE,
                     df = 1.0,
                     ...) {
  if (is.null(fast_tsne_path)) {
    stop("Please pass in path to FIt-SNE directory as fitsne.path.")
  }
  else {
    if (.Platform$OS.type == "unix") {
      fast_tsne_path <- file.path(fast_tsne_path, "bin", "fast_tsne")
    } else {
      fast_tsne_path <- file.path(fast_tsne_path, "bin", "FItSNE.exe")
    }
  }

  if (is.null(x = data_path)) {
    data_path <- tempfile(pattern = "fftRtsne_data_", fileext = ".dat")
  }
  if (is.null(x = result_path)) {
    result_path <- tempfile(pattern = "fftRtsne_result_", fileext = ".dat")
  }

  fast_tsne_path <- normalizePath(path = fast_tsne_path)
  if (!file_test(op = "-x", x = fast_tsne_path)) {
    stop("fast_tsne_path '", fast_tsne_path, "' does not exist or is not executable")
  }
  # check fast_tsne version
  ft.out <- suppressWarnings(expr = system2(command = fast_tsne_path, stdout = TRUE))
  if (grepl(pattern = "= t-SNE v1.1", x = ft.out[1])) {
    version_number <- "1.1.0"
  } else if (grepl(pattern = "= t-SNE v1.0", x = ft.out[1])) {
    version_number <- "1.0"
  } else {
    message("First line of fast_tsne output is")
    message(ft.out[1])
    stop("Our FIt-SNE wrapper requires FIt-SNE v1.X.X, please install the appropriate version from github.com/KlugerLab/FIt-SNE and have fast_tsne_path point to it if it's not in your path")
  }
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    return(abs(x = x - round(x = x)) < tol)
  }
  if (version_number == "1.0" && df != 1.0) {
    stop("This version of FIt-SNE does not support df!=1. Please install the appropriate version from github.com/KlugerLab/FIt-SNE")
  }
  if (!is.numeric(x = theta) || (theta < 0.0) || (theta > 1.0)) {
    stop("Incorrect theta.")
  }
  if (nrow(x = X) - 1 < 3 * perplexity) {
    stop("Perplexity is too large.")
  }
  if (!is.matrix(x = X)) {
    stop("Input X is not a matrix")
  }
  if (!(max_iter > 0)) {
    stop("Incorrect number of iterations.")
  }
  if (!is.wholenumber(x = stop_early_exag_iter) || stop_early_exag_iter < 0) {
    stop("stop_early_exag_iter should be a positive integer")
  }
  if (!is.numeric(x = exaggeration_factor)) {
    stop("exaggeration_factor should be numeric")
  }
  if (!is.wholenumber(x = dims) || dims <= 0) {
    stop("Incorrect dimensionality.")
  }
  if (search_k == -1) {
    if (perplexity > 0) {
      search_k <- n_trees * perplexity * 3
    } else if (perplexity == 0) {
      search_k <- n_trees * max(perplexity_list) * 3
    } else {
      search_k <- n_trees * K * 3
    }
  }
  nbody_algo <- ifelse(test = fft_not_bh, yes = 2, no = 1)
  if (is.null(load_affinities)) {
    load_affinities <- 0
  } else {
    if (load_affinities == "load") {
      load_affinities <- 1
    } else if (load_affinities == "save") {
      load_affinities <- 2
    } else {
      load_affinities <- 0
    }
  }
  knn_algo <- ifelse(test = ann_not_vptree, yes = 1, no = 2)
  f <- file(description = data_path, open = "wb")
  n <- nrow(x = X)
  D <- ncol(x = X)
  writeBin(object = as.integer(x = n), con = f, size = 4)
  writeBin(object = as.integer(x = D), con = f, size = 4)
  writeBin(object = as.numeric(x = theta), con = f, size = 8) # theta
  writeBin(object = as.numeric(x = perplexity), con = f, size = 8) # theta
  if (perplexity == 0) {
    writeBin(object = as.integer(x = length(x = perplexity_list)), con = f, size = 4)
    writeBin(object = perplexity_list, con = f)
  }
  writeBin(object = as.integer(x = dims), con = f, size = 4) # theta
  writeBin(object = as.integer(x = max_iter), con = f, size = 4)
  writeBin(object = as.integer(x = stop_early_exag_iter), con = f, size = 4)
  writeBin(object = as.integer(x = mom_switch_iter), con = f, size = 4)
  writeBin(object = as.numeric(x = momentum), con = f, size = 8)
  writeBin(object = as.numeric(x = final_momentum), con = f, size = 8)
  writeBin(object = as.numeric(x = learning_rate), con = f, size = 8)
  writeBin(object = as.integer(x = K), con = f, size = 4) # K
  writeBin(object = as.numeric(x = sigma), con = f, size = 8) # sigma
  writeBin(object = as.integer(x = nbody_algo), con = f, size = 4) # not barnes hut
  writeBin(object = as.integer(x = knn_algo), con = f, size = 4)
  writeBin(object = as.numeric(x = exaggeration_factor), con = f, size = 8) # compexag
  writeBin(object = as.integer(x = no_momentum_during_exag), con = f, size = 4)
  writeBin(object = as.integer(x = n_trees), con = f, size = 4)
  writeBin(object = as.integer(x = search_k), con = f, size = 4)
  writeBin(object = as.integer(x = start_late_exag_iter), con = f, size = 4)
  writeBin(object = as.numeric(x = late_exag_coeff), con = f, size = 8)
  writeBin(object = as.integer(x = nterms), con = f, size = 4)
  writeBin(object = as.numeric(x = intervals_per_integer), con = f, size = 8)
  writeBin(object = as.integer(x = min_num_intervals), con = f, size = 4)
  tX <- c(t(X))
  writeBin(object = tX, con = f)
  writeBin(object = as.integer(x = rand_seed), con = f, size = 4)
  if (version_number != "1.0") {
    writeBin(object = as.numeric(x = df), con = f, size = 8)
  }
  writeBin(object = as.integer(x = load_affinities), con = f, size = 4)
  if (!is.null(x = initialization)) {
    writeBin(object = c(t(x = initialization)), con = f)
  }
  close(con = f)
  if (version_number == "1.0") {
    flag <- system2(
      command = fast_tsne_path,
      args = c(data_path, result_path, nthreads)
    )
  } else {
    flag <- system2(
      command = fast_tsne_path,
      args = c(version_number, data_path, result_path, nthreads)
    )
  }
  if (flag != 0) {
    stop("tsne call failed")
  }
  f <- file(description = result_path, open = "rb")
  n <- readBin(con = f, what = integer(), n = 1, size = 4)
  d <- readBin(con = f, what = integer(), n = 1, size = 4)
  Y <- readBin(con = f, what = numeric(), n = n * d)
  Y <- t(x = matrix(Y, nrow = d))
  if (get_costs) {
    tmp <- readBin(con = f, what = integer(), n = 1, size = 4)
    costs <- readBin(con = f, what = numeric(), n = max_iter, size = 8)
    Yout <- list(Y = Y, costs = costs)
  } else {
    Yout <- Y
  }
  close(con = f)
  file.remove(data_path)
  file.remove(result_path)
  return(Yout)
}
