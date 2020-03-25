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
    scaled = scale(Hs[[i]], center=F, scale=T)
    
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

#' Perform fast and memory-efficient normalization operation on sparse matrix data.
#'
#' @param A Sparse matrix DGE.
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' ligerex = createLiger(list(Y,Z))
#' ligerex@var.genes = c(1,2,3,4)
#' ligerex = scaleNotCenter(ligerex)
#' }
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
  print('Assigning singletons')
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
  
  print("making edge file.")
  edge_file <- paste0("edge", id.number, fileext=".txt")
  # Make sure no scientific notation written to edge file
  saveScipen=options(scipen=10000)[[1]]
  write.table(x=edge,file=edge_file,sep="\t",row.names=F,col.names=F)
  options(scipen=saveScipen)
  output_file <- tempfile("louvain.out", fileext=".txt")
  if (is.null(random.seed)) {
    # NULL random.seed disallowed for this program.
    random.seed = 0
  }
  liger.dir <- system.file(package = "liger")
  ModularityJarFile <- paste0(liger.dir, "/java/ModularityOptimizer.jar")
  command <- paste("java -jar", ModularityJarFile, edge_file, output_file, modularity, R, 
                   algorithm, nstart,
                   iter.max, random.seed, as.numeric(print.mod), sep = " ")
  print ("Starting SLM")
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
# Set transpose = T if rownames of matrix1 should be referenced
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
#' normalized or scaled form. If matrices are in cell x gene format, set use.cols = T. 
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
#'  # liger object with factorization complete
#' ligerex
#' gene_values <- getGeneValues(ligerex@raw.data, 'MALAT1')
#' }

getGeneValues <- function(list, gene, use.cols = F, methylation.indices = NULL, log2scale = F,
                          scale.factor = 10000) {
  gene_vals <- unlist(lapply(seq_along(list), function(i) {
    mtx <- unname(list)[[i]]
    if (use.cols) {
      mtx <- t(mtx)
    } 
    if (gene %in% rownames(mtx)) {
      gene_vals_int <- mtx[gene, ]
    } else {
      gene_vals_int <- rep(list(0), ncol(mtx))
      names(gene_vals_int) <- colnames(mtx)
    }
    if (log2scale & !(i %in% methylation.indices)) {
      gene_vals_int <- log2(scale.factor * gene_vals_int + 1)
    } 
    return(gene_vals_int)
  }),
  use.names = T)
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



#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @importFrom dplyr %>%
#' @examples
#' \dontrun{
#' x <- 5 %>% sum(10)
#' }
#' @return return value of rhs function. 
NULL


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
#' @examples
#' \dontrun{
#' data(exprs)
#' rank_res <- rank_matrix(exprs)
#' }
#' @return List with 2 items
#' \itemize{
#' \item X_ranked - matrix of entry ranks
#' \item ties - list of tied group sizes
#' }
rank_matrix <- function(X) {
  UseMethod('rank_matrix')
}

#' @rdname rank_matrix
rank_matrix.dgCMatrix <- function(X) {
  Xr <- Matrix(X, sparse = TRUE)
  ties <- cpp_rank_matrix_dgc(Xr@x, Xr@p, nrow(Xr), ncol(Xr))
  return(list(X_ranked = Xr, ties = ties))
}

#' @rdname rank_matrix
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
#' @examples
#' \dontrun{
#' data(exprs)
#' data(y)
#' sumGroups_res <- sumGroups(exprs, y, 1)
#' sumGroups_res <- sumGroups(t(exprs), y, 2)
#' }
#' @return Matrix of groups by features
sumGroups <- function(X, y, MARGIN=2) {
  if (MARGIN == 2 & nrow(X) != length(y)) {
    stop('wrong dims')
  } else if (MARGIN == 1 & ncol(X) != length(y)) {
    stop('wrong dims') 
  }
  UseMethod('sumGroups')
}

#' @rdname sumGroups
sumGroups.dgCMatrix <- function(X, y, MARGIN=2) {
  if (MARGIN == 1) {
    cpp_sumGroups_dgc_T(X@x, X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1,
                        length(unique(y)))        
  } else {
    cpp_sumGroups_dgc(X@x, X@p, X@i, ncol(X), as.integer(y) - 1,
                      length(unique(y)))
  }
}

#' @rdname sumGroups
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
#' @examples
#' 
#' \dontrun{
#' data(exprs)
#' data(y)
#' nnz_res <- nnzeroGroups(exprs, y, 1)
#' nnz_res <- nnzeroGroups(t(exprs), y, 2)
#' }
#' @return Matrix of groups by features
nnzeroGroups <- function(X, y, MARGIN=2) {
  if (MARGIN == 2 & nrow(X) != length(y)) {
    stop('wrong dims')
  } else if (MARGIN == 1 & ncol(X) != length(y)) {
    stop('wrong dims')        
  }
  UseMethod('nnzeroGroups')
}

#' @rdname nnzeroGroups
nnzeroGroups.dgCMatrix <- function(X, y, MARGIN=2) {
  if (MARGIN == 1) {
    cpp_nnzeroGroups_dgc_T(X@p, X@i, ncol(X), nrow(X), as.integer(y) - 1,
                           length(unique(y)))        
  } else {
    cpp_nnzeroGroups_dgc(X@p, X@i, ncol(X), as.integer(y) - 1,
                         length(unique(y)))
  }
}

#' @rdname nnzeroGroups
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
