#' @importFrom plyr rbind.fill.matrix
NULL

#' The Analogizer Class
#'
#' The analogizer object is created from two or more single cell datasets.
#' The class provides functions for data preprocessing, integrative analysis, and visualization.
#'
#' @slot raw.data List of raw data matrices, one per experiment (genes by cells)
#' @slot norm.data List of normalized matrices
#' @slot scale.data List of scaled matrices
#' @slot var.genes Subset of informative genes shared across datasets to be used in matrix factorization
#' @slot H Cell loading factors (one matrix per dataset)
#' @slot H.norm Normalized cell loading factors (one matrix per dataset)
#' @slot W Shared gene loading factors
#' @slot V Dataset-specific gene loading factors
#' @slot tsne.coords Matrix of 2D coordinates obtained from running t-SNE on H.norm
#' @slot parameters List of parameters used in analysis
#'
#' @name analogizer
#' @rdname analogizer
#' @aliases analogizer-class
#' @exportClass analogizer
#' @importFrom Rcpp evalCpp
#' @useDynLib Analogizer

analogizer <- methods::setClass(
  "analogizer",
  slots = c(
    raw.data = "list",
    norm.data = "list",
    scale.data = "list",
    var.genes = "vector",
    H = "list",
    H.norm = "list",
    W = "matrix",
    V = "list",
    tsne.coords = "matrix",
    parameters = "list"
  )
)

Analogizer <- function(raw.data) {
  object <- new(
    Class = "analogizer",
    raw.data = raw.data
  )
  return(object)
}

#' Select a subset of informative genes from each dataset and combine them.
#'
#' @param object analogizer object. Should have already called normalize.
#' @param alphathresh alpha threshold
#' @param varthresh variance threshold
#' @param cex.use point size for plot
#' @param combine How to combine variable genes across experiments. Either "union" or "intersect"
#' @param keep.unique Keep genes that occur (i.e., there is a corresponding column in raw.data) only in one dataset?
#' @return analogizer object
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy = normalize(analogy)
#' analogy = selectGenes(analogy)
#' }

selectGenes = function(object,alphathresh=0.99,varthresh=0.1,cex.use=0.3,combine="union",keep.unique=F)
{
  genes.use = c()
  for (i in 1:length(object@raw.data))
  {
  trx_per_cell <- colSums(object@raw.data[[i]])
  gene_expr_mean <- rowMeans(object@norm.data[[i]])  # Each gene's mean expression level (across all cells)
  gene_expr_var  <- apply( object@norm.data[[i]], 1, var  )  # Each gene's expression variance (across all cells)
  nolan_constant <- mean((1 / trx_per_cell )) 
  alphathresh.corrected=alphathresh/dim(object@raw.data[[i]])[1]
  genemeanupper <- gene_expr_mean+qnorm(1-alphathresh.corrected/2)*sqrt(gene_expr_mean*nolan_constant/dim(object@raw.data[[i]])[2])
  genes.new=names(gene_expr_var)[which(gene_expr_var/nolan_constant> genemeanupper & log10(gene_expr_var) > log10(gene_expr_mean)+(log10(nolan_constant)+varthresh))]
  plot( log10(gene_expr_mean), log10(gene_expr_var), cex=cex.use)
  
  points(log10(gene_expr_mean[genes.new]),log10(gene_expr_var[genes.new]),cex=cex.use,col='green')
  abline(log10(nolan_constant),1,col='purple')
  
  legend("bottomright",paste0("Selected genes: ",length(genes.new)),pch=20,col='green')
  if (combine=="union")
  {
    genes.use = union(genes.use,genes.new)
    if (!keep.unique)
    {
      genes.use = genes.use[genes.use %in% rownames(object@raw.data)]
    }
  }
  if (combine=="intersection")
  {
    genes.use = intersect(genes.use,genes.new)
  }
  }
  object@var.genes = genes.use
  return(object)
}

#' Normalize datasets to column sums
#'
#' @param object analogizer object. 
#' @return analogizer object
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy = normalize(analogy)
#' }

normalize = function(object)
{
  object@norm.data = lapply(object@raw.data,function(x){sweep(x,2,colSums(x),"/")})
  return(object)
}

#' Scale genes to unit variance (but no mean centering).
#'
#' @param object analogizer object. Should normalize and select genes before calling.
#' @return analogizer object
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }

scaleNotCenter = function(object)
{
  object@scale.data = lapply(object@norm.data[object@var.genes,],function(x){scale(t(x),center=F,scale=T)})
  return(object)
}

#' show method for analogizer
#'
#' @param object analogizer object
#' @name show
#' @aliases show,analogizer-method
#' @docType methods
#' @rdname show-methods
#'
setMethod(
  f = "show",
  signature = "analogizer",
  definition = function(object) {
    cat(
      "An object of class",
      class(object),
      "\n",
      length(norm.data),
      "datasets.\n"
    )
    invisible(x = NULL)
  }
)

#' Optimize objective function using block coordinate descent (alternating nonnegative least squares)
#'
#' @param object analogizer object. Should normalize, select genes, and scale before calling.
#' @param k Inner dimension of factorization (number of factors)
#' @param lambda Regularization parameter. Larger values penalize dataset-specific effects more strongly.
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max_iters Maximum number of block coordinate descent iterations to perform
#' @param nrep Number of restarts to perform (NMF objectives are non-convex, so taking the best objective from multiple successive initializations is recommended)
#' @return list of length N+2 containing Hs, W, Vs, and run statistics
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(Y,Z)
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' analogy = optimize_als(analogy,k=2,nrep=1)
#' }
optimizeALS = function(object,k,lambda=5.0,thresh=1e-4,max_iters=100,nrep=20)
{
  E = lapply(object@scale.data,t)
  N = length(E)
  ns = sapply(E,nrow)

  g = ncol(E[[1]])
  W_m = matrix(0, k, g)
  V_m = lapply(1:N,function(i){matrix(0, k, g)})
  H_m = lapply(ns,function(n){matrix(0, n, k)})

  best_obj = Inf
  run_stats = matrix(0,nrow=nrep,ncol=2)
  for (i in 1:nrep)
  {
    start_time <- Sys.time()
    W = matrix(0, k, g)
    V = lapply(1:N,function(i){matrix(abs(runif(g * k,0,2)), k, g)})
    H = lapply(ns,function(n){matrix(abs(runif(n * k,0,2)), n, k)})

    delta = 1
    iters = 0
    pb = txtProgressBar(min=0,max=max_iters,style=3)
    sqrt_lambda = sqrt(lambda)
    obj0 = sum(sapply(1:N,function(i){norm(E[[i]]-H[[i]]%*%(W+V[[i]]),"F")^2}))+sum(sapply(1:N,function(i){lambda*norm(H[[i]]%*%V[[i]],"F")^2}))
    start_obj = obj0
    while(delta > thresh & iters < max_iters)
    {
      H = lapply(1:N,function(i){t(solve_nnls(rbind(t(W)+t(V[[i]]),sqrt_lambda*t(V[[i]])),rbind(t(E[[i]]),matrix(0,nrow=g,ncol=ns[i]))))})
      V = lapply(1:N,function(i){solve_nnls(rbind(H[[i]],sqrt_lambda*H[[i]]),rbind(E[[i]]-H[[i]]%*%W,matrix(0,nrow=ns[[i]],ncol=g)))})
      W = solve_nnls(rbind.fill.matrix(H),rbind.fill.matrix(lapply(1:N,function(i){E[[i]]-H[[i]]%*%V[[i]]})))
      obj = sum(sapply(1:N,function(i){norm(E[[i]]-H[[i]]%*%(W+V[[i]]),"F")^2}))+sum(sapply(1:N,function(i){lambda*norm(H[[i]]%*%V[[i]],"F")^2}))
      delta = abs(obj0-obj)/(mean(obj0,obj))
      obj0 = obj
      iters = iters + 1
      setTxtProgressBar(pb,iters)
    }
    setTxtProgressBar(pb,max_iters)
    if (obj<best_obj)
    {
      W_m = W
      H_m = H
      V_m = V
      best_obj = obj
    }
    end_time <- Sys.time()
    run_stats[i,1]=as.double(difftime(end_time,start_time,units="secs"))
    run_stats[i,2]=iters
    cat("\nConverged in",run_stats[i,1],"seconds,",iters,"iterations. Objective:",obj,"\n")
  }
  cat("\nBest objective:",best_obj,"\n")
  names(H_m)=names(V_m)=names(object@raw.data)
  object@H = H_m
  object@W = W_m
  object@V = V_m
  return(object)
}
