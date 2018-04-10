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
#' @slot clusters Joint cluster assignments for cells
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
    H.norm = "matrix",
    W = "matrix",
    V = "list",
    tsne.coords = "matrix",
    clusters= "factor",
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
#' @param capitalize Capitalize gene names to match homologous genes?
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

selectGenes = function(object,alphathresh=0.99,varthresh=0.1,cex.use=0.3,combine="union",keep.unique=F,capitalize=F)
{
  genes.use = c()
  for (i in 1:length(object@raw.data))
  {
  if (capitalize){rownames(object@raw.data)=toupper(rownames(object@raw.data))}
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
  }
  if (combine=="intersection")
  {
    genes.use = intersect(genes.use,genes.new)
  }
  }
  if (!keep.unique)
  {
    for (i in 1:length(object@raw.data))
    {
      genes.use = genes.use[genes.use %in% rownames(object@raw.data[[i]])] 
    }
  }
  object@var.genes = genes.use
  return(object)
}

#' Normalize raw datasets to column sums
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
  object@scale.data = lapply(object@norm.data,function(x){scale(t(x[object@var.genes,]),center=F,scale=T)})
  for (i in 1:length(object@scale.data))
  {
    object@scale.data[[i]][is.na(object@scale.data[[i]])]=0
  }
  return(object)
}

#' Quantile normalizes cell factors
#'
#' @param object analogizer object. Should run optimizeALS before calling.
#' @param quantiles Number of quantiles to use for quantile normalization
#' @param ref_dataset Name of dataset to use as a "reference" for normalization. By default,
#' the dataset with the largest number of cells is used.
#' 
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
quantile_norm = function(object,quantiles=50,ref_dataset=NULL)
{
  if (is.null(ref_dataset))
  {
    ns = sapply(object@scale.data,nrow)
    ref_dataset = names(object@scale.data)[which.max(ns)]
  }
  
  Hs_scaled = object@H
  for (i in 1:ncol(Hs_scaled[[1]]))
  {
    for (j in 1:length(Hs_scaled)){
      Hs_scaled[[j]][,i] = object@H[[j]][,i]/sqrt(sum(object@H[[j]][,i]^2))
    }
  }
  labels = list()
  for(i in 1:length(Hs_scaled)){
    knn_k=15
    knn = get.knn(object@H[[i]],knn_k)
    labels[[i]] = as.factor(apply(Hs_scaled[[i]],1,which.max))
    labels[[i]] = as.factor(t(apply(knn$nn.index,1,function(x){which.max(table(labels[[i]][x]))}))[1,])
  }
  
  object@clusters = as.factor(unlist(lapply(labels,as.character)))
  
  clusters = labels
  names(clusters)=names(object@H)
  dims = ncol(object@H[[ref_dataset]])
  
  Hs = object@H
  num_clusters = dims
    for (k in 1:length(Hs))
    {
      for (i in 1:dims)
      {
        for (j in 1:num_clusters)
        {
          if (sum(clusters[[ref_dataset]]==j)==0 | sum(clusters[[k]]==j)==0){next}
          if (sum(clusters[[k]]==j)==1){
            Hs[[k]][clusters[[k]]==j,i] = mean(Hs[[ref_dataset]][clusters[[ref_dataset]]==j,i])
            next
          }
          q2 = quantile(Hs[[k]][clusters[[k]]==j,i],seq(0,1,by=1/quantiles))
          q1 = quantile(Hs[[ref_dataset]][clusters[[ref_dataset]]==j,i],seq(0,1,by=1/quantiles))
          if (sum(q1)==0 | sum(q2)==0 | length(unique(q1)) < 2 | length(unique(q2)) < 2)
          {
            new_vals = rep(0,sum(clusters[[k]]==j))
          }
          else
          {
            warp_func = approxfun(q2,q1)  
            new_vals = warp_func(Hs[[k]][clusters[[k]]==j,i])
          }
          
          Hs[[k]][clusters[[k]]==j,i] = new_vals
        }
      }
    }
    object@H.norm = Reduce(rbind,Hs)
  return(object)
}

#' Perform t-SNE on the normalized cell factors to generate a 2D embedding for visualization
#'
#' @param object analogizer object. Should run quantile_norm before calling.
#' 
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
run_tSNE = function(object)
{
  object@tsne.coords = Rtsne(object@H.norm,pca=F)$Y
  return (object)
}

#' Plot t-SNE coordinates of aligned datasets, colored by dataset and by cluster
#'
#' @param object analogizer object. Should call run_tSNE before calling.
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }
plotByDatasetAndCluster = function(object)
{
  tsne_df = data.frame(object@tsne.coords)
  colnames(tsne_df)=c("tsne1","tsne2")
  tsne_df$Dataset = unlist(sapply(1:length(object@H),function(x){rep(names(object@H)[x],nrow(object@H[[x]]))}))
  tsne_df$Cluster = object@clusters
  print(plot_grid(ggplot(tsne_df,aes(x=tsne1,y=tsne2,color=Dataset))+geom_point(),ggplot(tsne_df,aes(x=tsne1,y=tsne2,color=Cluster))+geom_point()))
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
      length(object@norm.data),
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
  E = object@scale.data
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
