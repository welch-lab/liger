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
#' @slot snf Shared nearest factor matrix for use in clustering and alignment
#' @slot agg.data Data aggregated within clusters
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
    agg.data = "list",
    parameters = "list",
    snf = 'list'
  )
)

Analogizer <- function(raw.data, sparse_dcg=T) {
  object <- methods::new(
    Class = "analogizer",
    raw.data = raw.data
  )
  if (sparse_dcg) {
    raw.data = lapply(raw.data, function(x) {
      if (class(x)[1] == "dgTMatrix") {
        temp = summary(x)
        sparseMatrix(i=temp[,1],j=temp[,2],x=temp[,3])
      } else {
        Matrix(as.matrix(x), sparse = T)
      }
    })
    object@raw.data = raw.data
  } 
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
  if (capitalize){
    rownames(object@raw.data[[i]])=toupper(rownames(object@raw.data[[i]]))
    rownames(object@norm.data[[i]])=toupper(rownames(object@norm.data[[i]]))
  }
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
  if (class(object@raw.data[[1]])[1] == 'dgTMatrix' | class(object@raw.data[[1]])[1] == 'dgCMatrix') {
    object@norm.data = lapply(object@raw.data, Matrix.column_norm)
  } else {
    object@norm.data = lapply(object@raw.data,function(x){sweep(x,2,colSums(x),"/")})
  }
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

scaleNotCenter = function(object,cells=NULL)
{
  if(is.null(cells))
  {
    cells = lapply(1:length(object@raw.data),function(i){1:ncol(object@raw.data[[i]])})
  }
  object@scale.data = lapply(1:length(object@norm.data),function(i){
    scale(t(object@norm.data[[i]][object@var.genes,]),center=F,scale=T)
  })
  for (i in 1:length(object@scale.data))
  {
    object@scale.data[[i]][is.na(object@scale.data[[i]])]=0
  }
  object = removeMissingCells(object)
  return(object)
}

# helper function to remove cells not expressing any selected genes
# must have set scale.data first
removeMissingCells = function(object) {
  object@scale.data = lapply(seq_along(object@scale.data), function(x) {
    missing = which(rowSums(object@scale.data[[x]]) == 0)
    if (length(missing) > 0) {
      print(paste0('Removing cells not expressing selected genes in ', 
                   names(object@scale.data)[[x]],':'))
      print(rownames(object@scale.data[[x]])[missing])
      object@scale.data[[x]][-missing,]
    } else {
      object@scale.data[[x]]
    }
  })
  names(object@scale.data) = names(object@raw.data)
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
#' @importFrom FNN get.knn
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }
quantile_norm = function(object,quantiles=50,ref_dataset=NULL,min_cells=2)
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
      #Hs_scaled[[j]][,i] = object@H[[j]][,i]/sqrt(sum(object@H[[j]][,i]^2))
      Hs_scaled[[j]] = scale(Hs_scaled[[j]],scale=T,center=T)
    }
  }
  labels = list()
  for(i in 1:length(Hs_scaled)){
    #knn_k=15
    #knn = get.knn(object@H[[i]],knn_k)
    pct1 = apply(object@H[[i]],2,sum)/sum(apply(object@H[[i]],2,sum))
    pct2 = apply(object@H[[ref_dataset]],2,sum)/sum(apply(object@H[[ref_dataset]],2,sum))
    if (names(object@H)[i]==ref_dataset)
    {
      pct1 = apply(object@H[[i]],2,sum)/sum(apply(object@H[[i]],2,sum))
      pct2 = apply(object@H[[2]],2,sum)/sum(apply(object@H[[2]],2,sum))
    }
    use_these_factors = 1:ncol(object@H[[i]])#which(log(pct1/pct2) > -2)
    
    labels[[i]] = as.factor(use_these_factors[apply(Hs_scaled[[i]][,use_these_factors],1,which.max)])
    #labels[[i]] = as.factor(t(apply(knn$nn.index,1,function(x){which.max(table(labels[[i]][x]))}))[1,])
  }
  
  object@clusters = as.factor(unlist(lapply(labels,as.character)))
  names(object@clusters)=unlist(lapply(object@scale.data,rownames))
  
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
          if (sum(clusters[[ref_dataset]]==j) < min_cells | sum(clusters[[k]]==j) < min_cells){next}
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
#' @param rand.seed Random seed to make results reproducible
#' @param use.raw Use scaled data or factorization result?
#' @param dims.use Which factors to use for computing tSNE embedding
#' @return analogizer object
#' @importFrom Rtsne Rtsne
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }
run_tSNE<-function (object, rand.seed = 42,use.raw = F,dims.use = 1:ncol(object@H.norm),perplexity=30)
{
  set.seed(rand.seed)
  if (use.raw) {
    raw.data = do.call(rbind,object@H)
    object@tsne.coords = Rtsne(raw.data[,dims.use])$Y
  } else {
    object@tsne.coords = Rtsne(object@H.norm[,dims.use], pca = F,check_duplicates = F,perplexity=perplexity)$Y
  }
  rownames(object@tsne.coords) = rownames(object@H.norm)
  return(object)
}

#' Run UMAP on the normalized cell factors to generate a 2D embedding for visualization
#'
#' @param object analogizer object. Should run quantile_norm before calling.
#' @param rand.seed Random seed to make results reproducible
#' @param use.raw Use scaled data or factorization result?
#' @param dims.use Indices of factors to use 
#' @param k Number of dimensions to reduce to
#' @param distance Name of distance metric to use in defining fuzzy simplicial sets
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
run_umap<-function (object, rand.seed = 42, use.raw = F, dims.use = 1:ncol(object@H.norm),
                    k=2, distance = 'euclidean', n_neighbors = 10, min_dist = 0.1)
{
  if (!require("reticulate", quietly = TRUE)) {
    stop("Package \"reticulate\" needed for this function to work. Please install it.\n
         Also ensure Python package umap (PyPI name umap-learn) is installed in python \n
         version accesible to reticulate.",
         call. = FALSE)
  }
  UMAP<-import("umap")
  umapper = UMAP$UMAP(n_components=as.integer(k),metric = distance, n_neighbors = as.integer(n_neighbors),
                      min_dist = min_dist)
  Rumap = umapper$fit_transform
  if (use.raw) {
    raw.data = do.call(rbind,object@H)
    object@tsne.coords = Rumap(raw.data)
  } else {
    object@tsne.coords = Rumap(object@H.norm)
  }
  rownames(object@tsne.coords) = rownames(object@H.norm)
  return(object)
}

#' Plot t-SNE coordinates of aligned datasets, colored by dataset and by cluster
#'
#' @param object analogizer object. Should call run_tSNE before calling.
#' @param title Plot title
#' @param pt.size Controls size of points representing cells
#' @param text.size Controls size of plot text
#' @param do.shuffle Randomly shuffle points so that points from same dataset are not plotted one after the other.
#' @param axis.labels Vector of two strings to use as x and y labels respectively
#' @param return.plots Return ggplot objects instead of printing directly
#' @export
#' @importFrom cowplot rid
#' @importFrom ggplot2 ggplot geom_point aes
#' @importFrom dplyr %>%
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }
plotByDatasetAndCluster<-function(object,title=NULL,pt.size = 0.3,text.size = 3,do.shuffle = T,clusters=NULL,
                                  axis.labels = NULL, return.plots=F, do.legend=T,legend.size=5){
  tsne_df = data.frame(object@tsne.coords)
  colnames(tsne_df) = c("tsne1", "tsne2")
  tsne_df$Dataset = unlist(lapply(1:length(object@H), function(x) {
    rep(names(object@H)[x], nrow(object@H[[x]]))
  }))
  if (is.null(clusters))
  {
    clusters = object@clusters
  }
  tsne_df$Cluster = clusters[names(object@clusters)]
  if (do.shuffle) {
    idx = sample(1:nrow(tsne_df))
    tsne_df = tsne_df[idx,]
  }
  
  p1 = ggplot(tsne_df, aes(x = tsne1, y = tsne2,
                           color = Dataset)) + geom_point(size=pt.size)+ 
    guides(color = guide_legend(override.aes = list(size=legend.size)))
  
  centers <- tsne_df %>% dplyr::group_by(Cluster) %>% dplyr::summarize(tsne1 = median(x = tsne1),
                                                                tsne2 = median(x = tsne2))
  p2 = ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = Cluster)) + geom_point(size=pt.size) + 
          geom_text(data=centers,mapping = aes(label = Cluster),colour='black',size=text.size) + 
      guides(color = guide_legend(override.aes = list(size=legend.size)))
          
  if (!is.null(title)) {
    p1 = p1 + ggtitle(paste0(title,", dataset alignment"))
    p2 = p2 + ggtitle(paste0(title,", published clustering"))
  }
  if (!is.null(axis.labels)) {
    p1 = p1 + xlab(axis.labels[1]) + ylab(axis.labels[2])
    p2 = p2 + xlab(axis.labels[1]) + ylab(axis.labels[2])
  }
  if (!do.legend) {
    p1 = p1 + theme(legend.position = 'none')
    p2 = p2 + theme(legend.position = 'none')
  }
  if (return.plots) {
    return(list(p1, p2))
  } else {
    print(p1)
    print(p2)
  }
  
}

#' Plot comparison scatter plots of unaligned and aligned factor loadings 
#'
#' @param object analogizer object. Should call quantile_align_SNF before calling.
#' @param num_genes Number of genes to display for each factor
#' @param cells.highlight Names of specific cells to highlight in plot (black)
#' @param plot.tsne Plot t-SNE coordinates for each factor 
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }
factor_plots = function (object, num_genes = 10,cells.highlight = NULL, plot.tsne = F)
{
  k = ncol(object@H.norm)
  pb = txtProgressBar(min = 0, max = k, style = 3)
  
  W = t(object@W)
  rownames(W)= colnames(object@scale.data[[1]])
  Hs_norm = object@H.norm
  for (i in 1:k) {
    par(mfrow=c(2,1))
    top_genes.W = rownames(W)[order(W[,i],decreasing=T)[1:num_genes]]
    top_genes.W.string = paste0(top_genes.W,collapse=", ")
    factor_textstring = paste0("Factor",i)

    plot_title1 = paste(factor_textstring,'\n',top_genes.W.string,'\n')
    cols = rep("gray",times=nrow(Hs_norm))
    names(cols) = rownames(Hs_norm)
    cols.use = rainbow(length(object@H))
    
    for (cl in 1:length(object@H)) {
      cols[rownames(object@H[[cl]])] = rep(cols.use[cl],times=nrow(object@H[[cl]]))
    }
    if(!is.null(cells.highlight)) {
      cols[cells.highlight] = rep('black',times = length(cells.highlight))
      
    }
    plot(1:nrow(Hs_norm),do.call(rbind,object@H)[,i],cex=0.2,pch=20,
         col=cols,main=plot_title1,xlab="Cell",ylab="Raw H Score")
    legend("top",names(object@H),pch=20,col=cols.use,horiz=T,cex=0.75)
    plot(1:nrow(Hs_norm),object@H.norm[,i],pch=20,cex=0.2,
         col=cols,xlab="Cell",ylab = "Quantile_norm Score")
    if (plot.tsne) {
      par(mfrow = c(1,1))
      fplot(object@tsne.coords,object@H.norm[,i],title=paste0('Factor ',i))
    }
    setTxtProgressBar(pb, i)
  }
}

# Helper function for factor_plot
fplot = function(tsne,NMFfactor,title,cols.use=heat.colors(10),pt.size=0.7,pch.use=20) {
  data.cut=as.numeric(as.factor(cut(as.numeric(NMFfactor),breaks=length(cols.use))))
  data.col=rev(cols.use)[data.cut]
  plot(tsne[,1],tsne[,2],col=data.col,cex=pt.size,pch=pch.use,main=title)
  
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

rbindlist = function(mat_list)
{
  return(do.call(rbind,mat_list))
}

#' Optimize objective function using block coordinate descent (alternating nonnegative least squares)
#'
#' @param object analogizer object. Should normalize, select genes, and scale before calling.
#' @param k Inner dimension of factorization (number of factors)
#' @param lambda Regularization parameter. Larger values penalize dataset-specific effects more strongly.
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max_iters Maximum number of block coordinate descent iterations to perform
#' @param nrep Number of restarts to perform (NMF objectives are non-convex, so taking the best objective from multiple successive initializations is recommended)
#' @param rand.seed Random seed to allow reproducible results
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
optimizeALS = function(object,k,lambda=5.0,thresh=1e-4,max_iters=100,nrep=1,
                       H_init=NULL,W_init=NULL,V_init=NULL,rand.seed=1)
{
  # remove cells with no selected gene expression 
  object = removeMissingCells(object)
  E = object@scale.data
  N = length(E)
  ns = sapply(E,nrow)
  tmp = gc()
  g = ncol(E[[1]])
  set.seed(rand.seed)
  W_m = matrix(0, k, g)
  V_m = lapply(1:N,function(i){matrix(0, k, g)})
  H_m = lapply(ns,function(n){matrix(0, n, k)})
  tmp = gc()
  
  best_obj = Inf
  run_stats = matrix(0,nrow=nrep,ncol=2)
  for (i in 1:nrep)
  {
    start_time <- Sys.time()
    
    W = matrix(abs(runif(g * k,0,2)), k, g)  
    V = lapply(1:N,function(i){matrix(abs(runif(g * k,0,2)), k, g)})
    H = lapply(ns,function(n){matrix(abs(runif(n * k,0,2)), n, k)})  
    tmp = gc()
    
    if (!is.null(W_init))
    {
      W = W_init  
    }
    if (!is.null(V_init))
    {
      V = V_init  
    }
    if (!is.null(H_init))
    {
      H = H_init  
    }
    
    delta = 1
    iters = 0
    pb = txtProgressBar(min=0,max=max_iters,style=3)
    sqrt_lambda = sqrt(lambda)
    obj0 = sum(sapply(1:N,function(i){norm(E[[i]]-H[[i]]%*%(W+V[[i]]),"F")^2}))+sum(sapply(1:N,function(i){lambda*norm(H[[i]]%*%V[[i]],"F")^2}))
    start_obj = obj0
    tmp = gc()
    
    while(delta > thresh & iters < max_iters)
    {
      H = lapply(1:N,function(i){t(solve_nnls(rbind(t(W)+t(V[[i]]),sqrt_lambda*t(V[[i]])),rbind(t(E[[i]]),matrix(0,nrow=g,ncol=ns[i]))))})
      tmp = gc()
      
      V = lapply(1:N,function(i){solve_nnls(rbind(H[[i]],sqrt_lambda*H[[i]]),rbind(E[[i]]-H[[i]]%*%W,matrix(0,nrow=ns[[i]],ncol=g)))})
      tmp = gc()
      W = solve_nnls(rbindlist(H),rbindlist(lapply(1:N,function(i){E[[i]]-H[[i]]%*%V[[i]]})))
      tmp = gc()
      obj = sum(sapply(1:N,function(i){norm(E[[i]]-H[[i]]%*%(W+V[[i]]),"F")^2}))+sum(sapply(1:N,function(i){lambda*norm(H[[i]]%*%V[[i]],"F")^2}))
      tmp = gc()
      delta = abs(obj0-obj)/(mean(obj0,obj))
      obj0 = obj
      iters = iters + 1
      setTxtProgressBar(pb,iters)
    }
    setTxtProgressBar(pb,max_iters)
    if (iters==max_iters)
    {
      print("Warning: failed to converge within the allowed number of iterations. Re-running with a higher max_iter is recommended.")
    }
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
  object@H = H_m
  object@H=lapply(1:length(object@scale.data),function(i){rownames(object@H[[i]])=rownames(object@scale.data[[i]]);return (object@H[[i]])})
  names(object@H)=names(object@raw.data)
  object@W = W_m
  names(V_m)=names(object@raw.data)
  object@V = V_m
  # set parameter values
  object@parameters$lambda = lambda
  return(object)
}

#' Optimize objective function for new value of k. Uses an efficient strategy for updating that takes advantage of the information in the existing factorization.
#'
#' @param object analogizer object. Should normalize, select genes, and scale before calling.
#' @param k_new Inner dimension of factorization (number of factors)
#' @param lambda Regularization parameter. Larger values penalize dataset-specific effects more strongly.
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max_iters Maximum number of block coordinate descent iterations to perform
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
optimizeNewK = function(object,k_new,lambda=5.0,thresh=1e-4,max_iters=100,rand.seed=1)
{
  k = ncol(object@H[[1]])
  if (k_new == k)
  {
    return(object)
  }
  H = object@H
  W = object@W
  V = object@V
  
  if (k_new > k)
  {
    sqrt_lambda = sqrt(lambda)
    g = ncol(W)
    N = length(H)
    ns = sapply(H,nrow)
    W_new = matrix(abs(runif(g * k,0,2)), k_new-k, g)
    V_new = lapply(1:N,function(i){matrix(abs(runif(g * (k_new-k),0,2)), k_new-k, g)})
    H_new = lapply(ns,function(n){matrix(abs(runif(n * (k_new-k),0,2)), n, k_new-k)})
    H_new = lapply(1:N,function(i){t(solve_nnls(rbind(t(W_new)+t(V_new[[i]]),sqrt_lambda*t(V_new[[i]])),rbind(t(object@scale.data[[i]]-H[[i]]%*%(W+V[[i]])),matrix(0,nrow=g,ncol=ns[i]))))})
    V_new = lapply(1:N,function(i){solve_nnls(rbind(H_new[[i]],sqrt_lambda*H_new[[i]]),rbind(object@scale.data[[i]]-H[[i]]%*%(W+V[[i]])-H_new[[i]]%*%W_new,matrix(0,nrow=ns[[i]],ncol=g)))})
    W_new = solve_nnls(rbind.fill.matrix(H_new),rbind.fill.matrix(lapply(1:N,function(i){object@scale.data[[i]]-H[[i]]%*%(W+V[[i]])-H_new[[i]]%*%V_new[[i]]})))
    H = lapply(1:N,function(i){cbind(H[[i]],H_new[[i]])})
    V = lapply(1:N,function(i){rbind(V[[i]],V_new[[i]])})
    W = rbind(W,W_new)
  }
  else
  {
    deltas = rep(0,k)
    for (i in 1:length(object@H))
    {
      deltas = deltas + sapply(1:k,function(x){norm(H[[i]][,k] %*% t(W[k,]+V[[i]][k,]),"F")})
    }
    k.use = order(deltas,decreasing = T)[1:k_new]
    W = W[k.use,]
    H = lapply(H,function(x){x[,k.use]})
    V = lapply(V,function(x){x[k.use,]})
  }
  object = optimizeALS(object,k_new,H_init=H,W_init=W,V_init=V,nrep=1, rand.seed = rand.seed)
  return(object)
}

#' Optimize objective function for new data. Uses an efficient strategy for updating that takes advantage of the information in the existing factorization.
#'
#' @param object analogizer object. Should normalize, select genes, and scale before calling.
#' @param new.data list of raw data matrices (one or more). Each list entry should be named.
#' @param which.datasets list of datasets to append new.data to if add.to.existing is true. Otherwise, the most similar existing datasets for each entry in new.data.
#' @param add.to.existing add the new data to existing datasets or treat as totally new datasets (new Vs)?
#' @param lambda Regularization parameter. Larger values penalize dataset-specific effects more strongly.
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max_iters Maximum number of block coordinate descent iterations to perform
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
optimizeNewData = function(object,new.data,which.datasets,add.to.existing=T,lambda=5.0,thresh=1e-4,max_iters=100)
{
  if (add.to.existing)
  {
    for (i in 1:length(new.data))
    {
      print(dim(object@raw.data[[which.datasets[[i]]]]))
      object@raw.data[[which.datasets[[i]]]] = cbind(object@raw.data[[which.datasets[[i]]]],new.data[[i]])
      print(dim(object@raw.data[[which.datasets[[i]]]]))
    }
    object = normalize(object)
    object = scaleNotCenter(object)
    sqrt_lambda = sqrt(lambda)
    g = ncol(object@W)
    H_new = lapply(1:length(new.data),function(i){t(solve_nnls(rbind(t(object@W)+t(object@V[[which.datasets[[i]]]]),sqrt_lambda*t(object@V[[which.datasets[[i]]]])),rbind(t(object@scale.data[[which.datasets[[i]]]][colnames(new.data[[i]]),]),matrix(0,nrow=g,ncol=ncol(new.data[[i]])))))})    
    for (i in 1:length(new.data))
    {
      object@H[[which.datasets[[i]]]] = rbind(object@H[[which.datasets[[i]]]],H_new[[i]])
    }
  }
  else
  {
    old.names = names(object@raw.data)
    new.names = names(new.data)
    combined.names = c(old.names,new.names)
    for (i in 1:length(which.datasets))
    {
      object@V[[names(new.data)[i]]] = object@V[[which.datasets[[i]]]]
    }
    object@raw.data = c(object@raw.data,new.data)
    names(object@raw.data)=names(object@V)=combined.names
    object = normalize(object)
    object = scaleNotCenter(object)
    ns = lapply(object@raw.data,ncol)
    N = length(ns)
    g = ncol(object@W)
    sqrt_lambda = sqrt(lambda)
    for (i in 1:N)
    {
      print(ns[[i]])
      print(dim(object@raw.data[[i]]))
      print(dim(object@norm.data[[i]]))
      print(dim(object@scale.data[[i]]))
      print(dim(object@V[[i]]))
    }
    H_new = lapply(1:length(new.data),function(i){t(solve_nnls(rbind(t(object@W)+t(object@V[[new.names[i]]]),sqrt_lambda*t(object@V[[new.names[i]]])),rbind(t(object@scale.data[[new.names[i]]]),matrix(0,nrow=g,ncol=ncol(new.data[[i]])))))})
    object@H = c(object@H,H_new)
    names(object@H) = combined.names
  }
  k = ncol(object@H[[1]])
  object = optimizeALS(object,k,lambda,thresh,max_iters,H_init=object@H,W_init=object@W,V_init=object@V)
  return(object)
}

#' Optimize objective function for a subset of the data. Uses an efficient strategy for updating that takes advantage of the information in the existing factorization.
#'
#' @param object analogizer object. Should normalize, select genes, and scale before calling.
#' @param cell.subset list of cell names to retain from each dataset.
#' @param lambda Regularization parameter. Larger values penalize dataset-specific effects more strongly.
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max_iters Maximum number of block coordinate descent iterations to perform
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
optimizeSubset = function(object,cell.subset=NULL,cluster.subset=NULL,lambda=5.0,thresh=1e-4,max_iters=100,datasets.scale=NULL)
{
  if (is.null(cell.subset) & is.null(cluster.subset))
  {
    print("Please specify a cell subset or cluster subset.")
    return(object)
  }
  else if (is.null(cell.subset) & !is.null(cluster.subset))
  {
    cell.subset = lapply(1:length(object@scale.data),function(i){which(object@clusters[rownames(object@scale.data[[i]])] %in% cluster.subset)})
  }
  old_names = names(object@raw.data)
  H = object@H
  H = lapply(1:length(object@H),function(i){object@H[[i]][cell.subset[[i]],]})
  object@raw.data = lapply(1:length(object@raw.data),function(i){object@raw.data[[i]][,cell.subset[[i]]]})
  for (i in 1:length(object@norm.data))
  {
    object@norm.data[[i]] = object@norm.data[[i]][,cell.subset[[i]]]   
    if (names(object@norm.data)[i] %in% datasets.scale)
    {
      object@scale.data[[i]] = scale(t(object@norm.data[[i]][object@var.genes,]),scale=T,center=F)
      object@scale.data[[i]][is.na(object@scale.data[[i]])] = 0
    }
    else
    {
      object@scale.data[[i]] = t(object@norm.data[[i]][object@var.genes,])
    }
    print(dim(object@scale.data[[i]]))
  }
  
  names(object@raw.data)=names(object@norm.data)=names(object@H)=old_names
  k = ncol(H[[1]])
  object = optimizeALS(object,k=k,lambda=lambda,thresh=thresh,max_iters=max_iters,H_init=H,W_init=object@W,V_init=object@V,nrep=1)
  return(object)
}

#' Optimize objective function for new lambda value. Uses an efficient strategy for updating 
#' that takes advantage of the information in the existing factorization.
#'
#' @param object analogizer object. Should call optimizeALS before calling.
#' @param new_lambda Regularization parameter. Larger values penalize dataset-specific effects more strongly.
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max_iters Maximum number of block coordinate descent iterations to perform
#' @param rand.seed Random seed for reproducibility 
#' @return Analogizer object with optimized factorization values
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
optimizeNewLambda = function(object, new_lambda, thresh=1e-4, max_iters=100, rand.seed = 1) {
  k = ncol(object@H[[1]])
  H = object@H
  W = object@W
  if (new_lambda < object@parameters$lambda) {
    print('New lambda less than current lambda; new factorization may not be optimal.
          Re-optimization with optimizeAlS recommended instead.')
  }
  object = optimizeALS(object, k, lambda = new_lambda, thresh = thresh, max_iters = max_iters,
                       H_init=H, W_init=W, rand.seed = rand.seed)
  return(object)
}

#' Plot alignment and agreement for various test values of lambda. Can be used to select
#' appropriate value of lambda for factorization of particular dataset. 
#'
#' @param object analogizer object. Should normalize, select genes, and scale before calling.
#' @param k Number of factors for factorizations
#' @param lambda_test Vector of lambda values to test. If not given, use default set spanning 
#'                    0.25 to 60
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max_iters Maximum number of block coordinate descent iterations to perform
#' @param rand.seed Random seed for reproducibility 
#' @param num.cores Number of cores to use for optimizing factorizations in parallel
#' @param k2 Horizon parameter for SNF quantile alignment
#' @param ref_dataset Reference dataset for SNF quantile alignment
#' @param resolution Resolution for SNF quantile alignment
#' @param agree.method Reference dr.method for calculating distortion_metric
#' @param gen.new Do not use optimizeNewLambda in factorizations. Recommended to use
#'                when looking at only a small range of lambdas (ie. 1:7)
#' @param return_results Return matrix of alignment and agreement values 
#' @return Matrix of results or plot 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach "%dopar%"
#' @importFrom Hmisc minor.tick
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
lambdaSuggestion = function(object, k, lambda_test = NULL, rand.seed = 1, num.cores = 1, 
                            thresh = 1e-4, max_iters = 100, k2 = 500, ref_dataset=NULL, resolution = 1, 
                            agree.method='PCA', gen.new=F, return_results=F) {
  if (is.null(lambda_test)){
    lambda_test = c(seq(0.25, 1, 0.25), seq(2, 10, 1), seq(15, 60, 5))
  }
  registerDoParallel(cores = num.cores)
  print('This may take several minutes depending on number of values tested')
  print(paste('Optimizing initial factorization with lambda =', lambda_test[1]))
  object = optimizeALS(object,k=k,thresh = thresh, lambda = lambda_test[1], max_iters=max_iters,
                       nrep=1, rand.seed = rand.seed)
  data_matrix <- foreach(i=1:length(lambda_test), .combine = 'rbind') %dopar% {
    if (i != 1) {
      if (gen.new) {
        ob.test = optimizeALS(object, k=k, lambda = lambda_test[i], thresh = thresh, 
                              max_iters = max_iters, rand.seed = rand.seed)
      } else {
        ob.test = optimizeNewLambda(object, new_lambda = lambda_test[i], thresh = thresh, 
                                    max_iters = max_iters, rand.seed = rand.seed)
      }
    } else {
      ob.test = object
    }
    ob.test = quantile_align_SNF(ob.test, k2 = k2, resolution = resolution, ref_dataset = ref_dataset,
                                 id.number = i) 
    align = alignment_metric(ob.test)
    agree_unaligned = distortion_metric(ob.test, ndims = k, k = 15, dr_method = agree.method, 
                                        use_aligned = F)
    agree_aligned = distortion_metric(ob.test, ndims = k, k = 15, dr_method = agree.method, 
                                      use_aligned = T)
    c(align, agree_unaligned, agree_aligned)
  }
  # plot results on same plot 
  plot(lambda_test, data_matrix[,1], type='p', col='black', ylim=c(0,1), xlab = 'Lambda', 
       ylab='Value')
  minor.tick(nx=4, ny=2, tick.ratio=0.5)
  grid()
  lines(lambda_test, data_matrix[,1], col='black')
  lines(lambda_test, data_matrix[,2], col='blue')
  lines(lambda_test, data_matrix[,3], col='green')
  legend('bottomright', legend=c('Alignment', paste0('Agreement(',agree.method, ')'), 
                                 paste0('Agreement(',agree.method, '-aligned)')), 
         col = c('black', 'blue', 'green'), lty=1, cex=0.8)
  
  if (return_results) {
    data_matrix = cbind(lambda_test, data_matrix)
    return(data_matrix)
  }
}

#' Plot alignment and agreement for various test values of lambda. Can be used to select
#' appropriate value of lambda for factorization of particular dataset. 
#'
#' @param object analogizer object. Should normalize, select genes, and scale before calling.
#' @param k_test Set of factor numbers to test (default seq(5, 50, 5))
#' @param lambda Lambda for optimization
#' @param thresh Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
#' @param max_iters Maximum number of block coordinate descent iterations to perform
#' @param rand.seed Random seed for reproducibility 
#' @param num.cores Number of cores to use for optimizing factorizations in parallel
#' @param gen.new Do not use optimizeNewK in factorizations. 
#' @param return_results Return matrix of alignment and agreement values 
#' @return Matrix of results or plot 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach "%dopar%"
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
kSuggestion = function(object, k_test=seq(5, 50, 5), lambda=5, thresh=1e-4, max_iters=100, num.cores=1, 
                       rand.seed = 1, plot.metric='median', gen.new=F, plot.log2=T, return_results=F) {
  registerDoParallel(cores = num.cores)
  
  # optimize largest k value first to take advantage of efficient updating 
  print('This operation may take several minutes depending on number of values being tested')
  print(paste0('Optimizing initial factorization with largest test k=', k_test[length(k_test)]))
  object = optimizeALS(object,k=k_test[length(k_test)],thresh = thresh, lambda, max_iters=max_iters,
                       nrep=1, rand.seed = rand.seed)
  data_matrix <- foreach(i=length(k_test):1, .combine = 'rbind') %dopar% {
    if (i != length(k_test)) {
      if (gen.new) {
        ob.test = optimizeALS(object, k=k_test[i], lambda = lambda, thresh = thresh, 
                              max_iters = max_iters, rand.seed = rand.seed)
      } else {
        ob.test = optimizeNewK(object, k_new = k_test[i], lambda=lambda, thresh = thresh, 
                               max_iters = max_iters, rand.seed = rand.seed)
      }
    } else {
      ob.test = object
    }
    dataset_split = kl_divergence_uniform(ob.test)
    unlist(dataset_split)
  }
  data_matrix = data_matrix[nrow(data_matrix):1,]
  medians = apply(data_matrix, 1, median)
  
  # plot out results 
  max_lim = max(log2(k_test)) + 0.05
  min_lim = min(medians) - 0.05
  if (plot.log2) {
    plot(k_test, log2(k_test), type='p', col='green', ylim=c(min_lim, max_lim))
    points(k_test, medians, type='p', xlab='Number of factors', 
           ylab='Median KL divergence across combined data', col='black')
    legend('topleft', legend=(c('log2(k) (upper lim)', 'KL div')),
           col=c('green', 'black'), lty=1, cex=0.8)
  } else {
    plot(k_test, medians, type='p', xlab='Number of factors', 
         ylab='Median KL divergence across combined data', col='black')
  }
  lines(k_test, medians, col='black')
  
  if (return_results) {
    data_matrix = cbind(k_test, data_matrix)
    return(data_matrix)
  }
}

# helper function for calculating KL divergence from uniform distribution
# (related to Shannon entropy) for factorization
kl_divergence_uniform = function(object) 
{
  n_cells = sum(sapply(object@H, nrow))
  n_factors = ncol(object@H[[1]])
  dataset_list = list()
  for (i in 1:length(object@H)) {
    scaled = scale(object@H[[i]], center=F, scale=T)
    
    inflated = t(apply(scaled, 1, function(x) {
      replace(x, x == 0, 1e-20)
    }))
    inflated = inflated/rowSums(inflated)
    divs = apply(inflated, 1, function(x) {log2(n_factors) + sum(log2(x) * x)})
    dataset_list[[i]] = divs
  }
  return(dataset_list)
}

#' Find shared and dataset-specific markers using the factorization. 
#'
#' @param object analogizer object. Should call run_tSNE before calling.
#' @param dataset1 Name of first dataset
#' @param dataset2 Name of second dataset
#' @param factor.share.thresh Use only factors with a dataset specificity less than threshold
#' @param dataset_specificity Pre-calculated dataset specificity if available
#' @param log_fc_thresh Lower log-fold change threshold for differential expression in markers
#' @param umi_thresh Lower UMI threshold for markers
#' @param frac_thresh Lower threshold for fraction of cells expressing marker
#' @param pval_thresh Upper p-value threshold for Wilcoxon rank test for gene expression 
#' @param num_genes Max number of genes to report for each dataset
#' @param print.genes Print ordered markers passing logfc, umi and frac thresholds
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
get_factor_markers = function(object,dataset1=NULL,dataset2=NULL,factor.share.thresh=10,
                              dataset_specificity=NULL,log_fc_thresh=1,umi_thresh=30,frac_thresh=0,
                              pval_thresh=0.05,num_genes=30,print.genes=F)
{
  if (is.null(dataset1) | is.null(dataset2))
  {
    dataset1 = names(object@H)[1]
    dataset2 = names(object@H)[2]
  }
  if (is.null(num_genes))
  {
    num_genes = length(object@var.genes)
  }
  if (is.null(dataset_specificity)){
    dataset_specificity = calc_dataset_specificity(object)
  }
  factors.use = which(abs(dataset_specificity[[3]]) <= factor.share.thresh)
  
  if (length(factors.use)<2)
  {
    print(paste("Warning: only",length(factors.use),"factors passed the dataset specificity threshold."))
  }
  
  Hs_scaled = lapply(object@H, function(x) {scale(x, scale=T, center=T)})
  labels = list()
  for(i in 1:length(Hs_scaled)){
    labels[[i]] = factors.use[as.factor(apply(Hs_scaled[[i]][,factors.use],1,which.max))]
  }
  names(labels)=names(object@H)
  
  V1_matrices = list()
  V2_matrices = list()
  W_matrices = list()
  for (j in 1:length(factors.use))
  {
    i=factors.use[j]
    
    W = t(object@W)
    V1 = t(object@V[[dataset1]])
    V2 = t(object@V[[dataset2]])
    rownames(W)=rownames(V1)=rownames(V2)=object@var.genes
    # if not max factor for any cell in either dataset
    if(sum(labels[[dataset1]]==i) <= 1 | sum(labels[[dataset2]]==i) <= 1)
    {
      print(paste('Warning: factor', i,'did not appear as max in any cell in either dataset'))
      next
    }
    # filter genes by gene_count and cell_frac thresholds
    gene_info = list()
    for (dset in c(dataset1, dataset2)) {
      gene_info[[dset]]$gene_counts = rowSums(object@raw.data[[dset]][object@var.genes,labels[[dset]]==i])
      gene_info[[dset]]$cell_fracs = rowSums(object@raw.data[[dset]][object@var.genes,
                                                                     labels[[dset]]==i]>0)/sum(labels[[dset]]==i)
      gene_info[[dset]]$norm_counts = object@norm.data[[dset]][object@var.genes,labels[[dset]]==i]
      gene_info[[dset]]$mean = rowMeans(object@norm.data[[dset]][object@var.genes,labels[[dset]]==i])
      names(gene_info[[dset]]$gene_counts) = names(gene_info[[dset]]$mean) = object@var.genes
      rownames(gene_info[[dset]]$norm_counts) = object@var.genes
    }
    log2fc = log2(gene_info[[dataset1]]$mean / gene_info[[dataset2]]$mean)
    initial_filtered = object@var.genes[(gene_info[[dataset1]]$gene_counts > umi_thresh | 
                                           gene_info[[dataset2]]$gene_counts > umi_thresh) & 
                                          (gene_info[[dataset1]]$cell_fracs > frac_thresh | 
                                             gene_info[[dataset2]]$cell_fracs > frac_thresh)]
    filtered_genes_V1 = initial_filtered[log2fc[initial_filtered] > log_fc_thresh]
    filtered_genes_V2 = initial_filtered[(-1 * log2fc)[initial_filtered] > log_fc_thresh]
    
    W = pmin(W+V1,W+V2)
    V1 = V1[filtered_genes_V1,]
    V2 = V2[filtered_genes_V2,]
    
    if(length(filtered_genes_V1)==0)
    {
      top_genes_V1 = ""
    } else {
      top_genes_V1 = row.names(V1)[ order(V1[,i], decreasing=T )[1:num_genes] ]
      top_genes_V1 = top_genes_V1[!is.na(top_genes_V1)]
      top_genes_V1 = top_genes_V1[which(V1[top_genes_V1,i]>0)]
    }
    if(length(filtered_genes_V2)==0)
    {
      top_genes_V2 = ""
    } else {
      top_genes_V2 = row.names( V2 )[ order(V2[,i], decreasing=T )[1:num_genes] ]
      top_genes_V2 = top_genes_V2[!is.na(top_genes_V2)]
      top_genes_V2 = top_genes_V2[which(V2[top_genes_V2,i]>0)]
    }
    top_genes_W = row.names( W )[ order(W[,i], decreasing=T )[1:num_genes] ]
    top_genes_W = top_genes_W[!is.na(top_genes_W)]
    top_genes_W = top_genes_W[which(W[top_genes_W,i]>0)]
    
    if(print.genes)
    {
      print(paste("Factor",i))
      print(top_genes_V1)
      print(top_genes_W)
      print(top_genes_V2)
    }
    
    pvals = list() # order is V1, V2, W
    top_genes = list(top_genes_V1, top_genes_V2, top_genes_W)
    for (k in 1:length(top_genes)) {
      pvals[[k]] = sapply(top_genes[[k]],function(x){
        suppressWarnings(wilcox.test(as.numeric(gene_info[[dataset1]]$norm_counts[x,]),
                                     as.numeric(gene_info[[dataset2]]$norm_counts[x,]))$p.value)
      })
    }
    # bind values in matrices
    V1_matrices[[j]] = Reduce(cbind, list(rep(i, length(top_genes_V1)), top_genes_V1, 
                                          gene_info[[dataset1]]$gene_counts[top_genes_V1],
                                          gene_info[[dataset2]]$gene_counts[top_genes_V1],
                                          gene_info[[dataset1]]$cell_fracs[top_genes_V1],
                                          gene_info[[dataset2]]$cell_fracs[top_genes_V1],
                                          log2fc[top_genes_V1], pvals[[1]]))
    V2_matrices[[j]] = Reduce(cbind, list(rep(i, length(top_genes_V2)), top_genes_V2, 
                                          gene_info[[dataset1]]$gene_counts[top_genes_V2],
                                          gene_info[[dataset2]]$gene_counts[top_genes_V2],
                                          gene_info[[dataset1]]$cell_fracs[top_genes_V2],
                                          gene_info[[dataset2]]$cell_fracs[top_genes_V2],
                                          log2fc[top_genes_V2], pvals[[2]]))
    W_matrices[[j]] = Reduce(cbind, list(rep(i, length(top_genes_W)), top_genes_W, 
                                         gene_info[[dataset1]]$gene_counts[top_genes_W],
                                         gene_info[[dataset2]]$gene_counts[top_genes_W],
                                         gene_info[[dataset1]]$cell_fracs[top_genes_W],
                                         gene_info[[dataset2]]$cell_fracs[top_genes_W],
                                         log2fc[top_genes_W], pvals[[3]]))
  }
  V1_genes = data.frame(Reduce(rbind, V1_matrices), stringsAsFactors = F)
  V2_genes = data.frame(Reduce(rbind, V2_matrices), stringsAsFactors = F)
  W_genes = data.frame(Reduce(rbind, W_matrices), stringsAsFactors = F)
  df_cols = c('factor_num', 'gene', 'counts1', 'counts2', 'fracs1', 'fracs2', 'log2fc', 'p_value')
  output_list = list(V1_genes, W_genes, V2_genes)
  output_list = lapply(output_list, function(df) {
    colnames(df) = df_cols
    df = transform(df, factor_num = as.numeric(factor_num), gene=as.character(gene),
                   counts1 = as.numeric(counts1), counts2 = as.numeric(counts2), 
                   fracs1 = as.numeric(fracs1), fracs2 = as.numeric(fracs2), 
                   log2fc = as.numeric(log2fc), p_value = as.numeric(p_value))
    df[which(df$p_value < pval_thresh),]
  })
  names(output_list) = c(dataset1, 'shared', dataset2)
  output_list[['num_factors_V1']] = table(output_list[[dataset1]]$gene)
  output_list[['num_factors_V2']] = table(output_list[[dataset2]]$gene)
  return(output_list)
}

#' Word cloud plots coloring t-SNE points by their loading on specifed factors as well as the
#' most highly loading shared and dataset-specific genes
#'
#' @param object analogizer object. Should call run_tSNE before calling.
#' @param num_genes Number of genes to show in word clouds
#' @param min_size Size of smallest gene symbol in word cloud
#' @param max_size Size of largest gene symbol in word cloud
#' @param dataset1 Name of first dataset
#' @param dataset2 Name of second dataset
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid roundrectGrob
#' @importFrom grid gpar
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
plot_word_clouds = function(object,num_genes=30,min_size=1,max_size=4,dataset1=NULL,
                            dataset2=NULL,factor.share.thresh=10, log_fc_thresh=1,
                            umi_thresh=30,frac_thresh=0,pval_thresh=0.05)
{
  if (is.null(dataset1)|is.null(dataset2))
  {
    dataset1 = names(object@H)[1]
    dataset2 = names(object@H)[2]
  }
  
  H_aligned = object@H.norm
  W = t(object@W)
  V1 = t(object@V[[dataset1]])
  V2 = t(object@V[[dataset2]])
  W = pmin(W+V1,W+V2)
  
  dataset_specificity = calc_dataset_specificity(object)
  factors.use = which(abs(dataset_specificity[[3]]) <= factor.share.thresh)
  
  markers = get_factor_markers(object,factor.share.thresh=factor.share.thresh,
                               num_genes=num_genes, log_fc_thresh = log_fc_thresh,
                               umi_thresh = umi_thresh, frac_thresh = frac_thresh,
                               pval_thresh = pval_thresh, 
                               dataset_specificity = dataset_specificity)
  
  rownames(W)=rownames(V1)=rownames(V2)=object@var.genes
  loadings_list = list(V1, W, V2)
  names_list = list(dataset1, 'Shared', dataset2)
  tsne_coords = object@tsne.coords
  pb = txtProgressBar(min=0,max=length(factors.use),style=3)
  for (i in factors.use)
  {
    tsne_df = data.frame(H_aligned[,i],tsne_coords)
    factorlab = paste("Factor",i,sep="")
    colnames(tsne_df)=c(factorlab,"tSNE1","tSNE2")
    factor_ds = paste("Factor",i,"Dataset Specificity:",dataset_specificity[[3]][i])
    p1 = ggplot(tsne_df,aes_string(x="tSNE1",y="tSNE2",color=factorlab))+geom_point()+
      scale_color_gradient(low="yellow",high="red")+ggtitle(label=factor_ds)
    
    top_genes_V1 = markers[[1]]$gene[markers[[1]]$factor_num==i]
    top_genes_W = markers[[2]]$gene[markers[[2]]$factor_num==i]
    top_genes_V2 = markers[[3]]$gene[markers[[3]]$factor_num==i]
    
    top_genes_list = list(top_genes_V1, top_genes_W, top_genes_V2)
    plot_list = lapply(seq_along(top_genes_list), function(x) {
      top_genes = top_genes_list[[x]]
      gene_df = data.frame(genes=top_genes,
                           loadings=loadings_list[[x]][top_genes,i])
      if(length(top_genes)==0){gene_df = data.frame(genes=c("no genes"),loadings=c(1))}
      out_plot = ggplot(gene_df,aes(x = 1, y = 1, size = loadings, label = genes)) +
        geom_text_repel(force = 100,segment.color=NA) +
        scale_size(range = c(min_size, max_size), guide = FALSE) +
        scale_y_continuous(breaks = NULL) +
        scale_x_continuous(breaks = NULL) +
        labs(x = '', y = '')+ggtitle(label=names_list[[x]])+coord_fixed()
      return(out_plot)
    })
    
    p2 = (plot_grid(plotlist = plot_list,align="hv",nrow = 1)
          + draw_grob(roundrectGrob(x=0.33,y=0.5,width=0.67,height=0.70,
                                    gp = gpar(fill = "khaki1", col = "Black",alpha=0.5,lwd=2))) 
          + draw_grob(roundrectGrob(x=0.67,y=0.5,width=0.67,height=0.70,
                                    gp = gpar(fill = "indianred1", col = "Black",alpha=0.5,lwd=2))))
    print(plot_grid(p1,p2,nrow=2,align="h"))
    setTxtProgressBar(pb,i)
  }
} 

#' Calculate distortion statistic to quantify how much alignment distorts
#' the geometry of the original datasets.
#'
#' @param object analogizer object. Should call quantile_norm before calling.
#' @param dr_method Dimensionality reduction method to use for assessing pre-alignment geometry (either "PCA", "NMF", or "ICA"). Should call quantile_norm before calling.
#' @return alignment statistic
#' @importFrom FNN get.knn
#' @importFrom NNLM nnmf
#' @importFrom ica icafast
#' @importFrom irlba prcomp_irlba
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
distortion_metric = function(object,dr_method="PCA",ndims=40,k=10, use_aligned=TRUE, 
                             rand.seed=42, return_breakdown=FALSE)
{
  print(paste("Reducing dimensionality using",dr_method))
  dr = list()
  if (dr_method=="NMF")
  {
    set.seed(rand.seed)
    dr = lapply(object@scale.data,function(x){nnmf(x,k=ndims)$W})
  }
  else if(dr_method=="ICA")
  {
    set.seed(rand.seed)
    dr = lapply(object@scale.data,function(x){icafast(x,nc=ndims)$S})
  }
  else #PCA
  {
    dr = lapply(object@scale.data,function(x){suppressWarnings(prcomp_irlba(t(x),n = ndims,
                                                                            scale. = (colSums(x)>0),center=F)$rotation)})
    for (i in 1:length(dr)) 
    {
      rownames(dr[[i]]) = rownames(object@scale.data[[i]])
    }
  }
  ns = sapply(object@scale.data,nrow)
  n = sum(ns)
  jaccard_inds = c()
  distorts = c()
  
  for (i in 1:length(dr))
  {
    jaccard_inds_i = c()
    if (use_aligned)
    {
      original = object@H.norm[rownames(dr[[i]]),]
    } else 
    {
      original = object@H[[i]]
    }
    fnn.1 = get.knn(dr[[i]],k=k)
    fnn.2 = get.knn(original,k=k)
    jaccard_inds_i = c(jaccard_inds_i,sapply(1:ns[i],function(i){
      intersect = intersect(fnn.1$nn.index[i,],fnn.2$nn.index[i,])
      union = union(fnn.1$nn.index[i,],fnn.2$nn.index[i,])
      length(intersect)/length(union)
    }))
    jaccard_inds_i = jaccard_inds_i[is.finite(jaccard_inds_i)]
    jaccard_inds = c(jaccard_inds,jaccard_inds_i)
    
    distorts = c(distorts, mean(jaccard_inds_i))
  }
  if (return_breakdown) 
  {
    return(distorts)
  }
  return(mean(jaccard_inds))
}

#' Calculate alignment statistic to quantify how well-aligned two or more datasets are.
#'
#' @param object analogizer object. Should call quantile_norm before calling.
#' @param k Number of nearest neighbors to use in calculating alignment.
#' @param rand.seed Random seed for reproducibility
#' @param cells.use Vector of cells across all datasets to use in calculating score
#' @param clusters.use Names of clusters to use in calculating score 
#' @return alignment statistic
#' @importFrom FNN get.knn
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
alignment_metric<-function(object, k=NULL, rand.seed=1, cells.use = NULL, 
                           clusters.use = NULL, by_dataset=F)
{
  if (is.null(cells.use)) {
    cells.use = rownames(object@H.norm)
  } 
  if (!is.null(clusters.use)){
    cells.use = names(object@clusters)[which(object@clusters %in% clusters.use)]
  } 
  nmf_factors = object@H.norm[cells.use,]
  num_cells = length(cells.use)
  func_H = lapply(seq_along(object@H), function(x){
    cells.overlap = intersect(cells.use, rownames(object@H[[x]]))
    if (length(cells.overlap) > 0) {
      object@H[[x]][cells.overlap,]
    } else {
      warning(paste0('Selected subset eliminates dataset ', names(object@H)[x]), 
              immediate.=T)
      return(NULL)
    }
  })
  func_H = func_H[!sapply(func_H, is.null)]
  num_factors = ncol(object@H.norm)
  N = length(func_H)
  if (N == 1) {
    warning('Alignment null for single dataset', immediate. = T)
  }
  set.seed(rand.seed)
  min_cells = min(sapply(func_H, function(x){nrow(x)}))
  sampled_cells = unlist(lapply(1:N, function(x) {
    sample(rownames(func_H[[x]]), min_cells)
  }))
  max_k = length(sampled_cells)-1
  if (is.null(k)) {
    k = min(max(floor(0.01 * num_cells),10), max_k)
  } else if (k > max_k) {
    stop(paste0('Please select k <=', max_k))
  }
  knn_graph = get.knn(nmf_factors[sampled_cells, 1:num_factors], k)
  dataset = unlist(sapply(1:N, function(x) {
    rep(names(object@H)[x], nrow(func_H[[x]]))
  }))
  names(dataset)=cells.use
  dataset = dataset[sampled_cells]
  num_sampled = N*min_cells
  num_same_dataset = rep(k, num_sampled)
  
  for (i in 1:num_sampled) {
    inds = knn_graph$nn.index[i, ]
    num_same_dataset[i] = sum(dataset[inds] == dataset[i])
  }
  if (by_dataset) {
    alignments = c()
    for (i in 1:N) {
      start = 1+ (i-1)*min_cells
      end = i*min_cells
      alignment = 1 - ((mean(num_same_dataset[start:end]) - (k/N))/(k-k/N))
      alignments = c(alignments, alignment)
    }
    return(alignments)
  }
  return(1 - ((mean(num_same_dataset) - (k/N))/(k-k/N)))
}

alignment_metric_per_factor<-function(object,k=NULL)
{
  H1_scaled = scale(object@H[[1]],center=F,scale=T)
  H2_scaled = scale(object@H[[2]],center=F,scale=T)
  H_scaled = rbind(H1_scaled,H2_scaled)
  num_cells = nrow(H_scaled)
  num_factors = ncol(H_scaled)
  nmf_factors = H_scaled
  N = length(object@H)
  rownames(nmf_factors)=names(object@clusters)
  num_clusters = length(levels(object@clusters))
  
  dataset = unlist(sapply(1:N, function(x) {
    rep(names(object@H)[x], nrow(object@H[[x]]))
  }))
  names(dataset)=names(object@clusters)
  align_metrics = rep(0,num_clusters)
  for (i in 1:num_clusters)
  {
    cells_i = which(object@clusters==levels(object@clusters)[i])
    cell_names = names(object@clusters)[cells_i]
    num_cells = length(cells_i)
    if(is.null(k))
    {
      k = max(floor(0.05 * num_cells),10)
      print(k)    
    }
    knn_graph = get.knn(nmf_factors[, 1:num_factors], k)
    
    num_same_dataset = rep(0, num_cells)
    num_diff_dataset = rep(0, num_cells)
    for (j in 1:length(cells_i)) {
      inds = knn_graph$nn.index[cells_i[j], ]
      num_same_dataset[j] = sum(dataset[inds] == dataset[cells_i[j]])
      num_diff_dataset[j] = sum(dataset[inds] != dataset[cells_i[j]])
    }
    align_metrics[i] = sum(num_diff_dataset/(num_same_dataset+num_diff_dataset))
  }
  return(align_metrics)
}

#' Calculate alignment statistic for each cluster in aligned datasets. 
#'
#' @param object analogizer object. Should call quantile_norm before calling.
#' @param rand.seed Random seed for reproducibility
#' @param k Number of nearest neighbors in calculating alignment
#'          Can pass in single value or vector with same length as clusters
#' @return Vector of alignment statistics (with names of clusters)
#' @importFrom FNN get.knn
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
alignment_metric_per_cluster<-function(object, rand.seed=1, k=NULL, by_dataset=F)
{
  clusters = levels(object@clusters)
  if (typeof(k) == 'double') {
    if (length(k) == 1) {
      k = rep(k, length(clusters))
    } else if (length(k) != length(clusters)) {
      print('Length of k does not match length of clusters')
    }
  } 
  align_metrics = sapply(seq_along(clusters), function(x) {
    alignment_metric(object, k=k[x], rand.seed = rand.seed, 
                     clusters.use = clusters[x], 
                     by_dataset = by_dataset)
  })
  names(align_metrics) = levels(object@clusters)
  return(align_metrics)
}

#' Perform graph-based clustering (Louvain algorithm) using number of shared nearest neighbors (Jaccard index) as a distance metric.
#'
#' @param object analogizer object. Should call quantile_norm before calling.
#' @param res.param cluster resolution parameter
#' @param k.param nearest neighbor parameter for shared nearest neighbor graph construction
#' @param k.scale scale parameter for shared nearest neighbor graph construction
#' @return analogizer object with cluster assignments
#' @importFrom Seurat FindClusters
#' @importFrom Seurat CreateSeuratObject
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(Y,Z)
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' analogy = optimize_als(analogy,k=2,nrep=1)
#' analogy = quantile_norm(analogy)
#' analogy = clusterLouvainJaccard(object)
#' }
clusterLouvainJaccard = function(object,res.param=0.1,k.param=30)
{
  temp.seurat = CreateSeuratObject(t(Reduce(rbind,object@scale.data)))
  temp.seurat@scale.data = t(Reduce(rbind,object@scale.data))
  rownames(object@H.norm)=colnames(temp.seurat@scale.data)
  temp.seurat@dr$NMF=new(Class="dim.reduction",cell.embeddings=object@H.norm,key="NMF")
  temp.seurat <- FindClusters(object = temp.seurat, reduction.type = "NMF", dims.use = 1:ncol(object@H.norm),force.recalc=T,save.SNN = T,resolution=res.param,k.param=k.param)
  object@clusters = temp.seurat@ident
  return(object)
}

#' Aggregate gene-level measurements across cells within clusters to allow correlation across datasets
#'
#' @param object analogizer object. Should run quantile_norm and possibly clusterLouvainJaccard before calling.
#' @return analogizer object with agg.data
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(Y,Z)
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' analogy = optimize_als(analogy,k=2,nrep=1)
#' analogy = quantile_norm(analogy)
#' analogy = clusterLouvainJaccard(object)
#' }
aggregateByCluster = function(object)
{
  object@agg.data = list()
  for (i in 1:length(object@raw.data))
  {
    clusters_i = object@clusters[names(object@clusters)%in%colnames(object@raw.data[[i]])]
    temp = matrix(0,nrow(object@raw.data[[i]]),length(levels(object@clusters)))
    for (j in 1:length(levels(object@clusters)))
    {
      temp[,j]=rowSums(object@raw.data[[i]][,clusters_i==levels(object@clusters)[j]])
    }
    object@agg.data[[names(object@raw.data)[i]]] = temp
    rownames(object@agg.data[[i]])=rownames(object@raw.data[[i]])
    print(dim(object@agg.data[[i]]))
  }
  object@agg.data = lapply(object@agg.data,function(x){sweep(x,2,colSums(x),"/")})
  
  return(object)
}

plot_violin_summary = function(object,cluster,genes.use)
{
  pdf(filename)
  nmf_factors = data.frame(rbind(object@H[[1]],object@H[[2]]))
  colnames(nmf_factors)=c("tsne1","tsne2")
  nmf_factors$Protocol = as.factor(unlist(sapply(1:length(object@H),function(x){rep(names(object@H)[x],nrow(object@H[[x]]))})))
  if(is.null(clusters))
  {
    clusters = object@clusters
  }
  nmf_factors$Cluster = clusters[names(object@clusters)]
  for (i in 1:(ncol(nmf_factors)-2))
  {
    print(ggplot(nmf_factors,aes(x=Protocol,y=nmf_factors[,i],fill=Protocol))+geom_violin() + geom_jitter(aes(colour=Cluster),shape=16,position=position_jitter(0.4),size=0.6) + guides(colour = guide_legend(override.aes = list(size=4)))+labs(y=paste("Factor",i)))
    print(ggplot(nmf_factors,aes(x=Cluster,y=nmf_factors[,i],fill=Cluster))+geom_violin() + geom_jitter(aes(colour=Protocol),shape=16,position=position_jitter(0.4),size=0.6) + guides(colour = guide_legend(override.aes = list(size=4)))+labs(y=paste("Factor",i)))
  }
  dev.off()
}

#' Uses the factorization information to calculate a dataset-specificity score for each factor
#'
#' @param object analogizer object. Should run optimizeALS before calling.
#' @return List containing three elements. First two elements are the norm of each metagene factor for each dataset. 
#' Last element is the vector of dataset specificity scores.
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(Y,Z)
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' analogy = optimize_als(analogy,k=2,nrep=1)
#' analogy = quantile_norm(analogy)
#' analogy = clusterLouvainJaccard(object)
#' }
calc_dataset_specificity=function(object)
{
  k = ncol(object@H[[1]])
  pct1 = rep(0,k)
  pct2 = rep(0,k)
  for (i in 1:k)
  {
    pct1[i] = norm(as.matrix(object@V[[1]][i,]+object@W[i,]),"F") #norm(object@H[[1]][,i] %*% t(object@W[i,] + object@V[[1]][i,]),"F")
    pct2[i] = norm(as.matrix(object@V[[2]][i,]+object@W[i,]),"F") #norm(object@H[[2]][,i] %*% t(object@W[i,] + object@V[[2]][i,]),"F")
  }
  #pct1 = pct1/sum(pct1)
  #pct2 = pct2/sum(pct2)
  barplot(100*(1-(pct1/pct2)),xlab="Factor",
          ylab="Percent Specificity",main="Dataset Specificity of Factors",
          names.arg = 1:k, cex.names = 0.75, mgp = c(2, 0.5, 0)) # or possibly abs(pct1-pct2)
  return(list(pct1,pct2,100*(1-(pct1/pct2))))
}

#' Plot violin plots for expression of specified gene for each dataset
#'
#' @param object analogizer object. 
#' @param gene Gene for which to plot relative expression. 
#' @param methylation_indices Indices of datasets in object with methylation data. 
#' @param return.plots Return ggplot objects instead of printing directly
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggplot geom_point aes_string scale_color_gradient2 ggtitle
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }

plot_gene_violin = function(object, gene, methylation_indices=NULL, 
                            by.dataset=T, return.plots=F)
{
  gene_vals = c()
  for (i in 1:length(object@norm.data))
  {
    if (i %in% methylation_indices) {
      gene_vals = c(gene_vals,object@norm.data[[i]][gene,])
    } else {
      if (gene %in% rownames(object@norm.data[[i]]))
      {
        gene_vals_int = log2(10000*object@norm.data[[i]][gene,] + 1)
      }
      else
      {
        gene_vals_int = rep(list(0), ncol(object@norm.data[[i]]))
        names(gene_vals_int) = colnames(object@norm.data[[i]])
      }
      gene_vals = c(gene_vals, gene_vals_int)
    }
  }
  
  gene_df = data.frame(object@tsne.coords)
  rownames(gene_df)=names(object@clusters)
  gene_df$Gene = as.numeric(gene_vals[rownames(gene_df)])
  colnames(gene_df)=c("tSNE1","tSNE2","gene")
  gene_plots = list()
  for (i in 1:length(object@scale.data))
  {
    if (by.dataset) {
      gene_df.sub = gene_df[rownames(object@scale.data[[i]]),]
      gene_df.sub$Cluster = object@clusters[rownames(object@scale.data[[i]])]
      title = names(object@scale.data)[i]
    } else {
      gene_df.sub = gene_df
      gene_df.sub$Cluster = object@clusters
      title = 'All Datasets'
    }
    max_v = max(gene_df.sub["gene"], na.rm = T)
    min_v = min(gene_df.sub["gene"], na.rm = T)
    midpoint = (max_v - min_v) / 2
    plot_i = (ggplot(gene_df.sub,aes_string(x="Cluster",y="gene",fill="Cluster"))+
                geom_boxplot(position="dodge",width=0.4,outlier.shape=NA, alpha=0.7) +
                geom_violin(position="dodge", alpha=0.7)+
                ggtitle(title))
    gene_plots[[i]] = plot_i + theme(legend.position="none") + labs(y=gene)
    if (i == 1 & !by.dataset) {break}
  }
  if (return.plots) {
    return(gene_plots)
  } else {
    for (i in 1:length(gene_plots)) {
      print(gene_plots[[i]])
    }
  }
}

#' Plot t-SNE coordinates per dataset, colored by expression of specified gene
#'
#' @param object analogizer object. Should call run_tSNE before calling.
#' @param gene Gene for which to plot relative expression. 
#' @param methylation_indices Indices of datasets in object with methylation data. 
#' @param return.plots Return ggplot objects instead of printing directly
#' @param pt.size Point size for plots
#' @param points.only Remove axes when plotting t-sne coordinates
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggplot geom_point aes_string scale_color_gradient2 ggtitle
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }
plot_gene = function(object, gene, methylation_indices=NULL, 
                     return.plots=F,pt.size=0.1,min.clip=0,
                     max.clip=1,points.only=F)
{
  gene_vals = c()
  for (i in 1:length(object@norm.data))
  {
    if (i %in% methylation_indices) {
      tmp = object@norm.data[[i]][gene,]
      max_v = quantile(tmp,probs=max.clip, na.rm = T)
      min_v = quantile(tmp,probs=min.clip, na.rm = T)
      tmp[tmp < min_v & !is.na(tmp)]=min_v
      tmp[tmp > max_v & !is.na(tmp)]=max_v
      gene_vals = c(gene_vals,tmp)
      
    } else {
      if (gene %in% rownames(object@norm.data[[i]]))
      {
        gene_vals_int = log2(10000*object@norm.data[[i]][gene,] + 1)
      }
      else
      {
        gene_vals_int = rep(list(0), ncol(object@norm.data[[i]]))
        names(gene_vals_int) = colnames(object@norm.data[[i]])
      }
      gene_vals = c(gene_vals, gene_vals_int)
    }
  }
  
  gene_df = data.frame(object@tsne.coords)
  rownames(gene_df)=names(object@clusters)
  gene_df$Gene = as.numeric(gene_vals[rownames(gene_df)])
  colnames(gene_df)=c("tSNE1","tSNE2","gene")
  gene_plots = list()
  for (i in 1:length(object@norm.data))
  {
    gene_df.sub = gene_df[rownames(object@scale.data[[i]]),]
    max_v = max(gene_df.sub["gene"],na.rm = T)
    min_v = min(gene_df.sub["gene"],na.rm = T)
    
    midpoint = (max_v - min_v) / 2
    plot_i = (ggplot(gene_df.sub,aes_string(x="tSNE1",y="tSNE2",color="gene"))+geom_point(size=pt.size)+
                scale_color_gradient(low="yellow",high="red",
                                     limits=c(min_v, max_v)) + labs(y=gene) +
                ggtitle(names(object@scale.data)[i]))
    gene_plots[[i]] = plot_i
  }
  if (points.only)
  {
    for (i in 1:length(gene_plots)) {
      gene_plots[[i]] = gene_plots[[i]] + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                axis.title.x=element_blank(),
                                                axis.title.y=element_blank(),legend.position="none",
                                                panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                                panel.grid.minor=element_blank(),plot.background=element_blank(),plot.title=element_blank())
    }
  }
  if (return.plots) {
    return(gene_plots)
  } else {
    for (i in 1:length(gene_plots)) {
      print(gene_plots[[i]])
    }
  }
}

#' Plot expression of multiple genes, each on a separate page.
#'
#' @param object analogizer object. Should call run_tSNE before calling.
#' @param genes vector of gene names.
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggplot geom_point aes_string scale_color_gradient2 ggtitle
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }
plot_genes = function(object,genes)
{
  for (i in 1:length(genes))
  {
    print(genes[i])
    plot_gene(object,genes[i])
  }
}

#Function takes in a list of DGEs, with gene rownames and cell colnames, and merges them into a single DGE.
MergeSparseDataAll<-function (datalist,library.names=NULL) {
  
  #use summary to convert the sparse matrices a and b into three-column indexes where i are the row numbers, j are the column numbers, and x are the nonzero entries
  a = datalist[[1]]
  allGenes=rownames(a)
  if (!is.null(library.names)) {
    allCells=paste0(library.names[1],"_",colnames(a))
  } else {
    allCells = colnames(a)
  }
  as = summary(a)
  for (i in 2:length(datalist)) {
    b = datalist[[i]]
    
    bs= summary(b)
    
    # Now, alter the indexes so that the two 3-column matrices can be properly merged.  First, make the a and b column numbers non-overlapping.
    bs[,2] = bs[,2] + max(as[,2])
    
    #Next, change the row (gene) indexes so that they index on the union of the gene sets, so that proper merging can occur.
    
    allGenesnew=union(allGenes, rownames(b))
    if (!is.null(library.names)) {
      cellnames = paste0(library.names[i],"_",colnames(b))
    } else {
      cellnames = colnames(b)
    }
    allCells=c(allCells,cellnames)
    idx=match(allGenes,allGenesnew)
    newgenesa = idx[as[,1]]
    as[,1] = newgenesa
    idx=match(rownames(b),allGenesnew)
    newgenesb = idx[bs[,1]]
    bs[,1] = newgenesb
    
    #Now bind the altered 3-column matrices together, and convert into a single sparse matrix.
    as = rbind(as,bs)
    allGenes=allGenesnew
  }
  M=sparseMatrix(i=as[,1],j=as[,2],x=as[,3],dims=c(length(allGenes),length(allCells)),dimnames=list(allGenes,allCells))
  return(M)
}

#' Create a Seurat object containing the data from an Analogizer object.
#'
#' @param object analogizer object.
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }
AnalogizerToSeurat<-function(object, need.sparse=T)  {
  if (!require("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  nms = names(object@H)
  if (need.sparse) {
    object@raw.data = lapply(object@raw.data, function(x){Matrix(as.matrix(x), sparse=T)})
    object@norm.data = lapply(object@norm.data, function(x){Matrix(as.matrix(x), sparse=T)})
  }
  raw.data = MergeSparseDataAll(object@raw.data,nms)
  norm.data = MergeSparseDataAll(object@norm.data, nms)
  
  scale.data = do.call(rbind, object@scale.data)
  rownames(scale.data) = colnames(norm.data)
  inmf.obj = new(Class = "dim.reduction", gene.loadings = t(object@W),
                 cell.embeddings = object@H.norm, key = "iNMF")
  tsne.obj = new(Class = "dim.reduction", cell.embeddings = object@tsne.coords,
                 key = "tSNE_")
  rownames(tsne.obj@cell.embeddings) = rownames(scale.data)
  rownames(inmf.obj@cell.embeddings) = rownames(scale.data)
  colnames(tsne.obj@cell.embeddings) = paste0("tSNE_", 1:2)
  new.seurat = CreateSeuratObject(raw.data)
  new.seurat = NormalizeData(new.seurat)
  new.seurat@scale.data = t(scale.data)
  new.seurat@dr$tsne = tsne.obj
  new.seurat@dr$inmf = inmf.obj
  
  new.seurat= SetIdent(new.seurat,ident.use = as.character(object@clusters))
  return(new.seurat)
}

#' Perform fast and memory-efficient normalization operation on sparse matrix data.
#'
#' @param A Sparse matrix DGE.
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }
Matrix.column_norm <- function(A){
  if (class(A)[1] == "dgTMatrix") {
    temp = summary(A)
    A = sparseMatrix(i=temp[,1],j=temp[,2],x=temp[,3])
  }
  A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
  return(A)
}

# k is the number of nearest neighbors (for geometry preservation)
# sigma is the bandwidth parameter for the RBF kernel (geometry preservation)
# alpha sets the tradeoff between reconstructing geometry and aligning corresponding clusters
spectral_alignment=function(object,k=30,alpha=1,sigma=NULL,neigen=NULL)
{
  print("Building graph")
  #knn1 = get.knn(scale(object@H[[1]],center=F,scale=T),k=k)
  #knn2 = get.knn(scale(object@H[[2]],center=F,scale=T),k=k)
  #knn1 = get.knn(object@H[[1]],k=k)
  #knn2 = get.knn(object@H[[2]],k=k)
  obj = CreateSeuratObject(t(analogy@scale.data[[1]]))
  knn1 = BuildSNN(obj,distance.matrix=as.matrix(dist(analogy@H[[1]])),k.param=30)@snn
  obj = CreateSeuratObject(t(analogy@scale.data[[2]]))
  knn2 = BuildSNN(obj,distance.matrix=as.matrix(dist(analogy@H[[2]])),k.param=30)@snn
  rm(obj)
  adj_mat = replicate(length(object@clusters),object@clusters)
  adj_mat = (adj_mat==t(adj_mat))
  adj_mat = as.matrix(1*adj_mat)
  colnames(adj_mat)=rownames(adj_mat)
  adj_mat[rownames(object@scale.data[[1]]),rownames(object@scale.data[[1]])]=0
  adj_mat[rownames(object@scale.data[[2]]),rownames(object@scale.data[[2]])]=0
  
  if (is.null(sigma))
  {
    sigma=mean(as.matrix(dist(object@H[[1]])))
    #sigma=mean(as.matrix(dist(scale(object@H[[1]],center=F,scale=T))))  
  }
  
  adj_mat_dists = matrix(0,nrow = nrow(adj_mat),ncol=ncol(adj_mat))
  dimnames(adj_mat_dists)=dimnames(adj_mat)
  for(i in 1:nrow(object@scale.data[[1]]))
  {
    adj_mat_dists[rownames(object@scale.data[[1]])[i],rownames(object@scale.data[[1]])[knn1$nn.index[i,]]] = exp(-1*knn1$nn.dist[i,]/sigma)
  }
  for(i in 1:nrow(object@scale.data[[2]]))
  {
    adj_mat_dists[rownames(object@scale.data[[2]])[i],rownames(object@scale.data[[2]])[knn2$nn.index[i,]]] = exp(-1*knn2$nn.dist[i,]/sigma)
  }
  adj_mat = adj_mat/sum(adj_mat)*sum(adj_mat_dists)
  g = graph_from_adjacency_matrix(adj_mat_dists+alpha*adj_mat,mode="undirected",weighted=T)
  print("Finding eigenvectors of the Laplacian")
  if (is.null(neigen))
  {
    neigen = length(levels(object@clusters))
  }
  res = eigs(laplacian_matrix(g),k=neigen,which="LM",sigma=0)
  rownames(res$vectors)=names(V(g))
  object@H.norm = res$vectors[names(object@clusters),]
  return(object)
}

#SNN version
spectral_alignment=function(object,k=30,alpha=1,sigma=NULL,neigen=NULL)
{
  print("Building graph")
  #knn1 = get.knn(scale(object@H[[1]],center=F,scale=T),k=k)
  #knn2 = get.knn(scale(object@H[[2]],center=F,scale=T),k=k)
  #knn1 = get.knn(object@H[[1]],k=k)
  #knn2 = get.knn(object@H[[2]],k=k)
  obj = CreateSeuratObject(t(object@scale.data[[1]]))
  knn1 = BuildSNN(obj,distance.matrix=as.matrix(dist(object@H[[1]])),k.param=k)@snn
  obj = CreateSeuratObject(t(object@scale.data[[2]]))
  knn2 = BuildSNN(obj,distance.matrix=as.matrix(dist(object@H[[2]])),k.param=k)@snn
  rm(obj)
  adj_mat = replicate(length(object@clusters),object@clusters)
  adj_mat = (adj_mat==t(adj_mat))
  adj_mat = as.matrix(1*adj_mat)
  colnames(adj_mat)=rownames(adj_mat)
  adj_mat[rownames(object@scale.data[[1]]),rownames(object@scale.data[[1]])]=0
  adj_mat[rownames(object@scale.data[[2]]),rownames(object@scale.data[[2]])]=0
  
  if (is.null(sigma))
  {
    sigma=mean(as.matrix(dist(object@H[[1]])))
    #sigma=mean(as.matrix(dist(scale(object@H[[1]],center=F,scale=T))))  
  }
  
  adj_mat_dists = bdiag(knn1,knn2)
  dimnames(adj_mat_dists)=dimnames(adj_mat)
  rm(knn1)
  rm(knn2)
  
  adj_mat = adj_mat/sum(adj_mat)*sum(adj_mat_dists)
  g = graph_from_adjacency_matrix(adj_mat_dists+alpha*adj_mat,mode="undirected",weighted=T)
  print("Finding eigenvectors of the Laplacian")
  if (is.null(neigen))
  {
    neigen = length(levels(object@clusters))
  }
  res = eigs(laplacian_matrix(g),k=neigen,which="LM",sigma=0)
  rownames(res$vectors)=names(V(g))
  object@H.norm = res$vectors[names(object@clusters),]
  return(object)
}

#' Perform fast and memory-efficient data scaling operation on sparse matrix data.
#'
#' @param object Sparse matrix DGE.
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }
scaleNotCenter_sparse<-function (object, cells = NULL)
{
  if (is.null(cells)) {
    cells = lapply(1:length(object@raw.data), function(i) {
      1:ncol(object@raw.data[[i]])
    })
  }
  object@scale.data = lapply(1:length(object@norm.data), function(i) {
    scale(Sparse_transpose(object@norm.data[[i]][object@var.genes, ]), center = F,
          scale = T)
  })
  names(object@scale.data) = names(object@norm.data)
  for (i in 1:length(object@scale.data)) {
    object@scale.data[[i]][is.na(object@scale.data[[i]])] = 0
    rownames(object@scale.data[[i]]) = colnames(object@raw.data[[i]])
    colnames(object@scale.data[[i]]) = object@var.genes
  }
  return(object)
}

scalar1 <- function(x, center=F) {
  if (center) {
    return(scale(x)) 
  }
  return(x / sqrt(sum(x^2)))
}

Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

#' Builds a shared nearest factor graph to jointly cluster cells, then quantile normalizes corresponding clusters.
#'
#' @param object analogizer object. Should run optimizeALS before calling.
#' @param knn_k Number of nearest neighbors for within-dataset knn graph
#' @param k2 Horizon parameter for shared nearest factor graph. Distances to all but the k2 nearest neighbors are set to 0 (cuts down on memory usage for very large graphs).
#' @param prune.thresh Minimum allowed edge weight. Any edges below this are removed (given weight 0)
#' @param ref_dataset Name of dataset to use as a "reference" for normalization. By default,
#' the dataset with the largest number of cells is used.
#' @param min_cells Minimum number of cells to consider a cluster shared across datasets.
#' @param nstart Number of times to perform Louvain community detection with different random starts
#' @param quantiles Number of quantiles to use for quantile normalization
#' @param resolution Controls the number of communities detected. Higher resolution=more communities.
#' @param dims.use Indices of factors to use for shared nearest factor determination
#' @param dist.use Distance metric to use in calculating nearest neighbors
#' @param center Centers the data when scaling factors (useful for less sparse modalities like methylation data)
#' @param small.clust.thresh Extracts small clusters smaller than this threshold before regularr alignment
#' @param id.number Number to use for identifying edge file (when running in parallel)
#' @param print.mod Print modularity output from clustering algorithm
#' @param print_align_summary Print summary of clusters which did not align normally
#' 
#' @return analogizer object
#' @export
#' @importFrom RANN.L1 nn2
#' @importFrom FNN get.knnx
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }

quantile_align_SNF<-function(object,knn_k=20,k2=500,prune.thresh=0.2,ref_dataset=NULL,min_cells=2,
                             quantiles=50,nstart=10, resolution = 1, dims.use = 1:ncol(object@H[[1]]),
                             dist.use='CR', center=F, small.clust.thresh=0, 
                             id.number=NULL, print.mod=F, print_align_summary=T) {
  if (is.null(ref_dataset)) {
    ns = sapply(object@scale.data, nrow)
    ref_dataset = names(object@scale.data)[which.max(ns)]
  }
  if (is.null(id.number)) {
    set.seed(NULL)
    id.number = sample(1000000:9999999, 1)
  }
  # recompute SNF if parameters have changed
  if (!isTRUE(object@parameters[['knn_k']] == knn_k) |
      !isTRUE(object@parameters[['k2']] == k2) |
      !isTRUE(object@parameters[['dist.use']] == dist.use) |
      !isTRUE(object@parameters[['SNF_center']] == center) |
      !isTRUE(identical(object@parameters[['dims.use']], dims.use)) |
      !isTRUE(object@parameters[['small.clust.thresh']] == small.clust.thresh)) {
    print('Recomputing shared nearest factor space')
    snf = SNF(object,knn_k=knn_k,k2=k2, dist.use=dist.use, center = center, 
              dims.use=dims.use, small.clust.thresh=small.clust.thresh)
  } else {
    snf = object@snf
  }
  cell.names = unlist(lapply(object@scale.data,rownames))
  idents = snf$idents
  
  idents.rest = SLMCluster(edge = snf$out.summary, nstart = nstart, R = resolution, 
                           prune.thresh = prune.thresh, id.number = id.number,
                           print.mod=print.mod)
  names(idents.rest) = setdiff(cell.names,snf$cells.cl)
  # Especially when datasets are large, SLM generates a fair number of singletons.  
  # To assign these to a cluster, take the mode of the cluster assignments of within-dataset neighbors
  if (min(table(idents.rest)) == 1) {
    idents.rest = assign.singletons(object, idents.rest, center = center)
  }
  idents[names(idents.rest)] = as.character(idents.rest)
  idents = factor(idents)
  names(idents) = cell.names
  
  Hs = object@H
  cs = cumsum(c(0,unlist(lapply(object@H,nrow))))
  clusters = lapply(1:length(Hs),function(i){
    idx = cs[i] + 1:nrow(Hs[[i]])
    return(idents[idx])
  })
  names(clusters) = names(object@H)
  dims = ncol(object@H[[ref_dataset]])
  
  too.few = rep(list(c()), length(Hs))
  names(too.few) = names(Hs)
  for (k in 1:length(Hs)) {
    for (i in 1:dims) {
      for (j in levels(idents)) {
        if (sum(clusters[[ref_dataset]] == j, na.rm = T) < min_cells |
            sum(clusters[[k]] == j, na.rm = T) < min_cells) {
          too.few[[names(Hs)[k]]] = c(too.few[[names(Hs)[k]]], j)
          next
        }
        if (sum(clusters[[k]] == j, na.rm = T) == 1) {
          Hs[[k]][clusters[[k]] == j, i] = mean(Hs[[ref_dataset]][clusters[[ref_dataset]] ==
                                                                    j, i])
          too_few[[names(Hs)[k]]] = c(too_few[[names(Hs)[k]]], j)
          next
        }
        q2 = quantile(Hs[[k]][clusters[[k]] == j, i],
                      seq(0, 1, by = 1/quantiles), na.rm = T)
        q1 = quantile(Hs[[ref_dataset]][clusters[[ref_dataset]] ==
                                          j, i], seq(0, 1, by = 1/quantiles), na.rm = T)
        if (sum(q1) == 0 | sum(q2) == 0 | length(unique(q1)) <
            2 | length(unique(q2)) < 2) {
          new_vals = rep(0, sum(clusters[[k]] == j))
        }
        else {
          warp_func = approxfun(q2, q1)
          new_vals = warp_func(Hs[[k]][clusters[[k]] ==
                                         j, i])
        }
        if (anyNA(new_vals)) {
          stop('Select lower resolution; too many communities detected.')
        }
        Hs[[k]][clusters[[k]] == j, i] = new_vals
      }
    }
  }
  if (print_align_summary & length(unlist(too.few)) > 0) {
    print('Summary:')
    for (i in 1:length(Hs)) {
      print(paste('In dataset', names(Hs)[i], 'these clusters did not align normally (too few cells):'))
      print(unique(too.few[[names(Hs)[i]]]))
    }
  }
  object@H.norm = Reduce(rbind, Hs)
  object@clusters = idents
  object@snf = snf
  # set parameters
  object@parameters$ref_dataset = ref_dataset
  object@parameters$knn_k = knn_k
  object@parameters$k2 = k2
  object@parameters$prune.thresh = prune.thresh
  object@parameters$min_cells = min_cells
  object@parameters$dims.use = dims.use
  object@parameters$dist.use = dist.use
  object@parameters$SNF_center = center
  object@parameters$small.clust.thresh = small.clust.thresh
  object@parameters$resolution = resolution
  return(object)
}

SNF<-function (object, dims.use = 1:ncol(object@H[[1]]), knn_k = 20, 
               k2 = 500, dist.use = "CR", center = F, 
               small.clust.thresh=knn_k) 
{
  NN.maxes = do.call(rbind, lapply(1:length(object@H), function(i) {
    sc = scale(object@H[[i]], center = center, scale = T)
    maxes = factor(apply(sc[, dims.use], 1, which.max), levels = 1:ncol(sc))
    if (dist.use == "CR") {
      norm = t(apply(object@H[[i]][, dims.use], 1, scalar1))
      if (any(is.na(norm))) {
        stop("Unable to normalize loadings for all cells; some cells\n             
             loading on no selected factors.")
      }
      }
    else {
      norm = object@H[[i]][, dims.use]
    }
    knn.idx = get.knn(norm, knn_k, algorithm = dist.use)$nn.index
    t(apply(knn.idx, 1, function(q) {
      table(maxes[q])
    }))
  }))
  rownames(NN.maxes) = unlist(lapply(object@H, rownames))
  # extract small clusters 
  if (small.clust.thresh > 0) {
    print(paste0('Isolating small clusters with fewer than ', small.clust.thresh,
                 ' members'))
  }
  max.val = factor(apply(NN.maxes,1,which.max))
  names(max.val) = rownames(NN.maxes)
  idents = rep("NA",nrow(NN.maxes))
  names(idents) = rownames(NN.maxes)
  cl = levels(max.val)[which(table(max.val) < small.clust.thresh)]
  cells.cl = names(max.val)[which(max.val %in% cl)]
  idents[cells.cl]= paste0("F",as.character(max.val[cells.cl]))
  
  nn.obj <- nn2(NN.maxes[setdiff(rownames(NN.maxes),cells.cl),], k = k2)
  out.snn = 1 - (nn.obj$nn.dists/(2 * knn_k))
  out.summary = matrix(ncol = 3, nrow = (ncol(out.snn) * nrow(out.snn)))
  
  counter = 1
  for (i in 1:nrow(out.snn)) {
    for (j in 1:ncol(out.snn)) {
      out.summary[counter, ] = c(i, nn.obj$nn.idx[i, j], 
                                 out.snn[i, j])
      counter = counter + 1
    }
  }
  out.summary[out.summary[, 1] == out.summary[, 2], 3] = 0
  out.summary[, 1] = out.summary[, 1] - 1
  out.summary[, 2] = out.summary[, 2] - 1
  
  # idents returned here only contain values for small clusters
  return(list(cells.cl=cells.cl, idents=idents, 
              out.summary=out.summary))
  }

assign.singletons<-function(object,idents,k.use = 15, center=F) {
  singleton.clusters = names(table(idents))[which(table(idents)==1)]
  singleton.cells = names(idents)[which(idents %in% singleton.clusters)]
  if (length(singleton.cells)>0) {
    singleton.list = lapply(object@H,function(H){
      
      return(intersect(rownames(H),singleton.cells))
    })
    out = unlist(lapply(1:length(object@H),function(d){
      H = object@H[[d]]
      H = t(apply(H,1,scalar1, center=center))
      cells.use = singleton.list[[d]]
      
      if (length(cells.use)>0) {
        knn.idx = get.knnx(H,matrix(H[cells.use,],ncol=ncol(H)),k = k.use,algorithm="CR")$nn.index
        o= sapply(1:length(cells.use),function(i){
          ids = idents[knn.idx[i,]]
          Mode(ids[which(!(ids %in% singleton.clusters))])
        })
        
        return(o)
      } else {
        return(NULL)
      }
    }))
    names(out)= unlist(singleton.list)
    idx = names(idents) %in% names(out)
    idents[idx]<-out
    idents = droplevels(idents)
  }
  return(idents)
}


SLMCluster<-function(edge,prune.thresh=0.2,nstart=100,iter.max=10,algorithm=1,R=1,
                     ModularityJarFile="",random.seed=1, id.number=NULL, print.mod=F) {
  
  #This is code taken from Seurat for preparing data for modularity-based clustering.
  edge = edge[which(edge[,3]>prune.thresh),]
  ident=modclust(edge=edge,resolution=R,n.iter= iter.max,n.start=nstart, id.number = id.number,
                 random.seed=random.seed,ModularityJarFile=ModularityJarFile, print.mod=print.mod)
  return(ident)
}
modclust <- function(edge, modularity=1,resolution,n.start=100,n.iter=25,random.seed=1,
                     id.number=NULL, algorithm=1,ModularityJarFile, print.mod=F){
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
  analogizer.dir <- system.file(package = "Analogizer")
  ModularityJarFile <- paste0(analogizer.dir, "/java/ModularityOptimizer.jar")
  command <- paste("java -jar", ModularityJarFile, edge_file, output_file,modularity, resolution, algorithm, n.start,
                   n.iter, random.seed, as.numeric(print.mod), sep = " ")
  print ("Starting SLM")
  ret = system(command, wait = TRUE)
  if (ret != 0) {
    stop(paste0(ret, " exit status from ", command))
  }
  unlink(edge_file)
  ident.use <- factor(read.table(file = output_file, header = FALSE, sep = "\t")[, 1])
  
  return(ident.use)
}

#' Makes a riverplot to show how separate cluster assignments from two datasets map onto a joint clustering.
#'
#' @param object analogizer object. Should run quantile_align_SNF before calling.
#' @param cluster1 Cluster assignments for dataset 1. Note that cluster names should be distinct across datasets.
#' @param cluster2 Cluster assignments for dataset 2. Note that cluster names should be distinct across datasets.
#' @param cluster_consensus Optional external consensus clustering (to use instead of object clusters)
#' @param min.frac Minimum fraction of cluster for edge to be shown
#' @param min.cells Minumum number of cells for edge to be shown 
#' @param river.yscale y-scale to pass to riverplot
#' @param river.lty lty to pass to riverplot
#' @param river.node_margin node_margin to pass to riverplot
#' @param label.cex Size of text labels
#' @param label.col Color of text labels
#' @param lab.srt Angle of text labels
#' @param node.order Order of clusters in each set (list with three vectors of ordinal numbers)
#' @export
#' @importFrom riverplot makeRiver
#' @importFrom riverplot plot.riverplot
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy@var.genes = c(1,2,3,4)
#' analogy = scaleNotCenter(analogy)
#' }
make_riverplot<-function(object, cluster1,cluster2, cluster_consensus=NULL,min.frac = 0.05,min.cells=10,
                         river.yscale = 1,river.lty=0,river.node_margin = 0.1,label.cex = 1,
                         label.col='black',lab.srt = 0,node.order = "auto") {
  cluster1 = droplevels(cluster1)
  cluster2 = droplevels(cluster2)
  if (is.null(cluster_consensus)) {
    cluster_consensus=droplevels(object@clusters)
  }
  # Make cluster names unique if necessary
  if(length(intersect(levels(cluster1),levels(cluster2))) > 0 | 
     length(intersect(levels(cluster1),levels(cluster_consensus)))>0 | 
     length(intersect(levels(cluster2),levels(cluster_consensus)))>0)
  {
    print("Duplicate cluster names detected. Adding 1- and 2- to make unique names.")
    cluster1 = mapvalues(cluster1,from=levels(cluster1),to=paste("1",levels(cluster1),sep="-"))
    cluster2 = mapvalues(cluster2,from=levels(cluster2),to=paste("2",levels(cluster2),sep="-"))
  }
  # set node order 
  if (node.order == "auto") {
    tab.1 = table(cluster1,cluster_consensus[names(cluster1)])
    tab.1 = sweep(tab.1,1,rowSums(tab.1),"/")
    tab.2 = table(cluster2,cluster_consensus[names(cluster2)])
    tab.2 = sweep(tab.2, 1, rowSums(tab.2),"/")
    whichmax.1 = apply(tab.1,1,which.max)
    whichmax.2 = apply(tab.2,1,which.max)
    ord.1 = order(whichmax.1)
    ord.2 = order(whichmax.2)
    cluster1 = factor(cluster1,levels=levels(cluster1)[ord.1])
    cluster2 = factor(cluster2,levels=levels(cluster2)[ord.2])
  } else {
    if (is.list(node.order)) {
      cluster1 = factor(cluster1,levels=levels(cluster1)[node.order[[1]]])
      cluster_consensus = factor(cluster_consensus,levels=levels(cluster_consensus)[node.order[[2]]])
      cluster2 = factor(cluster2,levels=levels(cluster2)[node.order[[3]]])
    } 
  }
  cluster1 = cluster1[!is.na(cluster1)]
  cluster2 = cluster2[!is.na(cluster2)]
  nodes1 = levels(cluster1)[table(cluster1) > 0]
  nodes2 = levels(cluster2)[table(cluster2) > 0]
  nodes_middle = levels(cluster_consensus)[table(cluster_consensus) > 0]
  node_Xs = c(rep(1, length(nodes1)), rep(2, length(nodes_middle)), 
              rep(3, length(nodes2)))
  
  # first set of edges
  edge_list = list()
  for (i in 1:length(nodes1)) {
    temp = list()
    i_cells = names(cluster1)[cluster1 == nodes1[i]]
    for (j in 1:length(nodes_middle)) {
      if(length(which(cluster_consensus[i_cells] == nodes_middle[j]))/length(i_cells) > min.frac & 
         length(which(cluster_consensus[i_cells] == nodes_middle[j])) > min.cells) {
        temp[[nodes_middle[j]]] = sum(cluster_consensus[i_cells] == 
                                        nodes_middle[j])/length(cluster1)
      } 
    }
    edge_list[[nodes1[i]]] = temp
  }
  # second set of edges
  cluster3 = cluster_consensus[names(cluster2)]
  for (i in 1:length(nodes_middle)) {
    temp = list()
    i_cells = names(cluster3)[cluster3 == nodes_middle[i]]
    for (j in 1:length(nodes2)) {
      j_cells = names(cluster2)[cluster2 == nodes2[j]]
      if (length(which(cluster_consensus[j_cells] == nodes_middle[i]))/length(j_cells) > min.frac & 
          length(which(cluster_consensus[j_cells] == nodes_middle[i])) > min.cells) {
        if (!is.na(sum(cluster2[i_cells] == nodes2[j]))) {
          temp[[nodes2[j]]] = sum(cluster2[i_cells] == 
                                    nodes2[j])/length(cluster2)
        }
      }
    } 
    edge_list[[nodes_middle[i]]] = temp
  }
  # set cluster colors 
  node_cols = list()
  ggplotColors <- function(g) {
    d <- 360/g
    h <- cumsum(c(15, rep(d, g - 1)))
    hcl(h = h, c = 100, l = 65)
  }
  pal = ggplotColors(length(nodes1))
  for (i in 1:length(nodes1)) {
    node_cols[[nodes1[i]]] = list(col = pal[i],textcex=label.cex,textcol=label.col,srt= lab.srt)
  }
  pal = ggplotColors(length(nodes_middle))
  for (i in 1:length(nodes_middle)) {
    node_cols[[nodes_middle[i]]] = list(col = pal[i],textcex=label.cex,textcol=label.col,srt= lab.srt)
  }
  pal = ggplotColors(length(nodes2))
  for (i in 1:length(nodes2)) {
    node_cols[[nodes2[i]]] = list(col = pal[i],textcex=label.cex,textcol=label.col,srt= lab.srt)
  }
  # create nodes and riverplot object
  nodes = list(nodes1,nodes_middle,nodes2)
  node.limit= max(unlist(lapply(nodes,length)))
  
  node_Ys = lapply(1:length(nodes),function(i){
    seq(1,node.limit,by = node.limit/length(nodes[[i]]))  
    
  })
  rp = makeRiver(c(nodes1, nodes_middle, nodes2), edge_list, 
                 node_xpos = node_Xs, node_ypos = unlist(node_Ys),node_styles = node_cols)
  # prevent normal riverplot output being printed to console
  invisible(capture.output(riverplot(rp,yscale = river.yscale,lty= river.lty,
                                     node_margin = river.node_margin)))
  
}
#' Calculate adjusted Rand index for Analogizer clustering and external clustering.
#' Should run quantile_align_SNF first. 
#'
#' @param object analogizer object.
#' @param clusters.compare Clustering with which to compare (named vector)
#' @importFrom mclust adjustedRandIndex
#' @return Adjusted Rand index value
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy =
#' }
calcARI = function(object, clusters.compare) {
  if (length(clusters.compare) < length(object@clusters)) {
    print('Calculating ARI for subset of full cells')
  }
  return(adjustedRandIndex(object@clusters[names(clusters.compare)],
                           clusters.compare))
}

#' Calculate purity for Analogizer clustering and external clustering (base truth).
#' Should run quantile_align_SNF first. 
#'
#' @param object analogizer object.
#' @param classes.compare Clustering with which to compare (named vector)
#' @return Purity value
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy =
#' }
calcPurity = function(object, classes.compare) {
  if (length(classes.compare) < length(object@clusters)) {
    print('Calculating purity for subset of full cells')
  }
  clusters = object@clusters[names(classes.compare)]
  purity = sum(apply(table(classes.compare, clusters), 2, max)) / length(clusters)

  return(purity)
}

#' Construct an Analogizer object with a specified subset of cells or clusters.
#' Should only be called after tsne coordinates have been computed and alignment
#' performed. 
#'
#' @param object analogizer object.
#' @param clusters.use Clusters to extract
#' @param cells.use Vector of cell names to keep from any dataset
#' 
#' @return analogizer object
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' analogy = Analogizer(list(Y,Z))
#' analogy =
#' }
subsetAnalogizer<-function(object, clusters.use = NULL,cells.use = NULL) {
  if (!is.null(clusters.use)){
    cells.use = names(object@clusters)[which(object@clusters %in% clusters.use)]
    
  }
  raw.data = lapply(seq_along(object@raw.data),function(q){
    cells = intersect(cells.use,colnames(object@raw.data[[q]]))
    if (length(cells) > 0) {
      object@raw.data[[q]][,cells]
    } else {
      warning(paste0('Selected subset eliminates dataset ', names(object@raw.data)[q]))
      return(NULL)
    }
  })
  raw.data = raw.data[!sapply(raw.data,is.null)]
  nms = names(object@raw.data)[!sapply(raw.data, is.null)]
  a = Analogizer(raw.data)
  
  a@norm.data = lapply(1:length(a@raw.data),function(i){
    object@norm.data[[i]][,colnames(a@raw.data[[i]])]
    
  })
  a@scale.data = lapply(1:length(a@raw.data),function(i){
    object@scale.data[[i]][colnames(a@raw.data[[i]]),]
    
  })
  a@H = lapply(1:length(a@raw.data),function(i){
    object@H[[i]][colnames(a@raw.data[[i]]),]
  })
  a@clusters = object@clusters[unlist(lapply(a@H,rownames))]
  a@clusters = droplevels(a@clusters)
  a@tsne.coords = object@tsne.coords[names(a@clusters),]
  a@H.norm = object@H.norm[names(a@clusters),]
  a@W = object@W
  a@V = object@V
  a@var.genes = object@var.genes
  names(a@scale.data) = names(a@raw.data) = names(a@norm.data) = names(a@H) = nms
  return(a)
}
