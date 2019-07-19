# These are deprecated functions likely to be removed in future versions. 
# Documentation for these functions is incomplete.

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

#' Calculate alignment metric per factor.
#' 
#' @param object Analogizer object. Should run quantile_align_SNF before calling.
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
clusterLouvainJaccard = function(object, resolution = 0.1, k.param=30, n.iter = 10, n.start = 10,
                                 print.output = F, ...)
{
  temp.seurat = CreateSeuratObject(t(Reduce(rbind,object@scale.data)))
  temp.seurat@scale.data = t(Reduce(rbind,object@scale.data))
  rownames(object@H.norm)=colnames(temp.seurat@scale.data)
  temp.seurat@dr$NMF=new(Class="dim.reduction",cell.embeddings=object@H.norm,key="NMF")
  temp.seurat <- FindClusters(object = temp.seurat, reduction.type = "NMF", 
                              dims.use = 1:ncol(object@H.norm),force.recalc=T,
                              save.SNN = T,resolution=resolution,k.param=k.param,
                              n.iter = n.iter, n.start = n.start, print.output = print.output, ...)
  object@clusters = temp.seurat@ident
  return(object)
}

clusterLouvain = function(object, knn_k = 20, k2 = 500, prune.thresh = 0.2, 
                                         min_cells = 2, nstart = 10, resolution = 1, 
                                         dims.use = 1:ncol(object@H[[1]]), dist.use = "JAC", center = F, 
                                         small.clust.thresh = 0, id.number = NULL, print.mod = F, 
                                         print.align.summary = F)
{
  snf <- SNF(object,
             knn_k = knn_k, k2 = k2, dist.use = dist.use, center = center,
             dims.use = dims.use, small.clust.thresh = small.clust.thresh)
  cell.names <- unlist(lapply(object@scale.data, rownames))
  if(is(snf, "liger")){
    idents = snf@snf
  } else {
    idents <- snf[["idents"]]
  }
  
  idents.rest <- SLMCluster(
    edge = snf$out.summary, nstart = nstart, R = resolution,
    prune.thresh = prune.thresh, id.number = id.number,
    print.mod = print.mod
  )
  names(idents.rest) <- setdiff(cell.names, snf$cells.cl)
  # Especially when datasets are large, SLM generates a fair number of singletons.
  # To assign these to a cluster, take mode of the cluster assignments of within-dataset neighbors
  if (min(table(idents.rest)) == 1) {
    idents.rest <- assign.singletons(object, idents.rest, center = center)
  }
  idents[names(idents.rest)] <- as.character(idents.rest)
  idents <- factor(idents)
  names(idents) <- cell.names
  
  object@alignment.clusters <- idents
  object@clusters <- idents
  if(is(snf, "liger")){
    object@snf = snf@snf
  } else {
    object@snf <- snf
  }
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
#' 
#' @export
#' @examples
#' \dontrun{
#' Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
#' Z = matrix(c(1,2,3,4,5,6,7,6,5,4,3,2),nrow=4,byrow=T)
#' ligerex = createLiger(list(Y,Z))
#' ligerex@var.genes = c(1,2,3,4)
#' ligerex = scaleNotCenter(ligerex)
#' }
scaleNotCenter_sparse<-function (object, cells = NULL)
{
  print('This function has been deprecated and will soon be removed. Its functionality has been
         added to scaleNotCenter.')
  if (is.null(cells)) {
    cells = lapply(1:length(object@raw.data), function(i) {
      1:ncol(object@raw.data[[i]])
    })
  }
  object@scale.data = lapply(1:length(object@norm.data), function(i) {
    scale(sparse.transpose(object@norm.data[[i]][object@var.genes, ]), center = F,
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

newLiger <- function(raw.data, make.sparse = T, take.gene.union = F) {
  print('This function has been deprecated and will soon be removed. Please use createLiger instead.')
  return(createLiger(raw.data, make.sparse = make.sparse, take.gene.union = take.gene.union))
}