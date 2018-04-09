#' @importFrom plyr rbind.fill.matrix
NULL

#' Optimize objective function using block coordinate descent (alternating nonnegative least squares)
#'
#' @param E List of N data matrices. Samples (cells) are in rows, features (genes) in columns.
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
#' optimize_als(list(Y,Z),k=2,nrep=1)
#' }
optimize_als = function(E,k,lambda=5.0,thresh=1e-4,max_iters=100,nrep=20)
{
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
  return(list(H_m,W_m,V_m,run_stats))
}
