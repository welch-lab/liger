#setwd("/broad/macosko/jwelch/NMF/nnls/nmf")
#library(Rcpp)
#library(RcppArmadillo)
#
cells = as.matrix(read.csv("Hippocampus5000CellsScaled.csv",header=T,row.names=1))
nuclei = as.matrix(read.csv("Hippocampus5000NucleiScaled.csv",header=T,row.names=1))
#Sys.setenv(PKG_CXXFLAGS="-O3")
#sourceCpp("run_nmf.cpp")
Y = cells
source("/broad/macosko/jwelch/NMF/Analogizer.R")
#Y = matrix(c(1,2,3,4,5,6,7,8,9,10,11,12),nrow=4,byrow=T)
k = 3
m = nrow(Y)
n = ncol(Y)
#W <- matrix(abs(rnorm(m * k)), m, k)	
H <- matrix(abs(runif(k * n,min=min(Y),max=max(Y))), k, n)
res = solve_nnls(t(H),t(Y))

