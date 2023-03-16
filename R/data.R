#' liger object of PBMC subsample data with Control and Stimulated datasets
#' @format \linkS4class{liger} object with two datasets named by "ctrl" and
#' "stim".
#' @source https://www.nature.com/articles/nbt.4042
"pbmc"

#' liger object of PBMC subsample data with plotting information available
#' @description This data was generated from data \code{"pbmc"} with default
#' parameter integration pipeline, with setting \code{k = 20, maxIter = 10}
#' for \code{\link{optimizeALS}}. UMAP was generated with \code{minDist = 0.5}.
#' Only the 2nd and 3rd factors, the normalized expression of the top 50
#' variable genes, and the clustering label were kept for visualization example.
#' @format \linkS4class{liger} object with two datasets named by "ctrl" and
#' "stim".
#' @source https://www.nature.com/articles/nbt.4042
"pbmcPlot"
