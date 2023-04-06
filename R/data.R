#' liger object of PBMC subsample data with Control and Stimulated datasets
#' @format \linkS4class{liger} object with two datasets named by "ctrl" and
#' "stim".
#' @source https://www.nature.com/articles/nbt.4042
#' @references Hyun Min Kang and et. al., Nature Biotechnology, 2018
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
#' @references Hyun Min Kang and et. al., Nature Biotechnology, 2018
"pbmcPlot"

#' liger object of bone marrow subsample data with RNA and ATAC modality
#' @format \linkS4class{liger} object with two dataset named by "rna" and "atac"
#' @source https://www.nature.com/articles/s41587-019-0332-7
#' @references Jeffrey M. Granja and et. al., Nature Biotechnology, 2019
"bmmc"
