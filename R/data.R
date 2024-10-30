#' liger object of PBMC subsample data with Control and Stimulated datasets
#' @format \linkS4class{liger} object with two datasets named by "ctrl" and
#' "stim".
#' @source https://www.nature.com/articles/nbt.4042
#' @references Hyun Min Kang and et. al., Nature Biotechnology, 2018
"pbmc"

#' liger object of PBMC subsample data with plotting information available
#' @description This data was generated from data \code{"pbmc"} with default
#' parameter integration pipeline: normalize, selectGenes, scaleNotCenter,
#' runINMF, runCluster, runUMAP. To minimize the object size distributed with
#' the package, rawData and scaleData were removed. Genes are downsampled to
#' the top 50 variable genes, for smaller normData and \eqn{W} matrix.
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

#' Data frame for example marker DEG test result
#' @description
#' The data frame is the direct output of marker detection DEG test applied on
#' example dataset which can be loaded with \code{data("pbmc")}. The DEG test
#' was done with:
#' ```
#' defaultCluster(pbmc) <- pbmcPlot$leiden_cluster
#' deg.marker <- runMarkerDEG(
#'     pbmc,
#'     minCellPerRep = 5
#' )
#' ````
#' The result is for the marker detection test for 8 clusters in the dataset by
#' comparing each cluster against all other clusters.
#' @seealso [runMarkerDEG()]
#' @format data.frame object of 1992 rows with columns:
#' \itemize{
#' \item feature: gene names, 249 unique genes repeated 8 times for the tests
#' done for 8 clusters.
#' \item group: cluster names, 8 unique cluster names, dividing the tests.
#' \item logFC: log fold change of the gene expression between the cluster of
#' interest against all other clusters.
#' \item pval: p-value of the DEG test.
#' \item padj: adjusted p-value of the DEG test.
#' \item pct_in: percentage of cells in the cluster of interest expressing the
#' gene.
#' \item pct_out: percentage of cells in all other clusters expressing the gene.
#' }
"deg.marker"

#' Data frame for example pairwise DEG test result
#' @description
#' The data frame is the direct output of pairwise DEG test applied on example
#' dataset which can be loaded with \code{data("pbmc")}. The DEG test was done
#' with:
#' ```
#' defaultCluster(pbmc) <- pbmcPlot$leiden_cluster
#' degTest <- runPairwiseDEG(
#'     pbmc,
#'     groupTest = "stim",
#'     groupCtrl = "ctrl",
#'     variable1 = "dataset",
#'     splitBy = "defaultCluster"
#' )`
#' ```
#' The result is for the DEG test split for each cluster in the dataset, and
#' within each cluster, compare the cells from "stim" against the cells from
#' "ctrl".
#' @seealso [runPairwiseDEG()]
#' @format data.frame object of 1743 rows with columns:
#' \itemize{
#' \item feature: gene names, 249 unique genes repeated 7 times for the tests
#' done for 7 clusters. (1 less cluster than in \code{\link{deg.marker}} due to
#' too tiny sample size in the smallest cluster)
#' \item group: cluster names, 7 unique cluster names, dividing the tests.
#' \item logFC: log fold change of the gene expression between the condition of
#' interest against the control condition.
#' \item pval: p-value of the DEG test.
#' \item padj: adjusted p-value of the DEG test.
#' \item pct_in: percentage of cells in the condition of interest expressing the
#' gene.
#' \item pct_out: percentage of cells in the control condition expressing the
#' gene.
#' }
"deg.pw"
