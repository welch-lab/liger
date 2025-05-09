% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ATAC.R
\name{linkGenesAndPeaks}
\alias{linkGenesAndPeaks}
\title{Linking genes to putative regulatory elements}
\usage{
linkGenesAndPeaks(
  object,
  useDataset,
  pathToCoords,
  useGenes = NULL,
  method = c("spearman", "pearson", "kendall"),
  alpha = 0.05,
  verbose = getOption("ligerVerbose", TRUE),
  path_to_coords = pathToCoords,
  genes.list = useGenes,
  dist = method
)
}
\arguments{
\item{object}{A \linkS4class{liger} object, with datasets that is of
\linkS4class{ligerATACDataset} class in the \code{datasets} slot.}

\item{useDataset}{Name of one dataset, with both normalized gene expression
and normalized peak counts available.}

\item{pathToCoords}{Path tothe gene coordinates file, usually a BED file.}

\item{useGenes}{Character vector of gene names to be tested. Default
\code{NULL} uses all genes available in \code{useDataset}.}

\item{method}{Choose the type of correlation to calculate, from
\code{"spearman"}, \code{"pearson"} and \code{"kendall"}. Default
\code{"spearman"}}

\item{alpha}{Numeric, significance threshold for correlation p-value.
Peak-gene correlations with p-values below this threshold are considered
significant. Default \code{0.05}.}

\item{verbose}{Logical. Whether to show information of the progress. Default
\code{getOption("ligerVerbose")} or \code{TRUE} if users have not set.}

\item{path_to_coords, genes.list, dist}{\bold{Deprecated}. See Usage section
for replacement.}
}
\value{
A sparse matrix with peak names as rows and gene names as columns,
with each element indicating the correlation between peak i and gene j, 0 if
the gene and peak are not significantly linked.
}
\description{
Evaluate the relationships between pairs of genes and peaks
based on specified distance metric. Usually used for inferring the
correlation between gene expression and imputed peak counts for datasets
without the modality originally (i.e. applied to \code{\link{imputeKNN}}
result).
}
\examples{
\donttest{
if (requireNamespace("RcppPlanc", quietly = TRUE) &&
    requireNamespace("GenomicRanges", quietly = TRUE) &&
    requireNamespace("IRanges", quietly = TRUE) &&
    requireNamespace("psych", quietly = TRUE)) {
    bmmc <- normalize(bmmc)
    bmmc <- selectGenes(bmmc)
    bmmc <- scaleNotCenter(bmmc)
    bmmc <- runINMF(bmmc, miniBatchSize = 100)
    bmmc <- alignFactors(bmmc)
    bmmc <- normalizePeak(bmmc)
    bmmc <- imputeKNN(bmmc, reference = "atac", queries = "rna")
    corr <- linkGenesAndPeaks(
        bmmc, useDataset = "rna",
        pathToCoords = system.file("extdata/hg19_genes.bed", package = "rliger")
    )
}
}
}
\seealso{
\code{\link{imputeKNN}}
}
