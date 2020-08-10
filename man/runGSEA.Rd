% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/liger.R
\name{runGSEA}
\alias{runGSEA}
\title{Analyze biological interpretations of metagene}
\usage{
runGSEA(object, gene_sets = c(), mat_w = T, mat_v = 0, custom_gene_sets = c())
}
\arguments{
\item{object}{\code{liger} object.}

\item{gene_sets}{A list of the Reactome gene sets names to be tested. If not specified,
this function will use all the gene sets from the Reactome by default}

\item{mat_w}{This indicates whether to use the shared factor loadings 'W' (default TRUE)}

\item{mat_v}{This indicates which V matrix to be added to the analysis. It can be a numeric number or a list
of the numerics.}

\item{custom_gene_sets}{A named list of character vectors of entrez gene ids. If not specified,
this function will use all the gene symbols from the input matrix by default}
}
\value{
A list of matrices with GSEA analysis for each factor
}
\description{
Identify the biological pathways (gene sets from Reactome) that each metagene (factor) might belongs to.
}
\examples{
\dontrun{
Y <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 4, byrow = T)
Z <- matrix(c(1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2), nrow = 4, byrow = T)
X <- matrix(c(1, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6), nrow = 4, byrow = T)
ligerex <- createLiger(list(y_set = Y, z_set = Z, x_set = X))
ligerex <- normalize(ligerex)
# select genes
ligerex <- selectGenes(ligerex)
ligerex <- scaleNotCenter(ligerex)
ligerex <- optimizeALS(ligerex, k = 20)
ligerex <- quantil_norm(ligerex)
wilcox.results <- runGSEA(ligerex)
wilcox.results <- runGSEA(ligerex, mat_v = c(1, 2))
}

}