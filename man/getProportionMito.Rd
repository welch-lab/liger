% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/liger.R
\name{getProportionMito}
\alias{getProportionMito}
\title{Calculate proportion mitochondrial contribution}
\usage{
getProportionMito(object, use.norm = F)
}
\arguments{
\item{object}{\code{liger} object.}

\item{use.norm}{Whether to use cell normalized data in calculating contribution (default FALSE).}
}
\value{
Named vector containing proportion of mitochondrial contribution for each cell.
}
\description{
Calculates proportion of mitochondrial contribution based on raw or normalized data.
}
\examples{
\dontrun{
# liger object, factorization done
ligerex
ligerex@cell.data[["percent_mito"]] <- getProportionMito(ligerex)
}
}