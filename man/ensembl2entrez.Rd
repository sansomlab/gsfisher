% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotations.R
\name{ensembl2entrez}
\alias{ensembl2entrez}
\title{Translate ensembl gene ids to entrez gene ids}
\usage{
ensembl2entrez(ensembl_ids, ENSEMBL = NULL, species = c("mm", "hs"))
}
\arguments{
\item{ensembl_ids}{A vector of ensembl gene_ids}

\item{ENSEMBL}{A org.Xx.eg.db ENSEMBL bimap object.}

\item{species}{Either "mm" or "hs".}
}
\description{
Translate ensembl gene ids to entrez gene ids
}
