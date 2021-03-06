% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genesetFisherTests.R
\name{fisherTest}
\alias{fisherTest}
\title{Perform a single Fisher test for gene set enrichement}
\usage{
fisherTest(n, genesets, foreground_ids, background_ids, symbols)
}
\arguments{
\item{n}{A numeric scalar indicating the index of the gene set to test.}

\item{genesets}{A named list of gene sets.}

\item{foreground_ids}{A list of ENTREZ ids of interest
(e.g. significantly differentially expressed genes).}

\item{background_ids}{A list of background ENTREZ ids.
against which enrichment will be tested (i.e., the gene universe).}

\item{symbols}{A named vector comprising entrez_id to gene_symbol mappings.}
}
\description{
Perform a single Fisher test for gene set enrichement
}
\examples{
# TODO
\dontrun{
# TODO
}
}
\seealso{
\code{\link{runFisherTests}}.
}
