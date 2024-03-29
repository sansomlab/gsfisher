% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GMT.R
\name{runGMT}
\alias{runGMT}
\title{Run gene set enrichement on custom gene sets}
\usage{
runGMT(
  foreground_ids,
  background_ids,
  gene_id_type = c("entrez", "ensembl"),
  species = c("mm", "hs"),
  gmt_file,
  ...
)
}
\arguments{
\item{foreground_ids}{A list of ENTREZ ids of interest
(e.g. significantly differentially expressed genes).}

\item{background_ids}{A list of background ENTREZ ids.
against which enrichment will be tested (i.e., the gene universe).}

\item{gene_id_type}{Either "entrez" (default) or "ensembl".}

\item{species}{Species identifier (only "hs" or "mm" are supported).}

\item{gmt_file}{path to a GMT file.}
}
\description{
A wrapper function to run Fisher tests for enrichement from GMT files.
}
\examples{
gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
ann_gmt <- readGMT(gmtFile)
# Take 50\% of the first gene set as an example list of interest
foreground <- head(ann_gmt[[1]], length(ann_gmt[[1]]) / 2)
\dontrun{
# TODO
}
}
\seealso{
\code{\link{readGMT}},
\code{\link{runFisherTests}}.
}
\author{
Steve Sansom
}
