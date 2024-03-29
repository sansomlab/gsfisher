% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyseGenesets.R
\name{analyseGenesets}
\alias{analyseGenesets}
\title{Wrapped to run gene set enrichement on multiple annotation databases}
\usage{
analyseGenesets(
  foreground_ids,
  background_ids,
  gene_id_type = c("entrez", "ensembl"),
  kegg_pathways = NULL,
  gmt_files = c(),
  species = c("hs", "mm")
)
}
\arguments{
\item{foreground_ids}{A list of ENTREZ ids of interest
(e.g. significantly differentially expressed genes).}

\item{background_ids}{A list of background ENTREZ ids.
against which enrichment will be tested (i.e., the gene universe).}

\item{gene_id_type}{Either "entrez" (default) or "ensembl")}

\item{kegg_pathways}{A list of KEGG gene sets.}

\item{gmt_files}{A named list of paths to GMT files.}

\item{species}{Species identifier (only "hs" or "mm" are supported).}
}
\description{
A wrapper function to run Fisher tests for enrichement from GO categories,
KEGG pathways and GMT files
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
\code{\link{runGO}},
\code{\link{runKEGG}},
\code{\link{runGMT}}.
}
\author{
Steve Sansom
}
