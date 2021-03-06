% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotations.R
\name{fetchKEGG}
\alias{fetchKEGG}
\title{Fetch KEGG pathway annotations}
\usage{
fetchKEGG(species = c("hs", "mm"))
}
\arguments{
\item{species}{Species identifier (only "hs" or "mm" are supported).}
}
\value{
A list of two elements: "genesets" and "geneset_info".
}
\description{
Fetch KEGG pathway annotations
}
\examples{
\dontrun{
kegg_hs <- fetchKEGG(species="hs")
kegg_mm <- fetchKEGG(species="mm")
}
}
\author{
Steve Sansom
}
