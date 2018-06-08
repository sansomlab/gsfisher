
#' Fetch KEGG annotations into list of gene sets
#'
#' @param species Species identifier (only "hs" or "mm" are supported).
#'
#' @return A named list of gene sets.
#'
#' @importFrom limma getGeneKEGGLinks
#'
#' @export
#'
#' @examples
#' \dontrun{
#' kegg_hs <- fetchKeggGeneSets(species="hs")
#' kegg_mm <- fetchKeggGeneSets(species="mm")
#' }
fetchKeggGeneSets <- function(species=c("hs", "mm")){
    species <- match.arg(species)

    if (species == "hs") {
        message("Using human KEGG pathways ...")
        kegg_species <- "hsa"
        description_suffix <- " - Homo sapiens \\(human\\)"
    }
    if (species == "mm") {
        message("Using mouse KEGG pathways ...")
        kegg_species <- "mmu"
        description_suffix=" - Mus musculus \\(mouse\\)"
    }

    # Get table of mapping between ENTREZ GeneID and KEGG PathwayID
    geneset_table <- getGeneKEGGLinks(species.KEGG=kegg_species)

    dupIdx <- duplicated(geneset_table)
    if (any(dupIdx)) {
        geneset_table <- geneset_table[!dupIdx, ]
        message(sum(dupIdx), " duplicated records removed.")
    }

    # Split table into named list of gene sets
    geneset_list <- tapply(geneset_table$GeneID, geneset_table$PathwayID, "c")

    return(geneset_list)
}

#' Make a KEGGCollection from a list of gene sets
#'
#' @param x A list of KEGG pathway gene sets
#' (e.g. from \code{fetchKeggGeneSets}).
#'
#' @return A \code{\link{KEGGCollection}}
#'
#' @importFrom GSEABase getGmt KEGGCollection EntrezIdentifier
#' @export
#'
#' @note
#' This function is hard-coded to use Entrez gene identifiers,
#' which is expected by \emph{gsfisher}.
#' It won't fail if that is not the case, but this may impact downstream
#' analyses.
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- read.gmt(gmtFile)
#'
#' kc <- makeKeggCollectionFromList(ann_gmt)
makeKeggCollectionFromList <- function(x) {
    # Write list to a temporary file, to benefit from the getGmt function
    tmpFile <- tempfile()
    writeGMT(x, tmpFile)
    # Use the friendly getGmt function to handle the hard work
    gsc <- getGmt(
        tmpFile,
        collectionType=KEGGCollection(), geneIdType=EntrezIdentifier()
        )

    return(gsc)
}
