
#' Fetch ENSEMBL ids, ENTREZ ids and gene names using bioMart
#'
#' Fetches annotations for human ("hs") or mouse ("mm") from the
#' Ensembl BioMart.
#'
#' @param species Species identifier (only "hs" or "mm" are supported).
#' @param ensembl_version Version of the ensembl annotation to use,
#' passed to \code{biomaRt::useEnsembl}.
#' The default \code{NULL} uses the current annotation release.
#'
#' @return A data.frame of three columns:
#' \code{c("ensembl_id", "entrez_id", "gene_name")}
#'
#' @details
#' Only human "hs" and mouse "mm" are supported.
#' Annotation must contain the columns \code{"entrez_id"} and
#' \code{"gene_name"}.
#'
#' @export
#'
#' @importFrom biomaRt useEnsembl getBM
#'
#' @author Steve Sansom
#'
#' @examples
#' \dontrun{
#' ann_hs <- fetchAnnotations(species="hs")
#' ann_mm <- fetchAnnotations(species="mm")
#' }
fetchAnnotations <- function(species=c("hs", "mm"), ensembl_version=NULL) {

    species <- match.arg(species)

    if (species == "hs") {
        message("Using human biomart ...")
        dataset <- "hsapiens_gene_ensembl"
        namecol <- "external_gene_name"
    }
    if (species == "mm") {
        message("Using mouse biomart ...")
        dataset <- "mmusculus_gene_ensembl"
        namecol <- "mgi_symbol"
    }

    ensembl <- useEnsembl(biomart="ensembl", dataset=dataset, version=ensembl_version)

    annotation <- getBM(
        attributes=c("ensembl_gene_id", "entrezgene", namecol), mart=ensembl)

    colnames(annotation) <- c("ensembl_id", "entrez_id", "gene_name")
    return(annotation)
}

#' Fetch KEGG pathway annotations
#'
#' @param species Species identifier (only "hs" or "mm" are supported).
#'
#' @return A list of two elements: "genesets" and "geneset_info".
#'
#' @importFrom limma getGeneKEGGLinks getKEGGPathwayNames
#'
#' @export
#'
#' @author Steve Sansom
#'
#' @examples
#' \dontrun{
#' kegg_hs <- fetchKEGG(species="hs")
#' kegg_mm <- fetchKEGG(species="mm")
#' }
fetchKEGG <- function(species=c("hs", "mm")) {

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

    geneset_table <- getGeneKEGGLinks(species.KEGG=kegg_species)
    geneset_info <- getKEGGPathwayNames(species.KEGG=kegg_species)

    geneset_info$description <- gsub(description_suffix, "", geneset_info$Description)

    rownames(geneset_info) <- geneset_info$PathwayID

    genesets <- list()

    for (pathway_id in unique(geneset_table$PathwayID)) {
        genesets[[pathway_id]] <- as.character(
            geneset_table$GeneID[geneset_table$PathwayID == pathway_id])
    }

    result <- list(genesets=genesets, geneset_info=geneset_info)

    return(result)
}

#' Translate ENTREZ gene sets from human to mouse
#'
#' Uses ENSEMBL biomaRt for converting human MSigDB gene sets to mouse.
#'
#' @param GMT A named list of gene sets
#' (e.g., from \code{qusage::read.gmt}).
#' @param ensembl_version Version of the ensembl annotation to use,
#' passed to \code{biomaRt::useEnsembl}.
#' The default \code{NULL} uses the current annotation release.
#'
#' @importFrom biomaRt useEnsembl getBM
#'
#' @export
#'
#' @seealso
#' \code{\link{read.gmt}}
#' \code{\link[biomaRt]{useEnsembl}}
#'
#' @author Steve Sansom
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)
#' \dontrun{
#' mapENTREZhuman2mouse(ann_gmt, ensembl_version=NULL)
#' }
mapENTREZhuman2mouse <- function(GMT, ensembl_version=NULL) {
    # Catch before wasting time fetching from ENSEMBL
    if (missing(GMT)) {
        stop("GMT must be supplied")
    }

    humanEnsembl <- useEnsembl(biomart="ensembl",
                               dataset="hsapiens_gene_ensembl",
                               version=ensembl_version)
    mouseEnsembl <- useEnsembl(biomart="ensembl",
                               dataset="mmusculus_gene_ensembl",
                               version=ensembl_version)

    human2mouse <- getBM(
        attributes=c(
            'ensembl_gene_id',
            'mmusculus_homolog_ensembl_gene',
            'mmusculus_homolog_orthology_type'),
        mart = humanEnsembl)

    ## only map one:one orthologs.
    human2mouse <- human2mouse[human2mouse$mmusculus_homolog_orthology_type == "ortholog_one2one",]

    human2entrez <- getBM(attributes=c('ensembl_gene_id','entrezgene'),
                          mart = humanEnsembl)

    mouse2entrez <- getBM(attributes=c('ensembl_gene_id','entrezgene'),
                          mart = mouseEnsembl)

    mmGMT <- list()
    for(category  in names(GMT)) {
        human_entrez <- GMT[[category]]

        human_ensembl <- human2entrez$ensembl_gene_id[human2entrez$entrezgene %in% GMT[[category]]]
        mouse_homolog <- human2mouse$mmusculus_homolog_ensembl_gene[human2mouse$ensembl_gene_id %in% human_ensembl]
        mouse_entrez <- unique(mouse2entrez$entrezgene[mouse2entrez$ensembl_gene_id %in% mouse_homolog])
        mmGMT[[category]] <- mouse_entrez[!is.na(mouse_entrez)]
    }

    return(mmGMT)
}
