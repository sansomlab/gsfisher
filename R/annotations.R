#' Fetch ENSEMBL ids, ENTREZ ids and gene names using bioMart
#'
#' Fetches annotations for human ("hs") or mouse ("mm") from the
#' Ensembl BioMart.
#'
#' @param ensembl_host The address of the ensembl host
#' The default \code{NULL} uses www.ensembl.org
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
#' @importFrom biomaRt useEnsembl getBM listAttributes
#'
#' @author Steve Sansom
#'
#' @examples
#' \dontrun{
#' ann_hs <- fetchAnnotation(species="hs")
#' ann_mm <- fetchAnnotation(species="mm")
#' }
fetchAnnotation <- function(species=c("hs", "mm"),
                            ensembl_version=NULL,
                            ensembl_host=NULL)
{

    species <- match.arg(species)

    if (species == "hs") {
        message("Using human biomart ...")
        dataset <- "hsapiens_gene_ensembl"
        name_col <- "external_gene_name"
    }
    if (species == "mm") {
        message("Using mouse biomart ...")
        dataset <- "mmusculus_gene_ensembl"
        name_col <- "mgi_symbol"
    }

    if(is.null(ensembl_host))
        {
            message("Using default host")
            ensembl <- useEnsembl(biomart="ensembl",
                                  dataset=dataset,
                                  version=ensembl_version)
        } else {
            message("Using host: ", ensembl_host)
            ensembl <- useEnsembl(biomart="ensembl",
                                  host=ensembl_host,
                                  dataset=dataset,
                                  version=ensembl_version)
        }

    # Deal with inconsistent biomart identifier names
    attrib_names <- listAttributes(ensembl)$name

    if("entrezgene" %in% attrib_names)
    {
      entrez_col = "entrezgene"
    } else if("entrezgene_id" %in% attrib_names)
    {
      entrez_col = "entrezgene_id"
    } else { stop("Entrez identifier attribute not found") }

    message("Entrez identified attribute name set to: ", entrez_col)

    message("Retrieving annotation")

    annotation <- getBM(
        attributes=c("ensembl_gene_id", entrez_col, name_col), mart=ensembl)

    message("Annotation retrieved")

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
#' (e.g., from \code{readGMT}).
#' @param ensembl_version Version of the ensembl annotation to use,
#' passed to \code{biomaRt::useEnsembl}.
#' The default \code{NULL} uses the current annotation release.
#'
#' @importFrom biomaRt useEnsembl getBM
#'
#' @export
#'
#' @seealso
#' \code{\link{readGMT}}
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

    # Deal with inconsistent biomart identifier names
    attrib_names <- listAttributes(humanEnsembl)$name

    if("entrezgene" %in% attrib_names)
    {
      entrez_col = "entrezgene"
    } else if("entrezgene_id" %in% attrib_names)
    {
      entrez_col = "entrezgene_id"
    } else { stop("Entrez identifier attribute not found") }

    message("Entrez identified attribute name set to: ", entrez_col)


    human2entrez <- getBM(attributes=c('ensembl_gene_id', entrez_col),
                          mart = humanEnsembl)

    mouse2entrez <- getBM(attributes=c('ensembl_gene_id', entrez_col),
                          mart = mouseEnsembl)

    mmGMT <- list()
    for(category  in names(GMT)) {
        human_entrez <- GMT[[category]]

        human_ensembl <- human2entrez$ensembl_gene_id[human2entrez[[entrez_col]] %in% GMT[[category]]]
        mouse_homolog <- human2mouse$mmusculus_homolog_ensembl_gene[human2mouse$ensembl_gene_id %in% human_ensembl]
        mouse_entrez <- unique(mouse2entrez[[entrez_col]][mouse2entrez$ensembl_gene_id %in% mouse_homolog])
        mmGMT[[category]] <- mouse_entrez[!is.na(mouse_entrez)]
    }

    return(mmGMT)
}

# Get the org.Xx.eg.db SYMBOL bimap object.
#' @param species Species, "mm" or "hs"
#'
#' @export
#'
#' @importFrom org.Mm.eg.db org.Mm.egSYMBOL
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL
getSYMBOL <- function(species=c("mm","hs"))
{
  if(species=="hs") {
    message("getting human gene symbols")
    require(org.Hs.eg.db)
    SYMBOL <- org.Hs.egSYMBOL
  } else if (species=="mm") {
    message("getting mouse gene symbols")
    require(org.Mm.eg.db)
    SYMBOL <- org.Mm.egSYMBOL
  } else { stop("species not recognised") }
  SYMBOL
}

#' Get the org.Xx.eg.db GO2ALLEGS bimap object.
#' @param species Species, "mm" or "hs"
#'
#'
#' @export
getGO <- function(species=c("mm","hs"))
{
  require(AnnotationDbi)
  species <- match.arg(species)
  if (species == "hs") {
    require(org.Hs.eg.db)
    genesets <- as.list(org.Hs.egGO2ALLEGS)
  } else if (species == "mm") {
    require(org.Mm.eg.db)
    genesets <- as.list(org.Mm.egGO2ALLEGS)
  } else { stop("species not recognised") }
  genesets
}

#' Get the org.Xx.eg.db ENSEMBL bimap object
#' @param species Species, "mm" or "hs"
#' @export
getENSEMBL <- function(species=c("mm","hs"))
{
  require(AnnotationDbi)
  species <- match.arg(species)
  if(species=="hs") {
  require(org.Hs.eg.db)
  ENSEMBL <- org.Hs.egENSEMBL
  } else if (species=="mm") {
  require(org.Mm.eg.db)
  ENSEMBL <- org.Mm.egENSEMBL
  } else { stop("species not recognised") }
  ENSEMBL
}

#' Translate ensembl gene ids to entrez gene ids
#' @param ensembl_ids A vector of ensembl gene_ids
#' @param ENSEMBL A org.Xx.eg.db ENSEMBL bimap object.
#' @param species Either "mm" or "hs".
#'
#' @export
ensembl2entrez <- function(ensembl_ids, ENSEMBL=NULL,species=c("mm","hs"))
{

  if(is.null(ENSEMBL)) {
    species <- match.arg(species)
    ENSEMBL <- getENSEMBL(species)
  }

  entrez <- na.omit(as.vector(unlist(AnnotationDbi::mget(
    ensembl_ids, revmap(ENSEMBL),ifnotfound = NA))))

  entrez
}

#' Verify input gene id type and return entrez identifiers
#' @param gene_ids A character vector of gene identifiers
#' @param gene_id_type Either "ensembl" or "entrez"
#' @param species Either "mm" or "hs"
getEntrez <- function(gene_ids,gene_id_type=c("ensembl","entrez"), species=c("mm","hs"))
{
  gene_id_type <- match.arg(gene_id_type)
  species <- match.arg(species)

  ENSEMBL <- getENSEMBL(species)

  if(gene_id_type == "ensembl")
  {
    if(!startsWith(gene_ids[1],"ENS")) {
      stop("Ensembl gene_ids specified but not supplied")
    }
    entrez_ids <- ensembl2entrez(gene_ids, ENSEMBL)
  } else if(gene_id_type =="entrez") {
    if(startsWith(gene_ids[1],"ENS")) {
      stop('"entrez" gene ids specified but "ensembl" gene ids supplied')
    }
    entrez_ids <- gene_ids
  } else { stop('gene_id_type must be either "ensembl" or "entrez"') }
  entrez_ids
}
