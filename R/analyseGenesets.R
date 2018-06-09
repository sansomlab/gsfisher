#' Wrapped to run gene set enrichement on multiple annotation databases
#'
#' A wrapper function to run Fisher tests for enrichement from GO categories,
#' KEGG pathways and GMT files
#'
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes).
#' @param background_ids A list of background ENTREZ ids.
#' against which enrichment will be tested (i.e., the gene universe).
#' @param gmt_files A named list of paths to GMT files.
#' @param annotations A data.frame with at least columns "entrez_id" and "gene_name"
#' (see \code{fetchAnnotation}).
#' @param kegg_pathways A list of KEGG gene sets.
#' @param species Species identifier (only "hs" or "mm" are supported).
#'
#' @export
#'
#' @seealso
#' \code{\link{runGO}},
#' \code{\link{runKEGG}},
#' \code{\link{runGMT}}.
#'
#' @author Steve Sansom
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)

#' # Take 50% of the first gene set as an example list of interest
#' foreground <- head(ann_gmt[[1]], length(ann_gmt[[1]]) / 2)

#' \dontrun{
#' # Fetch annotations
#' ann_hs <- fetchAnnotation(species="hs")
#' kegg_hs <- fetchKEGG(species="hs")
#' background <- subset(ann_hs, !is.na(entrez_id), "entrez_id", drop=TRUE)
#' result <- analyseGenesets(
#'     foreground, background, ann_hs,
#'     kegg_hs, gmt_files=c(extdata=gmtFile), species="hs")
#' }
analyseGenesets <- function(
  foreground_ids, background_ids, 
  kegg_pathways=NULL, gmt_files=c(), species=c("hs", "mm")
){
  species <- match.arg(species)
  
  if (length(gmt_files)) {
    if (is.null(names(gmt_files))) {
      print(gmt_files)
      stop("names(gmt_files) cannot be NULL")
    }
    if (any(names(gmt_files) == "")) {
      print(gmt_files)
      stop("All elements of gmt_files must be named")
    }
  }
  
  SYMBOL <- getSYMBOL(species)
  
  results <- list()
  
  ## runGO
  message("Running GO ...")
  go_result <- runGO(foreground_ids, background_ids, genesets=NULL, SYMBOL, species=species)
  
  message(paste0( "- nrow GO: ", nrow(go_result) ))
  
  ## make separate tables for the different GO ontology types
  for (ontology in unique(go_result$ontology)) {
    temp <- go_result[go_result$ontology == ontology, ]
    temp$ontology <- NULL
    results[[paste("GO", ontology, sep=".")]] <- temp
  }
  
  ## runKEGG
  message("Running KEGG ...")
  results[["KEGG"]] <- runKEGG(
    foreground_ids, background_ids, kegg_pathways, SYMBOL, species)
  
  message(paste0( "- nrow KEGG:", nrow(results[["KEGG"]]) ))
  
  ## runGMTs
  message("Running GMT files ...")
  
  for (geneset_name in names(gmt_files)) {
    message("- Running ", geneset_name, " ...")
    results[[geneset_name]] <- runGMT(
      foreground_ids, background_ids,
      gmt_files[[geneset_name]],
      species,
      SYMBOL)
    
    message(paste0( "- nrow ", geneset_name, ": ", nrow(results[[geneset_name]]) ))
  }
  
  return(results)
}