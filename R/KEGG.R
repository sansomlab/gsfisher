#' Run gene set enrichement on KEGG pathways
#'
#' A wrapper function to run Fisher tests for enrichement of KEGG Pathways.
#'
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes).
#' @param background_ids A list of background ENTREZ ids.
#' against which enrichment will be tested (i.e., the gene universe).
#' @param keggData A list of KEGG gene sets.
#' @param species Species identifier (only "hs" or "mm" are supported).
#' @param gene_id_type Either "entrez" (default) or "ensembl".
#'
#' @export
#'
#' @seealso
#' \code{\link{readGMT}},
#' \code{\link{fetchKEGG}},
#' \code{\link{runFisherTests}}.
#'
#' @author Steve Sansom
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)

#' # Take 50% of the first gene set as an example list of interest
#' foreground <- head(ann_gmt[[1]], length(ann_gmt[[1]]) / 2)

#' \dontrun{
#' # TODO
#' }
runKEGG <- function(
  foreground_ids, 
  background_ids, 
  keggData = NULL, 
  SYMBOL=NULL,
  gene_id_type=c("entrez","ensembl"),
  species=c("hs", "mm")
){
  species <- match.arg(species)
  
  gene_id_type <- match.arg(gene_id_type)
  
  foreground_ids <- getEntrez(foreground_ids, gene_id_type, species)
  background_ids <- getEntrez(background_ids, gene_id_type, species)
  
  if (is.null(keggData)) 
  { 
    message("getting KEGG pathways")
    keggData <- fetchKEGG(species=species) 
  }
  
  if(is.null(SYMBOL))
  {
    message("getting symbols")
    SYMBOL <- getSYMBOL(species) 
  }
  
  result_table <- runFisherTests(
    keggData$genesets, foreground_ids, background_ids, SYMBOL)
  
  result_table$description <- keggData$geneset_info[result_table$geneset_id, "description"]
  
  first_cols <- c("geneset_id", "description")
  other_cols <- colnames(result_table)[!colnames(result_table) %in% first_cols]
  
  result_table <- result_table[,c(first_cols, other_cols)]
  
  return(result_table)
}

#' Run KEGG analysis on a multi-sample results table.
#' @param results A multi-sample results table
#' @param species Either "mm" or "hs"
#' @param background_ids A vector of ENSEMBL gene ids. If NULL, taken from results.
#' @param sample_col The column in results that indicates the sample
#' @param gene_id_col The column containing the gene identifiers
#' @param gene_id_type Either "entrez" (default) or "ensembl".
#' @param keggData A list of KEGG gene sets.
#' (see \code{fetchKegg}).
#' @param p_col The column containing the p-values to use
#' @param p_threshold The significance threshold.
#' 
#' @export
runKEGG.all <- function(results=NULL,
                        species=c("mm","hs"),
                        background_ids=NULL, 
                        sample_col="cluster",
                        gene_id_col="gene_id",
                        gene_id_type=c("entrez","ensembl"),
                        keggData=NULL,
                        p_col="p_val_adj",
                        p_threshold=0.1)
{
  begin <- TRUE
  
  gene_id_type <- match.arg(gene_id_type)
  species <- match.arg(species)
  
  if(is.null(keggData)) {
  keggData <- fetchKEGG(species)
  }
  SYMBOL <- getSYMBOL(species)
  
  for(sample in unique(results[[sample_col]]))
  {
    message("working on sample:", sample)
    data <- results[results[[sample_col]]==sample,]
    
    if(is.null(background_ids))
    {
      background <- data[[gene_id_col]]
    } else {
      background <- background_ids
    }
    
    foreground <- data[[gene_id_col]][data[[p_col]] <= p_threshold]
    
    tmp <- runKEGG(foreground_ids = foreground,
                 background_ids = background,
                 keggData=keggData,
                 SYMBOL=SYMBOL,
                 gene_id_type = gene_id_type,
                 species = species)
    
    tmp[[sample_col]] <- sample 
    
    if(begin) { 
      res <- tmp
      begin <- FALSE
    } else {
      res <- rbind(res,tmp)
    }
  }
  res
}



