#' Run gene set enrichement on Gene Ontology categories
#'
#' A wrapper function to run Fisher's test for enrichement
#' on Gene Ontology (GO) categories.
#' Depends on the bioconductor \code{org.xx.eg.db} and \code{GO.db} packages.
#'
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes).
#' @param background_ids A list of background ENTREZ ids.
#' against which enrichment will be tested (i.e., the gene universe).
#' @param annotations A data.frame with at least columns "entrez_id" and "gene_name"
#' (see \code{fetchAnnotation}).
#' @param species Species identifier (only "hs" or "mm" are supported).
#' @param gene_id_type Either "entrez" (default") or "ensembl".
#'
#' @export
#'
#' @seealso
#' \code{\link{runFisherTests}}.
#'
#' @author Steve Sansom
#'
#' @examples
#' # TODO
#
#' \dontrun{
#' # TODO
#' }
runGO <- function(
  foreground_ids, 
  background_ids, 
  genesets=NULL, 
  SYMBOL=NULL, 
  gene_id_type=c("entrez","ensembl"), 
  species=c("hs","mm")
){
  require(GO.db)
  species <- match.arg(species)

  gene_id_type <- match.arg(gene_id_type)
  
  foreground_ids <- getEntrez(foreground_ids, gene_id_type, species)
  background_ids <- getEntrez(background_ids, gene_id_type, species)
  
  ## Get the gene sets
  if(is.null(genesets))
  {
    message("getting genesets")
    genesets <- getGO(species)
  } 
  
  if(is.null(SYMBOL))
  {
    message("getting symbols")
    SYMBOL <- getSYMBOL(species) 
  }
  
  ## run the fisher tests
  result_table <- runFisherTests(
    genesets, foreground_ids, background_ids, SYMBOL)
  
  ## retrieve additional info about the GO categories
  ## and add to the results table
  
  message("annotating the results")
  
  go_info <- AnnotationDbi::select(
    GO.db, keys=result_table$geneset_id,
    columns=c("TERM", "ONTOLOGY"), keytype="GOID")
  
  ## ensure congruence
  rownames(go_info) <- go_info$GOID
  
  result_table$description <- go_info[result_table$geneset_id, "TERM"]
  result_table$ontology <- go_info[result_table$geneset_id, "ONTOLOGY"]
  
  first_cols <- c("geneset_id", "description", "ontology")
  other_cols <- colnames(result_table)[!colnames(result_table) %in% first_cols]
  
  result_table <- result_table[,c(first_cols, other_cols)]
  
  return(result_table)
}

#' Run GO analysis on a multi-sample results table.
#' @param results A multi-sample results table
#' @param species Either "mm" or "hs"
#' @param background_ids A vector of gene ids. If NULL, taken from results.
#' @param gene_id_col The name of the column containing the gene identifiers.
#' @param gene_id_type Either "entrez" (default) or "ensembl".
#' @param sample_col The column in results that indicates the sample
#' @param p_col The column containing the p-values to use
#' @param p_threshold The significance threshold.
#' 
#' @export
runGO.all <- function(results=NULL,
                      species=c("mm","hs"),
                      background_ids=NULL,
                      gene_id_col="gene_id",
                      gene_id_type=c("entrez","ensembl"),
                      sample_col="cluster", 
                      p_col="p_val_adj",
                      p_threshold=0.1)
{
  begin <- TRUE
  
  species <- match.arg(species)
  gene_id_type <- match.arg(gene_id_type)
  
  genesets <- getGO(species)
  SYMBOL <- getSYMBOL(species)
  
  for(sample in unique(results[[sample_col]]))
  {
    message("working on ", sample_col,": ", sample)
    data <- results[results[[sample_col]]==sample,]

    if(is.null(background_ids))
    {
      background <-data[[gene_id_col]]
    } else {
      background <- background_ids
    }
    
    foreground <- data[[gene_id_col]][data[[p_col]] <= p_threshold]
    
    tmp <- runGO(foreground_ids = foreground,
                 background_ids = background,
                 genesets=genesets,
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