#' Filter a genesets results table
#' 
#' @param results_table A \emph{gofisher} multi-sample results table.
#' @param adjust_pvalues Logical value indicating whether to adjust p-values
#' across all sample.
#' @param pvalue_threshold P-value threshold.
#' @param min_odds_ratio Only genesets with odds ratios equal to or greater will be considered.
#' @param min_foreground_genes Minimum number of over-represent genes in foreground.
#' @param max_genes_geneset Maximum number of genes in the full geneset (used to remove generic categories).
#' @param padjust_method Correction method. Can be abbreviated
#' (see \code{\link{p.adjust.methods}}).
#'
#' @importFrom stats p.adjust
#'
#' @export
#'
#' @seealso
#'
#' @author Steve Sansom

filterGenesets <- function(results_table,
                           min_foreground_genes=3,
                           max_genes_geneset=500,
                           min_odds_ratio=2,
                           p_col="p.val",
                           adjust_pvalues=FALSE,
                           padjust_method="BH",
                           pvalue_threshold=0.05
)
{
  ## (Adjust) and filter on p-values before other filters are applied
  if (adjust_pvalues) {
    ## compute FDR accross all samples
    results_table$p.adj <- p.adjust(results_table[[p_col]], method=padjust_method)
    
    message("filtering on adjusted p values")
    results_table <- results_table[
      results_table$p.adj < pvalue_threshold & !is.na(results_table$p.adj), ]
  } else {
    message("filtering on unadjusted p values")
    results_table <- results_table[
      results_table[[p_col]] < pvalue_threshold & !is.na(results_table[[p_col]]), ]
  }
  
  ## Filter based on number of genes in the foreground
  results_table <- results_table[
    results_table$n_fg >= min_foreground_genes & !is.na(results_table$n_fg),]
  
  ## Filter based on number of genes in the geneset
  results_table <- results_table[
    results_table$n_set <= max_genes_geneset,]
  
  ## Filter based on odds ratio
  results_table <- results_table[
    results_table$odds.ratio >= min_odds_ratio & !is.na(results_table$odds.ratio),]
  
  ## Check for gene sets after filtering
  if (nrow(results_table) == 0) {
    stop("No gene set passed filters")
  }
  
  results_table
}