#' Cluster a set of genesets by gene membership
#' 
#' @param results_table The table of geneset results
#' @param genes_col The name of the column containing comma-separated genenames
#' @param dist_method The name of the distance method to pass to "dist"
#' @param hclust_method The name of the clustering method to pass to "hclust"
#' 
#' @export
#' 
clusterGenesetsByGenes <- function(results_table,
                                   genes_col = "gene_names",
                                   dist_method = "manhattan",
                                   hclust_method = "ward.D2")
{
  
  # get the per-geneset gene lists
  gene_lists <- sapply(results_table[[genes_col]], function(x) strsplit(x, ","))
  names(gene_lists) <- results_table$description
  
  all_genes <- unique(unlist(gene_lists))
  
  # compute the gene by geneset membership matrix
  membership_lists <- lapply(gene_lists, function(x) as.numeric(all_genes %in% x))
  
  membership_matrix <- t(matrix(unlist(membership_lists), 
                                ncol = length(all_genes), 
                                byrow = TRUE))
  
  colnames(membership_matrix) <- names(membership_lists)
  rownames(membership_matrix) <- all_genes
  
  # cluster by gene membership
  geneset_dist <- dist(t(membership_matrix), method="manhattan")
  geneset_clust <- hclust(geneset_dist, method="ward.D2")

  geneset_clust
}
