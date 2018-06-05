#' Draw a heatmap of enriched genesets in all samples
#' @param results_table A gofisher multi-sample results table 
#' @param max_rows The maximum number of rows to show in the heatmap
#' @param adjust_pvalues Adjust p-values across all clusters
#' @param pvalue_threshold P-value threshold
#' @param maxl  The maximum length for geneset names
#' @param show_common Show genesets enriched in all samples.
#' @param sample_id_col The column of the results_table containing the sample id
sampleEnrichmentHeatmap <- function(results_table, 
                                   max_rows=50,
                                   min_genes=2,
                                   adjust_pvalues=FALSE,
                                   padjust_method="BH",
                                   pvalue_threshold=0.1,
                                   maxl=45,
                                   show_common=TRUE,
                                   sample_id_col="cluster")
{
  
  total_n_clust <- length(unique(results_table$cluster))
  
  
  ## compute FDR accross all clusters (before filtering!)
  results_table$p.adj <- p.adjust(contents$p.val, method=padjust_method)
  
  ## Filter the genesets
  results_table <- results_table[results_table$n_fg >= min_genes,]
  
  if(adjust_pvalues)
  {
    results_table <- results_table[results_table$p.adj < pvalue_threshold,]
  } else {
    results_table <- results_table[results_table$p.val < pvalue_threshold,]
  }
  
  if(!show_common) 
  {
    results_table <- results_table[results_table$n_clust < total_n_clust,]
  }
  
  ## Check for genesets after filtering
  if(nrow(results_table)==0)
  {
    stop("No genesets pass filters")
  }
    
  ## Compute the number of clusters in which the geneset is enriched
  id_tab <- table(results_table$geneset_id)
  results_table$n_clust <- id_tab[results_table$geneset_id]
  
  ## Sort by p value
  results_table <- results_table[order(results_table$cluster,results_table$p.val),]
  
  # set the maximum number of output rows
  ntake <- round(max_rows/total_n_clust)
  
  # Identify the genesets to show in the heatmap
  hmap_genesets <- c()
  
  for(clust in unique(as.character(results_table[[sample_id_col]])))
  {
  
    temp <- results_table[results_table[[sample_id_col]]==clust,]
    nrows <- nrow(temp)
    if(nrows==0) {
      print(paste0("no enriched genesets found for cluster: ",clust))
      next }
  
    temp <- temp[1:min(nrows,ntake),]
  
    hmap_genesets <- c(hmap_genesets,temp$geneset_id)
  }

  if(!"description" %in% colnames(results_table)) 
  { 
    results_table$description <- results_table$geneset_id 
  }
  
  take <- c("geneset_id", "description", "odds.ratio", "cluster")
  temp <- results_table[results_table$geneset_id %in% unique(hmap_genesets),take]
  xx <- temp$description
  xx[is.na(xx)] <- "n/a"
  xx[nchar(xx)>maxl] <- paste0(strtrim(xx[nchar(xx)>maxl],maxl),"'")
  temp$description <- xx
  
  lu <- unique(temp[,c("geneset_id","description")])
  lu$description <- make.unique(lu$description)
  rownames(lu) <- lu$geneset_id
  
  dd <- dcast(temp, geneset_id~cluster, value.var="odds.ratio")
  rownames(dd) <- lu[dd$geneset_id,"description"]
  dd$geneset_id <- NULL
  
  ## deal with missing clusters
  for(clust in as.character(c(first:last)))
  {
    if(!clust %in% colnames(dd))
    {
      print("adding column")
      print(clust)
      dd[[clust]] <- 0
    }
  }
  
  ## convert to matrix and deal with missing values
  dd <- as.matrix(dd)
  dd[is.na(dd)] <- 0
  dd[is.infinite(dd)] <- max(dd[!is.infinite(dd)])
  
  m <- dd
  
  ramp_colors <- c("blue","darkblue","black","yellow","red")
  
  nbreaks=50 # number of graduations
  rm <- range(m)
  rm <- c(-2.5,2.5)
  breaks=seq(rm[1],rm[2],diff(rm)/nbreaks) 
  colors=colorRampPalette(ramp_colors)(nbreaks)
  
  if(nrow(m)==1)
    {
      print("Warning, matrix only has one row. It will be duplicated so that a heatmap can be made.")
      m <- rbind(m,m)
    }
    
    ## enforce column order
    mcols <- colnames(m)
    mcols <- mcols[order(mcols)]
    m <- m[,mcols]
    
    par(cex.main=0.7)
    heatmap.2(m,
              col=colors,
              breaks=breaks,
              scale="row",
              Colv=F,
              mar=c(4,30),
              trace="none",
              key.title = "",
              key.par=list(mar=c(4.4,0.3,4.4,0.3)),
              density.info=c("none"),
              lwid = c(1,8),
              lhei = c(1,4),
              key.xlab = "row-scaled\nodds-ratio",
              key.ylab = "",
              main=paste(geneset,"pathways",sep=" "),
              xlab=opt[[sample_id_col]]type,
              cexRow = 0.85,
              cexCol = 1.3
    )
}
  