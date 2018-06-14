#' Utility function to clean up desriptions for plots
#' 
#' @param descriptions Vector of descriptions to clean
#' @param remove Vector of string patterns to delete 
#' @param maxl  Max length of the descriptions
#' 
#' @export
formatDescriptions <- function(xx,
                               remove=c(),
                               maxl=45)
{
  for(rmstr in remove)
  {
  xx <- gsub(rmstr, "", xx)
  }
  xx <- gsub("_", " ", xx)
  
  # covert all upper case strings to all lower case. 
  xx <- sapply(xx, function(x) if(x == toupper(x)) { tolower(x) } else { x }) 

  xx[is.na(xx)] <- "n/a"
  xx[nchar(xx)>maxl] <- paste0(strtrim(xx[nchar(xx)>maxl], maxl), "'")
  xx
}

#' Draw a heatmap of enriched genesets in all samples
#'
#' @param results_table A \emph{gofisher} multi-sample results table.
#' @param max_rows Maximum number of rows to show in the heatmap.
#' @param adjust_pvalues Logical value indicating whether to adjust p-values
#' across all sample.
#' @param pvalue_threshold P-value threshold.
#' @param min_odds_ratio Only genesets with odds ratios equal to or greater will be considered.
#' @param maxl Maximum length for gene set names (trimmed otherwise).
#' @param show_common Logical value indicating whether to show gene sets
#' enriched in all samples.
#' @param min_genes Minimum number of genes to filter gene sets.
#' @param padjust_method Correction method. Can be abbreviated
#' (see \code{\link{p.adjust.methods}}).
#' @param title Title of the heat map plot.
#' @param sample_id_col The column of the results_table containing the sample id
#' @param sample_ids A vector of sample identifiers to report. Can be used to subset the frame, or to add missing columns.
#'
#' @importFrom gplots heatmap.2
#' @importFrom reshape2 dcast
#' @importFrom grDevices colorRampPalette
#' @importFrom stats p.adjust
#' @importFrom graphics par
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
#' @seealso
#' \code{\link{readGMT}}
#' \code{\link{p.adjust}}
#'
#' @author Steve Sansom
sampleEnrichmentHeatmap <- function(
                                    results_table,
                                    max_rows=50,
                                    min_genes=2,
                                    min_odds_ratio=1.5,
                                    p_col="p.val",
                                    adjust_pvalues=FALSE,
                                    padjust_method="BH",
                                    pvalue_threshold=0.1,
                                    maxl=45,
                                    show_common=TRUE,
                                    sample_id_col="cluster",
                                    sample_ids=NULL,
                                    title="Enriched gene sets"
                                    )
{
  
  results_table[[sample_id_col]] <- as.character(results_table[[sample_id_col]])

  # handle the sample identifiers
  if(!is.null(sample_ids))
  {
    sample_ids <- as.character(sample_ids)
    results_table <- results_table[results_table[[sample_id_col]] %in% sample_ids,] 
  } else {
    sample_ids <- unique(results_table[[sample_id_col]])
  }
   
  total_n_sample <- length(sample_ids)

    ## (Adjust) and filter on p-values before other filters are applied
    if (adjust_pvalues) {
        ## compute FDR accross all samples
        message("filtering on adjusted p values")
        results_table$p.adj <- p.adjust(results_table[[p_col]], method=padjust_method)
        results_table <- results_table[
            results_table$p.adj < pvalue_threshold & !is.na(results_table$p.adj), ]
    } else {
        message("filtering on unadjusted p values")
        results_table <- results_table[
            results_table[[p_col]] < pvalue_threshold & !is.na(results_table[[p_col]]), ]
    }

    if (!show_common) {
        results_table <- results_table[results_table$n_sample < total_n_sample, ]
    }
    
    ## Filter based on number of genes in the foreground
    results_table <- results_table[
      results_table$n_fg >= min_genes & !is.na(results_table$n_fg),]
    
    ## Filter based on odds ratio
    results_table <- results_table[
      results_table$odds.ratio >= min_odds_ratio & !is.na(results_table$odds.ratio),]

    ## Check for gene sets after filtering
    if (nrow(results_table) == 0) {
        stop("No gene set passed filters")
    }

    ## Compute the number of samples in which each gene set is enriched
    id_tab <- table(results_table$geneset_id)
    results_table$n_sample <- id_tab[results_table$geneset_id]

    ## Sort by p value
    results_table <- results_table[order(results_table[[sample_id_col]],results_table[[p_col]]),]

    # set the maximum number of output rows
    ntake <- round(max_rows / total_n_sample)

    # Identify the genesets to show in the heatmap
    hmap_genesets <- c()

    for ( sample in unique(results_table[[sample_id_col]])) {

        temp <- results_table[results_table[[sample_id_col]] == sample, ]
        nrows <- nrow(temp)
        if (nrows==0) {
            print(paste0("no enriched genesets found for", sample_id_col, " : ", sample))
            next
        }

        temp <- temp[1:min(nrows,ntake), ]

        hmap_genesets <- c(hmap_genesets, temp$geneset_id)
    }

    if (!"description" %in% colnames(results_table)) {
        results_table$description <- results_table$geneset_id
    }

    take <- c("geneset_id", "description", "odds.ratio", sample_id_col)
    temp <- results_table[results_table$geneset_id %in% unique(hmap_genesets), take]

    # process the term description
    xx <- temp$description
    xx <- formatDescriptions(xx, c("REACTOME_", "BIOCARTA_"), maxl)
    temp$description <- xx
    
    lu <- unique(temp[,c("geneset_id", "description")])
    lu$description <- make.unique(lu$description)
    rownames(lu) <- lu$geneset_id

    # unmelt the data
    dd <- dcast(temp, geneset_id ~ get(sample_id_col), value.var="odds.ratio")
    rownames(dd) <- lu[dd$geneset_id, "description"]
    dd$geneset_id <- NULL

    ## deal with missing samples
    for (sample in sample_ids) {
        if (!sample %in% colnames(dd)) {
            message("adding column: ", sample)
            dd[[sample]] <- NA
        }
    }
    
    ## convert to matrix and deal with infinite values
    dd <- as.matrix(dd)
    # dd[is.na(dd)] <- 0
    dd[is.infinite(dd)] <- 256

    m <- dd

    ## heatmap.2 requires a matrix with a least 2 rows and 2 columns.
    if (nrow(m) == 1 ) {
        warning(
            "Matrix only has one row. ",
            "It will be duplicated so that a heatmap can be made.")
        m <- rbind(m, m)
    }

    ## this should almost never happen.
    if (ncol(m) == 1) {
        warning(
            "Matrix only has one column. ",
            "It will be duplicated so that a heatmap can be made.")
        m <- cbind(m,m)
        scale="none"
    }
  
    ## enforce column order
    cnames <- colnames(m)
    typed_names <- type.convert(cnames, as.is=T)
    cnames <- cnames[order(typed_names)]
    m <- m[, cnames]
    
    ## compute the row dendrogram
    mtmp <- m
    mtmp[is.na(mtmp)] <- 0

    mdist <- dist(mtmp)
    mclust <- hclust(mdist)
    mrowDen <- as.dendrogram(mclust)

    # log2 transform the odds ratios
    mtrans <- log2(mtmp+1)

    # return the NA values to the matrix
    mtrans[is.na(m)] <- NA

    nbreaks <- 50 # number of graduations
    rm <- range(m)

    # specify the color range as 
    # minimum = minimum non-NA value
    # maximum = equivalent to an odds ratio of 256.
    rm <- c(min(log2(min_odds_ratio + 1)), 8)
    breaks=seq(rm[1], rm[2], diff(rm) / nbreaks)
    
    colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))

    if(adjust_pvalues)
    {
        ylab = paste0("over-represented\ngenesets (",
                      padjust_method," adjusted p < ",
                      pvalue_threshold, ")")
    } else {
        ylab = paste0("over-represented\ngenesets (",
                      "p < ",pvalue_threshold, ")")
    }

    op <- par() # backup current session settings
    par(cex.main=0.7)#, family="narrow")
    hm <- heatmap.2(
        mtrans,
        col=colors,
        breaks=breaks,
        scale="none",
        Rowv=mrowDen,
        Colv=F,
        margins=c(4, 22),
        trace="none",
        key.title = "",
        #keysize =0.4,
        colsep=c(0:(ncol(mtrans)+1)),
        rowsep=c(0:(nrow(mtrans)+1)),
        sepcolor="grey",
        sepwidth=c(0.01,0.01),
        key.par=list(mar=c(3.7, 0.3,0.6, 0.3)),
        density.info=c("none"),
        lwid = c(1, 8),
        lhei = c(1, 6),
        key.xlab = "odds-ratio\nlog2(n+1)",
        key.ylab = "",
        main=title,
        ylab=ylab,
        xlab=sample_id_col,
        cexRow = 0.85,
        cexCol = 1.3,
        na.color="lightgrey"
    )
    par(op) # restore original session settings
    return(hm)
}
