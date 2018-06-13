
#' Read gene sets from a GMT file
#'
#' @param filepath The location of the GMT file.
#' @details
#'
#' @export
#'
#' @author Steve Sansom
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)
readGMT <- function(filepath) {

    conn = file(filepath, "r")
    genesets <- list()

    while ( TRUE ) {
        line = readLines(conn, n=1)
        if ( length(line) == 0 ) {
            break
        } else {
            contents <- strsplit(line, "\t")[[1]]
            genesets[[contents[1]]] <- contents[3:length(contents)]
        }
    }

    close(conn)

    return(genesets)
}

#' Write gene sets to a GMT file
#'
#' @param geneset A list of gene sets  (e.g. from \code{readGMT}).
#' @param outfile Path to the output file.
#'
#' @seealso
#' \code{\link{readGMT}}
#'
#' @export
#'
#' @author Steve Sansom
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)


#' outfile <- tempfile()
#' writeGMT(ann_gmt, outfile)
writeGMT <- function(geneset, outfile) {
    conn <- file(outfile)

    lines <- c()
    for (category in names(geneset)) {
        contents <- c(category, "gsFisher.writeGMT", geneset[[category]])
        line <- paste(contents, collapse="\t")
        lines <- c(lines,line)
    }

    writeLines(lines, conn)
    status <- close(conn)

    return(status)
}

#' Run gene set enrichement on custom gene sets
#'
#' A wrapper function to run Fisher tests for enrichement from GMT files.
#'
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes).
#' @param background_ids A list of background ENTREZ ids.
#' against which enrichment will be tested (i.e., the gene universe).
#' @param gmt_file path to a GMT file.
#' @param gene_id_type Either "entrez" (default) or "ensembl".
#'
#'
#' @export
#'
#' @seealso
#' \code{\link{readGMT}},
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
runGMT <- function(
    foreground_ids,
    background_ids,
    gmt_file,
    species=c("mm","hs"),
    gene_id_type=c("entrez","ensembl"),
    SYMBOL=NULL)
{
    ## Get the gene sets
    gmt <- readGMT(gmt_file)
    species <- match.arg(species)

    gene_id_type <- match.arg(gene_id_type)

    foreground_ids <- getEntrez(foreground_ids, gene_id_type, species)
    background_ids <- getEntrez(background_ids, gene_id_type, species)

    if(is.null(SYMBOL))
    {
        message("getting symbols")
        SYMBOL <- getSYMBOL(species)
    }

    ## Run the fisher tests
    result_table <- runFisherTests(gmt, foreground_ids, background_ids, SYMBOL)

    return(result_table)
}


#' Run GMT analysis on a multi-sample results table.
#' @param results A multi-sample results table
#' @param species Either "mm" or "hs"
#' @param background_ids A vector of ENSEMBL gene ids. If NULL, taken from results.
#' @param sample_col The column in results that indicates the sample
#' @param gmt_file A GMT file.
#' @param gene_id_col The column containing the gene identifiers
#' @param gene_id_type Either "entrez" (default) or "ensembl".
#' @param p_col The column containing the p-values to use
#' @param p_threshold The significance threshold.
#'
#' @export
runGMT.all <- function(results=NULL,
                       species=c("mm","hs"),
                       background_ids=NULL,
                       sample_col="cluster",
                       gene_id_col="gene_id",
                       gene_id_type=c("entrez","ensembl"),
                       gmt_file=NULL,
                       p_col="p_val_adj",
                       p_threshold=0.1)
{
    begin <- TRUE

    species <- match.arg(species)
    gene_id_type <- match.arg(gene_id_type)

    if(is.null(gmt_file)) {
        stop("GMT file must be specified")
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

        tmp <- runGMT(foreground_ids = foreground,
                      background_ids = background,
                      gmt_file=gmt_file,
                      gene_id_type=gene_id_type,
                      SYMBOL=SYMBOL,
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


