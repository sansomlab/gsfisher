#' Perform a single Fisher test for gene set enrichement
#'
#' @param n A numeric scalar indicating the index of the gene set to test.
#' @param genesets A named list of gene sets.
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes).
#' @param background_ids A list of background ENTREZ ids.
#' against which enrichment will be tested (i.e., the gene universe).
#' @param annotations A data.frame with at least columns "entrez_id" and "gene_name"
#' (see \code{fetchAnnotation}).
#'
#' @importFrom stats fisher.test
#'
#' @seealso
#' \code{\link{fetchAnnotation}},
#' \code{\link{runFisherTests}}.
#'
#' @export
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)

#' # Take 50% of the first gene set as an example list of interest
#' index <- 1
#' foreground <- head(ann_gmt[[index]], length(ann_gmt[[index]]) / 2)

#' \dontrun{
#' # Fetch annotations
#' ann_hs <- fetchAnnotation(species="hs")
#' background <- subset(ann_hs, !is.na(entrez_id), "entrez_id", drop=TRUE)
#' # Demonstrate significant enrichment for the first gene set
#' result <- fisherTest(index, ann_gmt, foreground, background, ann_hs)
#' }
fisherTest <- function(
    n, genesets, foreground_ids, background_ids, symbols){
  
    set <- genesets[[n]]
    n_set <- length(set)

    ## get the set sizes
    lfg <- length(foreground_ids)
    lbg <- length(background_ids)

    ## get the intersection with the foreground set
    in_fg <- intersect(set, foreground_ids)
    nfg <- length(in_fg)

    ## if there is not intersection, skip
    if(nfg==0) { return(NA) }

    ## get the gene symbols of the gene set members in the foreground set
    ##in_fg_names <- unique(annotations$gene_name[annotations$entrez_id %in% in_fg])
    in_fg_names <- as.vector(symbols[in_fg])
    #  NA #na.omit(as.vector(unlist(
     # AnnotationDbi::mget(in_fg, SYMBOL, ifnotfound = NA))))

    ## get the intersection with the background set
    nbg <- length(intersect(set, background_ids))

    ## build the contingency table
    nnfg <- lfg - nfg
    nnbg <- lbg - nbg

    #
    # gene set is white ball, not gene set is black
    # and interprets the situation under the null hypothesis as taking a
    # sample of size a+c from an urn with a+b white balls and c+d black balls.
    #
    #                  | seen | not seen
    #  in gene set      |   a   | b
    #  not in gene set  |   c   | d

    ct <- matrix(c(
        nfg, nbg - nfg,
        nnfg, lbg - lfg - (nbg - nfg)
        ),
        ncol=2, byrow=TRUE)

    ## run one sided fishers' exact test
    ft <- fisher.test(ct, alternative="g")

    ## build the result
    result <- c(
        names(genesets)[n], # set_name
        nfg / lfg, # foreground frequency
        nbg / lbg, # background frequency
        nfg,
        nbg,
        n_set,
        ft$p.value,
        ft$estimate,
        paste(in_fg, collapse=","),
        paste(in_fg_names, collapse=","))

    ## return results vector
    result
}

#' Run a set of Fisher tests for gene set enrichement
#'
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes)
#' @param background_ids A list of background ENTREZ ids.
#' against which enrichment will be tested (i.e., the gene universe)
#' @param named_geneset_list List of named gene sets.
#' @param annotations A data.frame with at least columns "entrez_id" and "gene_name"
#' (see \code{fetchAnnotation}).
#'
#' @importFrom org.Mm.eg.db org.Mm.egSYMBOL
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL
#'
#' @export
#'
#' @seealso
#' \code{\link{readGMT}},
#' \code{\link{fetchAnnotation}},
#' \code{\link{fisherTest}}.
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
#' background <- subset(ann_hs, !is.na(entrez_id), "entrez_id", drop=TRUE)
#' result <- runFisherTests(ann_gmt, foreground, background, ann_hs)
#' }
runFisherTests <- function(
    named_geneset_list,
    foreground_ids,
    background_ids,
    SYMBOL
){
    
    ## ensure we are working with character vectors
    ## remove missing values
    ## ensure unique
    foreground_ids <- unique(as.character(foreground_ids[!is.na(foreground_ids)]))
    background_ids <- unique(as.character(background_ids[!is.na(background_ids)]))
  
    ## ensure background set contains the foreground_ids
    background_ids <- unique(c(foreground_ids, background_ids))
    
    symbols <- as.vector(unlist(AnnotationDbi::mget(foreground_ids, SYMBOL, ifnotfound = NA)))
    names(symbols) <- foreground_ids
    
    message("running the fisher tests")
    ## annotations must contain columns
    ## entrez_id, gene_name
    result <- lapply(
        seq_along(names(named_geneset_list)),
        "fisherTest",
        genesets=named_geneset_list,
        foreground_ids=foreground_ids,
        background_ids=background_ids,
        symbols=symbols)
    
    message("fisher tests complete")
    
    ## remove the empty results (test not performed)
    ## this may have implications for multiple testing correction
    result <- result[!sapply(
        result,
        function(x) { all(is.na(x)) }
        )]

    numeric_columns <- c(
        "fg_freq", "bg_freq", "n_fg", "n_bg", "n_set", "p.val", "odds.ratio")

    headings <- c(
        "geneset_id", numeric_columns, "entrez_ids", "gene_names")


    if (length(result) == 0) {
        result_table=data.frame(
            matrix(
                vector(), 0, length(headings),
                dimnames=list(c(), headings)),
            stringsAsFactors=FALSE)
    } else {
        ## build a results table
        result_table <- as.data.frame(
            matrix(unlist(result), ncol=length(result[[1]]), byrow=TRUE),
            stringsAsFactors=FALSE)

        names(result_table) <- headings

        ## reorder
        result_table <- result_table[
            , c("geneset_id", numeric_columns, "gene_names","entrez_ids")]

        ## make the numeric columns numeric..
        for (numcol in numeric_columns) {
            result_table[,numcol] <- as.numeric(result_table[, numcol])
        }

        result_table <- result_table[order(result_table$p.val), ]

    }

    ## Note that it is left to the user to adjust p-values as appropriate

    return(result_table)
}

