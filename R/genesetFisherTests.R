#' Perform a single Fisher test for gene set enrichement
#'
#' @param n A numeric scalar indicating the index of the gene set to test.
#' @param genesets A named list of gene sets.
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes).
#' @param background_ids A list of background ENTREZ ids.
#' against which enrichment will be tested (i.e., the gene universe).
#' @param symbols A named vector comprising entrez_id to gene_symbol mappings.
#'
#' @importFrom stats fisher.test
#'
#' @seealso
#' \code{\link{runFisherTests}}.
#'
#' @export
#'
#' @examples
#' # TODO
#' \dontrun{
#' # TODO
#' }
fisherTest <- function(
    n, genesets, foreground_ids, background_ids, symbols){
  
    ## This is the inner loop...
  
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
    in_fg_names <- as.vector(symbols[in_fg])

    ## get the intersection with the background set
    nbg <- length(intersect(set, background_ids))
    
    ##if(nfg >= nbg) { 
    ##  stop("There must be more genes in the background than the foreground (2)") }
    
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
#' @export
#'
#' @seealso
#' \code{\link{readGMT}},
#' \code{\link{fisherTest}}.
#'
#' @author Steve Sansom
#'
#' @examples
#' # TODO

#' \dontrun{
#' # TODO
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
    
    nfg <- length(foreground_ids)
    nbg <- length(background_ids)
    
    ## Sanity check
    if (nfg == 0) { stop("No foreground genes") }
    if (nbg == 0) { stop("No background genes") }
    
    if (nfg >= nbg) { 
      stop("The number of foreground genes must be less than the number of background genes")
    }
      
    fg_percent = round(nfg/nbg*100,2)
    
    message("There are: ",length(foreground_ids)," unique foreground ids")
    message("There are: ",length(background_ids)," unique background ids")

    if(fg_percent < 25) { 
    message("The foreground ids represent ", fg_percent,"% of the gene universe")
    } else {
    message("!! The foreground ids represent ", fg_percent,"% of the gene universe!!",
            " Is the background correctly specified?")
    }
    
    symbols <- as.vector(unlist(AnnotationDbi::mget(foreground_ids, SYMBOL, ifnotfound = NA)))
    names(symbols) <- foreground_ids
    
    message("running the fisher tests")
    #  n, genesets, foreground_ids, background_ids, symbols
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

