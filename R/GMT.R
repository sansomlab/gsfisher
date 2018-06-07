
#' Read gene sets from a GMT file
#'
#' @param filepath The location of the GMT file.
#'
#' @details
#' Deprecated. Use \code{qusage::read.gmt} instead.
#'
#' @export
#'
#' @seealso
#' \code{\link[qusage]{read.gmt}}
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
#' @param geneset A list of gene sets  (e.g. from \code{qusage::read.gmt}).
#' @param outfile Path to the output file.
#'
#' @seealso
#' \code{\link{read.gmt}}
#'
#' @export
#'
#' @author Steve Sansom
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- read.gmt(gmtFile)
#'
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
