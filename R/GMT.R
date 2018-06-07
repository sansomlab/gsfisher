
#' Read genesets from a GMT file
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



#' Write genesets to a GMT file
#'
#' @param geneset A geneset list (e.g. from readGMT)
#' @param outfile The name of the outfile
#'
#' @export
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
