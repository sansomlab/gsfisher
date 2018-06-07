
#' Perform a single Fisher test for gene set enrichement
#'
#' @param n A numeric scalar indicating the index of the gene set to test.
#' @param genesets A named list of gene sets.
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes).
#' @param background_ids A list of background ENTREZ ids.
#' against which enrichment will be tested (i.e., the gene universe).
#' @param annotations A data.frame with at least columns "entrez_id" and "gene_name"
#' (see \code{fetchAnnotations}).
#'
#' @importFrom stats fisher.test
#'
#' @seealso
#' \code{\link{fetchAnnotations}}
#'
#' @export
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)
#'
#' # Take 50% of the first gene set as an example list of interest
#' index <- 1
#' foreground <- head(ann_gmt[[index]], length(ann_gmt[[index]]) / 2)
#'
#' \dontrun{
#' # Fetch annotations
#' ann_hs <- fetchAnnotations(species="hs")
#' background <- subset(ann_hs, !is.na(entrez_id), "entrez_id", drop=TRUE)
#' # Demonstrate significant enrichment for the first gene set
#' result <- fisherTest(index, ann_gmt, foreground, background, ann_hs)
#' }
fisherTest <- function(
    n, genesets, foreground_ids, background_ids, annotations
){

    ## ensure we are working with character vectors
    ## remove missing values
    ## ensure unique
    foreground_ids <- unique(as.character(foreground_ids[!is.na(foreground_ids)]))
    background_ids <- unique(as.character(background_ids[!is.na(background_ids)]))

    ## ensure background set contains the foreground_ids
    background_ids <- unique(c(foreground_ids, background_ids))

    set_name <- names(genesets)[n]
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
    in_fg_names <- unique(annotations$gene_name[annotations$entrez_id %in% in_fg])

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

    ## get frequencies
    fg_freq <- nfg / lfg
    bg_freq <- nbg / lbg

    ## build the result
    result <- c(
        set_name,
        fg_freq,
        bg_freq,
        nfg,
        nbg,
        n_set,
        ft$p.value,
        ft$estimate,
        paste(in_fg, collapse=","),
        paste(in_fg_names, collapse=","))

    ## return results vector
    return(as.vector(result))
}

#' Run a set of Fisher tests for gene set enrichement
#'
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes)
#' @param background_ids A list of background ENTREZ ids.
#' against which enrichment will be tested (i.e., the gene universe)
#' @param named_geneset_list List of named gene sets.
#' @param annotations A data.frame with at least columns "entrez_id" and "gene_name"
#' (see \code{fetchAnnotations}).
#'
#' @export
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)
#'
#' # Take 50% of the first gene set as an example list of interest
#' foreground <- head(ann_gmt[[1]], length(ann_gmt[[1]]) / 2)
#'
#' \dontrun{
#' # Fetch annotations
#' ann_hs <- fetchAnnotations(species="hs")
#' background <- subset(ann_hs, !is.na(entrez_id), "entrez_id", drop=TRUE)
#' result <- runFisherTests(ann_gmt, foreground, background, ann_hs)
#' }
runFisherTests <- function(
    named_geneset_list,
    foreground_ids,
    background_ids,
    annotations
){
    ## annotations must contain columns
    ## entrez_id, gene_name
    result <- lapply(
        seq_along(names(named_geneset_list)),
        fisherTest,
        genesets=named_geneset_list,
        foreground_ids=foreground_ids,
        background_ids=background_ids,
        annotations=annotations)

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
#' (see \code{fetchAnnotations}).
#' @param species Species identifier (only "hs" or "mm" are supported).
#'
#' @importFrom AnnotationDbi select
#' @importFrom GO.db GO.db
#' @importFrom org.Mm.eg.db org.Mm.egGO2ALLEGS
#' @importFrom org.Hs.eg.db org.Hs.egGO2ALLEGS
#'
#' @export
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)
#'
#' # Take 50% of the first gene set as an example list of interest
#' foreground <- head(ann_gmt[[1]], length(ann_gmt[[1]]) / 2)
#'
#' \dontrun{
#' # Fetch annotations
#' ann_hs <- fetchAnnotations(species="hs")
#' background <- subset(ann_hs, !is.na(entrez_id), "entrez_id", drop=TRUE)
#' result <- runGO(foreground, background, ann_hs, "hs")
#' }
runGO <- function(
    foreground_ids, background_ids, annotations, species=c("hs","mm")
){
    species <- match.arg(species)

    ## Sanity check
    if (length(foreground_ids) == 0) { stop("No foreground genes") }

    ## Get the gene sets
    if (species == "hs") {
        genesets <- as.list(org.Hs.egGO2ALLEGS)
    }
    if (species == "mm") {
        genesets <- as.list(org.Mm.egGO2ALLEGS)
    }

    ## run the fisher tests
    result_table <- runFisherTests(
        genesets, foreground_ids, background_ids, annotations)

    ## retrieve additional info about the GO categories
    ## and add to the results table

    go_info <- select(
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


#' Run gene set enrichement on KEGG pathways
#'
#' A wrapper function to run Fisher tests for enrichement of KEGG Pathways
#'
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes).
#' @param background_ids A list of background ENTREZ ids.
#' against which enrichment will be tested (i.e., the gene universe).
#' @param annotations A data.frame with at least columns "entrez_id" and "gene_name"
#' (see \code{fetchAnnotations}).
#' @param keggData A list of KEGG gene sets.
#' @param species Species identifier (only "hs" or "mm" are supported).
#'
#' @export
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)
#'
#' # Take 50% of the first gene set as an example list of interest
#' foreground <- head(ann_gmt[[1]], length(ann_gmt[[1]]) / 2)
#'
#' \dontrun{
#' # Fetch annotations
#' ann_hs <- fetchAnnotations(species="hs")
#' kegg_hs <- fetchKEGG(species="hs")
#' background <- subset(ann_hs, !is.na(entrez_id), "entrez_id", drop=TRUE)
#' result <- runKEGG(foreground, background, kegg_hs, ann_hs, "hs")
#' }
runKEGG <- function(
    foreground_ids, background_ids, keggData = NULL, annotations,
    species=c("hs", "mm")
){
    species <- match.arg(species)

    if (is.null(keggData)) { keggData <- fetchKEGG(species=species) }

    result_table <- runFisherTests(
        keggData$genesets, foreground_ids, background_ids, annotations)

    result_table$description <- keggData$geneset_info[result_table$geneset_id, "description"]

    first_cols <- c("geneset_id", "description")
    other_cols <- colnames(result_table)[!colnames(result_table) %in% first_cols]

    result_table <- result_table[,c(first_cols, other_cols)]

    return(result_table)
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
#' @param annotations A data.frame with at least columns "entrez_id" and "gene_name"
#' (see \code{fetchAnnotations}).
#'
#' @export
#'
#' @author Steve Sansom and Kevin-Rue-Albrecht
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)
#'
#' # Take 50% of the first gene set as an example list of interest
#' foreground <- head(ann_gmt[[1]], length(ann_gmt[[1]]) / 2)
#'
#' \dontrun{
#' # Fetch annotations
#' ann_hs <- fetchAnnotations(species="hs")
#' background <- subset(ann_hs, !is.na(entrez_id), "entrez_id", drop=TRUE)
#' result <- runGMT(foreground, background, gmtFile, ann_hs)
#' }
runGMT <- function(foreground_ids, background_ids, gmt_file, annotations) {
    ## Get the gene sets
    gmt <- readGMT(gmt_file)

    ## Run the fisher tests
    result_table <- runFisherTests(gmt, foreground_ids, background_ids, annotations)

    return(result_table)
}

#' A wrapper function to run Fisher tests for enrichement from GO categories, KEGG pathways and GMT files
#'
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes).
#' @param background_ids A list of background ENTREZ ids.
#' against which enrichment will be tested (i.e., the gene universe).
#' @param gmt_files A named list of paths to GMT files.
#' @param annotations A data.frame with at least columns "entrez_id" and "gene_name"
#' (see \code{fetchAnnotations}).
#' @param kegg_pathways A list of KEGG gene sets.
#' @param species Species identifier (only "hs" or "mm" are supported).
#'
#' @export
#'
#' @author Steve Sansom and Kevin-Rue-Albrecht
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)
#'
#' # Take 50% of the first gene set as an example list of interest
#' foreground <- head(ann_gmt[[1]], length(ann_gmt[[1]]) / 2)
#'
#' \dontrun{
#' # Fetch annotations
#' ann_hs <- fetchAnnotations(species="hs")
#' kegg_hs <- fetchKEGG(species="hs")
#' background <- subset(ann_hs, !is.na(entrez_id), "entrez_id", drop=TRUE)
#' result <- analyseGenesets(
#'     foreground, background, ann_hs,
#'     kegg_hs, gmt_files=c(extdata=gmtFile), species="hs")
#' }
analyseGenesets <- function(
    foreground_ids, background_ids, annotations,
    kegg_pathways=NULL, gmt_files=c(), species=c("hs", "mm")
){
    species <- match.arg(species)

    if (length(gmt_files)) {
        if (is.null(names(gmt_files))) {
            print(gmt_files)
            stop("names(gmt_files) cannot be NULL")
        }
        if (any(names(gmt_files) == "")) {
            print(gmt_files)
            stop("All elements of gmt_files must be named")
        }
    }

    results <- list()

    ## runGO
    message("Running GO ...")
    go_result <- runGO(foreground_ids, background_ids, annotations, species=species)

    message(paste0( "- nrow GO: ", nrow(go_result) ))

    ## make separate tables for the different GO ontology types
    for (ontology in unique(go_result$ontology)) {
        temp <- go_result[go_result$ontology == ontology, ]
        temp$ontology <- NULL
        results[[paste("GO", ontology, sep=".")]] <- temp
    }

    ## runKEGG
    message("Running KEGG ...")
    results[["KEGG"]] <- runKEGG(
        foreground_ids, background_ids, kegg_pathways, annotations, species)

    message(paste0( "- nrow KEGG:", nrow(results[["KEGG"]]) ))

    ## runGMTs
    message("Running GMT files ...")

    for (geneset_name in names(gmt_files)) {
        message("- Running ", geneset_name, " ...")
        results[[geneset_name]] <- runGMT(
            foreground_ids, background_ids,
            gmt_files[[geneset_name]],
            annotations)

        message(paste0( "- nrow ", geneset_name, ": ", nrow(results[[geneset_name]]) ))
    }

    return(results)
}
