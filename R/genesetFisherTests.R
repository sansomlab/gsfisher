
#' Fetch ENSEMBL ids, ENTREZ ids and gene names using bioMart
#'
#' Fetches annotations for human ("hs") or mouse ("mm") from the
#' Ensembl BioMart.
#'
#' @param species Species identifiers (only "hs" or "mm" are supported).
#' @param ensembl_version Version of the ensembl annotation to use,
#' passed to \code{biomaRt::useEnsembl}.
#' The default \code{NULL} uses the current annotation release.
#'
#' @return A data.frame of three columns:
#' \code{c("ensembl_id", "entrez_id", "gene_name")}
#'
#' @details
#' Only human "hs" and mouse "mm" are supported.
#' Annotation must contain the columns \code{"entrez_id"} and
#' \code{"gene_name"}.
#'
#' @export
#'
#' @importFrom biomaRt useEnsembl getBM
#'
#' @examples
#' \dontrun{
#' ann_hs <- fetchAnnotations(species="hs")
#' ann_mm <- fetchAnnotations(species="mm")
#' }
fetchAnnotations <- function(species=c("hs", "mm"), ensembl_version=NULL) {

    species <- match.arg(species)

    if (species == "hs") {
        message("Using human biomart")
        dataset <- "hsapiens_gene_ensembl"
        namecol <- "external_gene_name"
    }
    if (species == "mm") {
        message("Using mouse biomart")
        dataset <- "mmusculus_gene_ensembl"
        namecol <- "mgi_symbol"
    }

    ensembl <- useEnsembl(biomart="ensembl", dataset=dataset, version=version)

    annotation <- getBM(
        attributes=c("ensembl_gene_id", "entrezgene", namecol), mart=ensembl)

    colnames(annotation) <- c("ensembl_id", "entrez_id", "gene_name")
    return(annotation)
}

#' Fetch KEGG pathway annotations
#'
#' @param species The species code, only "hs" or "mm" are supported.
#'
#' @return A list of two elements: "genesets" and "geneset_info".
#'
#' @importFrom limma getGeneKEGGLinks getKEGGPathwayNames
#'
#' @export
#'
#' @examples
#' \dontrun{
#' kegg_hs <- fetchKEGG(species="hs")
#' kegg_mm <- fetchKEGG(species="mm")
#' }
fetchKEGG <- function(species=c("hs", "mm")) {

    species <- match.arg(species)

    if (species == "hs") {
        kegg_species <- "hsa"
        description_suffix <- " - Homo sapiens \\(human\\)"
    }
    if (species == "mm") {
        kegg_species <- "mmu"
        description_suffix=" - Mus musculus \\(mouse\\)"
    }

    geneset_table <- getGeneKEGGLinks(species.KEGG=kegg_species)
    geneset_info <- getKEGGPathwayNames(species.KEGG=kegg_species)

    geneset_info$description <- gsub(description_suffix, "", geneset_info$Description)

    rownames(geneset_info) <- geneset_info$PathwayID

    genesets <- list()

    for (pathway_id in unique(geneset_table$PathwayID)) {
        genesets[[pathway_id]] <- as.character(
            geneset_table$GeneID[geneset_table$PathwayID == pathway_id])
    }

    result <- list(genesets=genesets, geneset_info=geneset_info)

    return(result)
}

#' Translate ENTREZ gene sets from human to mouse
#'
#' Uses ENSEMBL biomaRt for converting human MSigDB gene sets to mouse.
#'
#' @param GMT A named list of gene sets
#' (e.g., from \code{qusage::read.gmt}).
#' @param ensembl_version Version of the ensembl annotation to use,
#' passed to \code{biomaRt::useEnsembl}.
#' The default \code{NULL} uses the current annotation release.
#'
#' @importFrom biomaRt useEnsembl getBM
#'
#' @export
#'
#' @seealso
#' \code{\link[qusage]{read.gmt}}
#' \code{\link[biomaRt]{useEnsembl}}
#'
#' @examples
#' gmtFile <- system.file(package = "gsfisher", "extdata", "kegg_hs.gmt")
#' ann_gmt <- readGMT(gmtFile)
#' \dontrun{
#' mapENTREZhuman2mouse(ann_gmt, ensembl_version=NULL)
#' }
mapENTREZhuman2mouse <- function(GMT, ensembl_version=NULL) {
    # Catch before wasting time fetching from ENSEMBL
    if (missing(GMT)) {
        stop("GMT must be supplied")
    }

    humanEnsembl <- useEnsembl(biomart="ensembl",
                               dataset="hsapiens_gene_ensembl",
                               version=ensembl_version)
    mouseEnsembl <- useEnsembl(biomart="ensembl",
                               dataset="mmusculus_gene_ensembl",
                               version=ensembl_version)

    human2mouse <- getBM(
        attributes=c(
            'ensembl_gene_id',
            'mmusculus_homolog_ensembl_gene',
            'mmusculus_homolog_orthology_type'),
        mart = humanEnsembl)

    ## only map one:one orthologs.
    human2mouse <- human2mouse[human2mouse$mmusculus_homolog_orthology_type == "ortholog_one2one",]

    human2entrez <- getBM(attributes=c('ensembl_gene_id','entrezgene'),
                          mart = humanEnsembl)

    mouse2entrez <- getBM(attributes=c('ensembl_gene_id','entrezgene'),
                          mart = mouseEnsembl)

    mmGMT <- list()
    for(category  in names(GMT)) {
        human_entrez <- GMT[[category]]

        human_ensembl <- human2entrez$ensembl_gene_id[human2entrez$entrezgene %in% GMT[[category]]]
        mouse_homolog <- human2mouse$mmusculus_homolog_ensembl_gene[human2mouse$ensembl_gene_id %in% human_ensembl]
        mouse_entrez <- unique(mouse2entrez$entrezgene[mouse2entrez$ensembl_gene_id %in% mouse_homolog])
        mmGMT[[category]] <- mouse_entrez[!is.na(mouse_entrez)]
    }

    return(mmGMT)
}

#' Perform a single Fisher test for gene set enrichement
#'
#' @param n A numeric scalar indicating the index of the gene set to test.
#' @param genesets A named list of gene sets.
#' @param foreground_ids A list of ENTREZ ids of interest
#' (e.g. significantly differentially expressed genes)
#' @param background_ids A list of background ENTREZ ids
#' against which enrichment will be tested (i.e., the gene universe)
#' @param annotation A data.frame with columns "entrez_id" and "gene_name"
#' (see \code{fetchAnnotations})
#'
#' @importFrom stats fisher.test
#'
#' @seealso
#' \code{\link{fetchAnnotations}}
#'
#' @export
fisherTest <- function(
    n, genesets, foreground_ids, background_ids, annotation
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
    in_fg_names <- unique(annotation$gene_name[annotation$entrez_id %in% in_fg])

    ## get the intersection with the background set
    nbg <- length(intersect(set, background_ids))

    ## build the contingency table
    nnfg <- lfg - nfg
    nnbg <- lbg - nbg

    #################
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
        fg_freq, bg_freq,
        nfg, nbg, n_set,
        ft$p.value, ft$estimate,
        paste(in_fg, collapse=","),
        paste(in_fg_names, collapse=","))

    ## return results vector
    return(as.vector(result))
}

#' Run a set of Fisher tests for gene set enrichement
#'
#' @param foreground_ids the list of entrez ids of interest (e.g. significantly differentially expressed genes)
#' @param background_ids the list of entrez ids againt which enrichment will be tested (i.e. the gene universe)
#' @param named_geneset_list List of named gene sets.
#' @param annotation a dataframe with columns "entrez_id" and "gene_name" (see fetchAnnotations)
#'
#' @export
runFisherTests <- function(
    named_geneset_list,
    foreground_ids,
    background_ids,
    annotation
){
    ## annotation must contain columns
    ## entrez_id, gene_name
    result <- lapply(
        seq_along(names(named_geneset_list)),
        fisherTest,
        genesets=named_geneset_list,
        foreground_ids=foreground_ids,
        background_ids=background_ids,
        annotation=annotation)

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


#' A wrapper function to run Fisher tests for enrichement of Gene Ontology (GO) categories.
#' Depends on the bioconductor org.xx.eg.db and GO.db libraries
#'
#' @param foreground_ids the list of entrez ids of interest (e.g. significantly differentially expressed genes)
#' @param background_ids the list of entrez ids againt which enrichment will be tested (i.e. the gene universe)
#' @param annotation a dataframe with columns "entrez_id" and "gene_name" (see fetchAnnotations)
#' @param species The species code, only "hs" or "mm" are supported.
#'
#' @importFrom AnnotationDbi select
#' @importFrom GO.db GO.db
#' @importFrom org.Mm.eg.db org.Mm.egGO2ALLEGS
#' @importFrom org.Hs.eg.db org.Hs.egGO2ALLEGS
#'
#' @export
runGO <- function(
    foreground_ids, background_ids, annotation, species=c("hs","mm")
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
        genesets, foreground_ids, background_ids, annotation)

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


#' A wrapper function to run Fisher tests for enrichement of KEGG Pathways.
#'
#' @param foreground_ids The list of entrez ids of interest (e.g. significantly differentially expressed genes).
#' @param background_ids The list of entrez ids againt which enrichment will be tested (i.e. the gene universe).
#' @param annotation A dataframe with columns "entrez_id" and "gene_name" (see fetchAnnotations).
#' @param keggData List of KEGG gene sets.
#' @param species The species code, only "hs" or "mm" are supported.
#'
#' @export
runKEGG <- function(
    foreground_ids, background_ids, keggData = NULL, annotation, species=c("hs", "mm")
){
    species <- match.arg(species)

    if (is.null(keggData)) { keggData <- fetchKEGG(species=species) }

    result_table <- runFisherTests(
        keggData$genesets, foreground_ids, background_ids, annotation)

    result_table$description <- keggData$geneset_info[result_table$geneset_id, "description"]

    first_cols <- c("geneset_id", "description")
    other_cols <- colnames(result_table)[!colnames(result_table) %in% first_cols]

    result_table <- result_table[,c(first_cols, other_cols)]

    return(result_table)
}



#' A wrapper function to run Fisher tests for enrichement from GMT files
#'
#' @param foreground_ids The list of entrez ids of interest (e.g. significantly differentially expressed genes).
#' @param background_ids The list of entrez ids againt which enrichment will be tested (i.e. the gene universe).
#' @param gmt_file The location of the GMT file.
#' @param annotation A dataframe with columns "entrez_id" and "gene_name" (see fetchAnnotations).
#'
#' @export
runGMT <- function(foreground_ids, background_ids, gmt_file, annotation) {
    ## Get the gene sets
    gmt <- readGMT(gmt_file)

    ## Run the fisher tests
    result_table <- runFisherTests(gmt, foreground_ids, background_ids, annotation)

    return(result_table)
}

#' A wrapper function to run Fisher tests for enrichement from GO categories, KEGG pathways and GMT files
#'
#' @param foreground_ids The list of entrez ids of interest (e.g. significantly differentially expressed genes).
#' @param background_ids The list of entrez ids againt which enrichment will be tested (i.e. the gene universe).
#' @param gmt_files A named list of GMT file locations.
#' @param annotation A dataframe with columns "entrez_id" and "gene_name" (see fetchAnnotations).
#' @param kegg_pathways List of KEGG gene sets.
#' @param species The species code, only "hs" or "mm" are supported.
#'
#' @export
analyseGenesets <- function(
    foreground_ids, background_ids, annotation,
    kegg_pathways=NULL, gmt_files=c(), species="hs"
){
    results <- list()

    ## runGO
    message("Running GO ...")
    go_result <- runGO(foreground_ids, background_ids, annotation, species=species)

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
        foreground_ids, background_ids, kegg_pathways, annotation, species)

    message(paste0( "- nrow KEGG:", nrow(results[["KEGG"]]) ))

    ## runGMTs
    message("Running GMT files ...")

    for (geneset_name in names(gmt_files)) {
        results[[geneset_name]] <- runGMT(
            foreground_ids, background_ids,
            gmt_files[[geneset_name]],
            annotation)

        message(paste0( "- nrow ", geneset_name, ": ", nrow(results[[geneset_name]]) ))
    }

    return(results)
}
