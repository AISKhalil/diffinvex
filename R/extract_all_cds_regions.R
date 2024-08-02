#' Get CDS Information for a Specific Gene
#'
#' This function retrieves the coding sequence (CDS) information for a specified gene. It supports two types of annotation databases: Ensembl and UCSC. Depending on the selected database, it fetches the appropriate transcript and CDS information, potentially selecting the most-expressed or longest transcript if multiple options are available.
#'
#' @param HGNC_symbol A character string specifying the HGNC gene symbol for which CDS information is to be retrieved.
#' @param toolDirectory A character string specifying the path to the directory containing reference files and tools.
#' @param refGenome A character string specifying the reference genome (default is "hg19").
#' @param annotationDb A character string specifying the annotation database to use. Choices are "ensemble" and "UCSC" (default is "ensemble").
#'
#' @return A `GRanges` object containing the reduced CDS regions for the specified gene. The output is in UCSC coordinate style.
#'
#' @details
#' \itemize{
#'   \item For the Ensembl database:
#'     \item Fetches gene data based on the HGNC symbol.
#'     \item Selects the most-expressed transcript if available; otherwise, chooses the longest transcript.
#'     \item Retrieves the CDS for the selected transcript(s) and returns the reduced CDS regions.
#'   \item For the UCSC database:
#'     \item Fetches gene data using Entrez ID derived from the HGNC symbol.
#'     \item Retrieves all transcripts associated with the gene.
#'     \item Retrieves the CDS for the selected transcript(s) and returns the reduced CDS regions.
#' }
#'
#' @export
get_cds_info <- function(HGNC_symbol,toolDirectory, refGenome="hg19", annotationDb="ensemble"){
  #
  txPerGeneMethod=1
  txGenePairsFile     <- paste0(toolDirectory,"/data/GeneData/geneInfo/TREG_gene_transcript_pair.csv")
  txGenePairs         <- read.csv(txGenePairsFile)
  txGenePairs         <- data.frame(txGenePairs$Gene_name, sub("[.].*","",txGenePairs$Transcript_id))
  colnames(txGenePairs) <- c("Gene_name","Transcript_id")  
  #
  libs   <- import_annotation_db(refGenome, annotationDb)
  #
  hs     <- libs[[1]]
  genome <- libs[[2]]
  txdb   <- libs[[3]]
  #
  # 
  if(annotationDb == "ensemble") 
    {
        ##fetching gene data
        ## gene-data ##
        geneGR <- genes(hs, columns="gene_id", filter = GeneNameFilter(HGNC_symbol))
        #      
        ## check if the provided gene-name is present in the database or more than one ENSEMBL ID are presented
        if(is_empty(geneGR)){
            return(NULL)
        } else if(length(geneGR) > 1){
            geneGR <- geneGR[width(geneGR) == max(width(geneGR))]
        }
        #
        seqlevelsStyle(geneGR) <- "UCSC"
        ##
        ## transcripts ##
        geneID  <- geneGR$gene_id
        geneTxs <- transcripts(hs, columns="tx_id", filter = GeneIdFilter(geneID))
        #
        # choosing one-transcript per gene
        if(txPerGeneMethod == 1){
            # most-expressed transcript
            mExTx <- txGenePairs$Transcript_id[txGenePairs$Gene_name == HGNC_symbol]
            # longest transcript
            lnTx  <- geneTxs[width(geneTxs) == max(width(geneTxs))]
            #
            if((length(mExTx) == 1)&&(mExTx %in% geneTxs$tx_id)){
                geneTxs <- geneTxs[mExTx]
            }
            else{
                geneTxs <- lnTx
            }
        }
        ##
        ## cds ##
        geneCdsByTx <- cdsBy(hs, "tx", filter = TxNameFilter(geneTxs$tx_id))
        cdsReduced <- GenomicRanges::reduce(unlist(geneCdsByTx))
        seqlevelsStyle(cdsReduced) <- "UCSC"
        ##
        ##
        ## output
        return(cdsReduced)
    }
    ###############################
    ############################### 
    else if(annotationDb == "UCSC")
    {
      ##fetching gene data
      # gene-data
      hsSub <- ensembldb::select(hs,columns = c("SYMBOL", "ENSEMBL", "ENTREZID", "UCSCKG"), keys = keys(hs))
      geneInfo <- hsSub[hsSub$SYMBOL == HGNC_symbol,]
      geneNo <- geneInfo$ENTREZID[1]
      ##
      ##
      # transcript ##
      txs  <- transcriptsBy(txdb, by="gene")
      geneTxs <- txs[[as.character(geneNo)]]$tx_name
      ##
      ##
      ## cds ##
      geneCdsByTx <- cdsBy(txdb, "tx", use.names=TRUE)[geneTxs]
      cdsReduced  <- GenomicRanges::reduce(unlist(geneCdsByTx))
      ##
      ##
      return(cdsReduced)
    }
}



#' Retrieve Gene Names from an Annotation Database
#'
#' This function retrieves gene names from a specified annotation database. It supports two databases: Ensembl and UCSC. The function extracts gene symbols from the database and returns them as a vector of HGNC symbols.
#'
#' @param refGenome A character string specifying the reference genome (default is "hg19").
#' @param annotationDb A character string specifying the annotation database to use. Choices are "ensemble" and "UCSC" (default is "ensemble").
#'
#' @return A character vector containing gene names (HGNC symbols) retrieved from the specified annotation database.
#'
#' @details
#' \itemize{
#'   \item For the Ensembl database:
#'     \item Retrieves genes using the `genes` function and extracts gene names (HGNC symbols).
#'     \item Converts gene names to UCSC coordinate style.
#'   \item For the UCSC database:
#'     \item Retrieves genes and extracts unique gene identifiers.
#'     \item Uses these identifiers to query the Ensembl database for corresponding HGNC symbols.
#' }
#'
#' @export
getGeneNames <- function(refGenome="hg19", annotationDb="ensemble") {
  #
  libs <- import_annotation_db(refGenome, annotationDb)
  #
  hs     <- libs[[1]]
  genome <- libs[[2]]
  txdb   <- libs[[3]]
  #
  if(annotationDb == "ensemble") 
  {    
    genes <- genes(hs, columns=c("gene_name"))
    seqlevelsStyle(genes) <- "UCSC"
    #
    genes_HGNC_symbol <- genes$gene_name
  }
  else if(annotationDb == "UCSC")
  {
    genes <- genes(txdb, single.strand.genes.only=FALSE)
    #
    genesNos <- unique(data.frame(genes)$group_name)
    hsSub <- ensembldb::select(hs,columns = c("SYMBOL", "ENSEMBL"), keys = genesNos)
    genes_HGNC_symbol <- unique(hsSub$SYMBOL)
  }
  #
  return(genes_HGNC_symbol)
}



#' Extract All CDS Regions for Genes from Annotation Database
#'
#' This function retrieves coding sequence (CDS) regions for all genes listed in a specified annotation database. It generates and returns a consolidated GRanges object containing the CDS regions for the entire gene list.
#'
#' @param toolDirectory A character string specifying the path to the directory containing tools and reference files required for fetching CDS information.
#' @param refGenome A character string specifying the reference genome (default is "hg19").
#' @param annotationDb A character string specifying the annotation database to use. Choices are "ensemble" and "UCSC" (default is "ensemble").
#'
#' @return A GRanges object containing the sorted and reduced CDS regions for all genes from the specified annotation database.
#'
#' @details
#' \itemize{
#'   \item Retrieves a list of gene names using `getGeneNames`.
#'   \item For each gene, fetches CDS regions using `get_cds_info`.
#'   \item Consolidates all CDS regions into a single data frame.
#'   \item Converts the data frame to a GRanges object, sorts, and reduces it to remove redundant or overlapping regions.
#' }
#'
#' @export
extract_all_cds_regions <- function(toolDirectory, refGenome="hg19", annotationDb="ensemble") {
  #
  geneList <- getGeneNames(refGenome,annotationDb)
  #
  targetInfo <- lapply(geneList, function(geneID){
        print(geneID)
        targetGR <- get_cds_info(geneID, toolDirectory, refGenome, annotationDb)
        targetDF <- as.data.frame(targetGR)
        targetDF <- targetDF[,c("seqnames","start","end")]
        return(targetDF)})
  #
  targetInfoDF <- do.call(rbind, targetInfo)
  targetInfoGR <- GRanges(targetInfoDF)
  targetInfoGR <- GenomicRanges::sort(GenomicRanges::reduce(targetInfoGR))
  #
  return(targetInfoGR)
}