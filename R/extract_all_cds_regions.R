get_cds_info = function(HGNC_symbol,toolDirectory, refGenome="hg19", annotationDb="ensemble"){
  #
  txPerGeneMethod=1
  txGenePairsFile     <- paste0(toolDirectory,"/references/GeneData/geneInfo/TREG_gene_transcript_pair.csv")
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
###
###
###
getGeneNames = function(refGenome="hg19", annotationDb="ensemble") {
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
###
###
###
extract_all_cds_regions = function(toolDirectory, refGenome="hg19", annotationDb="ensemble") {
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
  return(targetInfoGR)}
###
###
###