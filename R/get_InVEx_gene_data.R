get_InVEx_gene_data = function(HGNC_symbol,
                              genomic_tracks,
                              toolDirectory,                                      
                              refGenome="hg19",
                              neighbors_window=25000,
                              outputDirectory,
                              diffInVEx_mode=1,
                              diffInVEx_BW=50){  
  ################
  ## parameters ##
  ################
  minBackgroundWidth    <- diffInVEx_BW*1000
  maxDistanceFromTarget <- round(diffInVEx_BW*1000/2)
  minTargetWidth        <- 256
  matchingByCpGs        <- 1
  matchingByMappability <- 0
  #
  neighbors_window <- maxDistanceFromTarget #max(neighbors_window, maxDistanceFromTarget)
  ##
  ###################
  ## genomic files ##
  ###################
  annotationDb        <- "ensemble"
  txPerGeneMethod     <- 1 #1: most GTEx expressed transcript, 2: collapse all transcripts per gene, 3: test each transcript separatelly (to be implemented)
  txGenePairsFile     <- paste0(toolDirectory,"/references/GeneData/geneInfo/TREG_gene_transcript_pair.csv")
  txGenePairs         <- read.csv(txGenePairsFile)
  txGenePairs         <- data.frame(txGenePairs$Gene_name, sub("[.].*","",txGenePairs$Transcript_id))
  colnames(txGenePairs) <- c("Gene_name","Transcript_id")
  ##
  ## genomic tracks ##
  mappRegions         <- genomic_tracks[[1]]
  mappScores          <- genomic_tracks[[2]]
  filterRegions       <- genomic_tracks[[3]]
  filterCondition     <- genomic_tracks[[4]]
  filterInclusion     <- genomic_tracks[[5]]
  BlackListed_regions <- genomic_tracks[[6]]
  Microsatellites     <- genomic_tracks[[7]] 
  CTCF_regions  	  <- genomic_tracks[[8]]
  backgroundHotspots  <- genomic_tracks[[9]]
  uCpGs               <- genomic_tracks[[10]]
  promoterMotifs  	  <- genomic_tracks[[11]]
  shmRegions	      <- genomic_tracks[[12]] 
  allTargetRegions	  <- genomic_tracks[[13]] 
  ##
  ##################
  ## Output files ##
  ##################
  logFileName <- paste(outputDirectory, "/",HGNC_symbol, "_DiffInVEx.log", sep = "")
  logFile     <- file(logFileName, open="w")
  cat(HGNC_symbol, file = logFile, sep="\n", append = TRUE)
  ##
  ## gene information output files ##
  geneDB_Dir <- paste0(toolDirectory,"/references/GeneData/geneDB_",refGenome,"_",diffInVEx_BW,"Kb")
  #
  if (!dir.exists(geneDB_Dir)){
        dir.create(geneDB_Dir)
  }
  geneFile   <- paste0(geneDB_Dir,"/",HGNC_symbol,"_diffInVEx.RDS")
  ##
  ##
  ######################################
  ### Gene (exon&intron) information ###
  ######################################
  # In DiffInVEx: grExons are the coding regions of genes (cds), without the 5´utr and 3´útr
  #               grIntrons are the regions between coding regions (cds)
  #               grFlanks are the regions on both sides of coding regions (cds), including the 5´utr, 3´útr and promoters
  ###
  ###
  ###
  # getting grExons, grIntrons, and grFlanks
  geneData <- get_gene_info(HGNC_symbol, refGenome, annotationDb, txPerGeneMethod, txGenePairs)
  #
  if(is.null(geneData)){
    cat(".... does not exist in the annotation database \n")
    cat("Does not exist in the annotation database", file = logFile, sep="\n", append = TRUE)
    close(logFile)
    return(NULL)
  } else{
    cat(".....", file = logFile, sep="\n", append = TRUE)   
  }
  ####
  ####
  ####
  ####
  # Define exons, introns, and flanks
  # "important to create Granges without specified genome and no seqlengths, otherwise gaps will include extra flanks"
  grExons  <- GRanges(geneData)
  #
  # Exons: expanding each exon by 5bp as these bps have evidence for selection, 
  # except first and last exon (next to 5´UTR and 3´UTR)
  start(grExons) <- start(grExons) - 5
  end(grExons)   <- end(grExons) + 5
  start(grExons[which(start(grExons) == min(start(grExons)))]) <- start(grExons[which(start(grExons) == min(start(grExons)))]) + 5
  end(grExons[which(end(grExons) == max(end(grExons)))])       <- end(grExons[which(end(grExons) == max(end(grExons)))]) - 5
  grExons <- GenomicRanges::reduce(grExons)
  #
  # Introns
  grIntrons <- gaps(grExons)
  grIntrons <- grIntrons[-1] #to remove the flank from the chromosome start
  #
  # Flanks
  cdsStart  <- min(start(grExons))
  cdsEnd    <- max(end(grExons))
  cdsStrand <- as.character(strand(grExons)[1])
  cdsChr	<-  as.character(seqnames(grExons)[1])
  #
  chromLength <- (as.data.frame(getChromInfoFromUCSC(refGenome)) %>% dplyr::filter(chrom == cdsChr) %>% dplyr::select(size))$size
  distanceToTelo <- 100
  #
  grFlanksDF <- data.frame(seqnames = cdsChr, 
                           start  = c(max(cdsStart - neighbors_window,distanceToTelo), cdsEnd + 1), 
                           end    = c(cdsStart -1, min(cdsEnd + neighbors_window,chromLength - distanceToTelo)), 
                           strand = cdsStrand)
  grFlanks <- GRanges(grFlanksDF)
  #
  # nts GRanges with distance to target and context information
  grExons   <- getObjBps(grExons, refGenome, 0)
  grIntrons <- getObjBps(grIntrons, refGenome, 0)
  grFlanks  <- getObjBps(grFlanks, refGenome, 1)
  ####
  ####
  ####
  ####
  ## Filtering exons, introns, and flanks from gene's promoters and neighbor genes exons/promoters ##
  #
  # 1) Promoter NYTTCCG hypermutated motifs
  grExons    <- gr_setdiff(grExons, promoterMotifs)
  grIntrons  <- gr_setdiff(grIntrons, promoterMotifs)
  grFlanks   <- gr_setdiff(grFlanks, promoterMotifs)
  #
  #
  # 2) Excluding neighboring exons from flanking/intronic regions
  grIntrons <- gr_setdiff(grIntrons, allTargetRegions)
  grFlanks  <- gr_setdiff(grFlanks, allTargetRegions)        
  #
  cat(paste0("Exons width   = ", sum(width(grExons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Introns width = ", sum(width(grIntrons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Flanks width  = ", sum(width(grFlanks))), file = logFile, sep="\n", append = TRUE)
  #
  #
  # Filtering exons, introns, and flanks from low-quality regions: 
  # non-unique mapped regions, microsatellites, CTCF-enriched regions, and additional user-defined regions
  #
  # 3) Excluding regions that are not uniquely-mapped
  grExons   <- gr_intersect(grExons, mappRegions)
  grIntrons <- gr_intersect(grIntrons, mappRegions)
  grFlanks  <- gr_intersect(grFlanks, mappRegions)
  #
  cat("#", file = logFile, sep="\n", append = TRUE)
  cat(paste0("Exons width   (post mappability-filtering) = ", sum(width(grExons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Introns width (post mappability-filtering) = ", sum(width(grIntrons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Flanks width  (post mappability-filtering) = ", sum(width(grFlanks))), file = logFile, sep="\n", append = TRUE)
  #
  # 4) Excluding microsatellites longer than 6 nts
  grExons   = gr_setdiff(grExons, Microsatellites)  
  grIntrons = gr_setdiff(grIntrons, Microsatellites)
  grFlanks  = gr_setdiff(grFlanks, Microsatellites)
  #
  cat("#", file = logFile, sep="\n", append = TRUE)  
  cat(paste0("Exons width   (post microsatellites-filtering) = ", sum(width(grExons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Introns width (post microsatellites-filtering) = ", sum(width(grIntrons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Flanks width  (post microsatellites-filtering) = ", sum(width(grFlanks))), file = logFile, sep="\n", append = TRUE)
  #  
  # 5) intersect with CTCF hypermutated regions
  grExons   = gr_setdiff(grExons, CTCF_regions)
  grIntrons = gr_setdiff(grIntrons, CTCF_regions)  
  grFlanks  = gr_setdiff(grFlanks, CTCF_regions)
  #
  cat("#", file = logFile, sep="\n", append = TRUE)
  cat(paste0("Exons width   (post CTCF-filtering) = ", sum(width(grExons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Introns width (post CTCF-filtering) = ", sum(width(grIntrons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Flanks width  (post CTCF-filtering) = ", sum(width(grFlanks))), file = logFile, sep="\n", append = TRUE)    
  #
  # 6) excluding blacklisted regions
  grExons   = gr_setdiff(grExons, BlackListed_regions)
  grIntrons = gr_setdiff(grIntrons, BlackListed_regions)
  grFlanks  = gr_setdiff(grFlanks, BlackListed_regions)
  #
  cat("#", file = logFile, sep="\n", append = TRUE)  
  cat(paste0("Exons width   (post blackListed-filtering) = ", sum(width(grExons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Introns width (post BlackListed-filtering) = ", sum(width(grIntrons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Flanks width  (post BlackListed-filtering) = ", sum(width(grFlanks))), file = logFile, sep="\n", append = TRUE)      
  #
  # 7) filtering background hotspots
  grIntrons <- gr_setdiff(grIntrons,backgroundHotspots)
  grFlanks  <- gr_setdiff(grFlanks,backgroundHotspots)
  #
  cat("#", file = logFile, sep="\n", append = TRUE)  
  cat(paste0("Introns width (No Background hotspots) = ", sum(width(grIntrons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Flanks width  (No Background hotspots) = ", sum(width(grFlanks))), file = logFile, sep="\n", append = TRUE)      
  #
  # 8) filtering SHM regions
  grExons   = gr_setdiff(grExons, shmRegions)
  grIntrons = gr_setdiff(grIntrons, shmRegions)
  grFlanks  = gr_setdiff(grFlanks, shmRegions)
  #
  cat("#", file = logFile, sep="\n", append = TRUE)  
  cat(paste0("Exons width   (post SHM-filtering) = ", sum(width(grExons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Introns width (post SHM-filtering) = ", sum(width(grIntrons))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Flanks width  (post SHM-filtering) = ", sum(width(grFlanks))), file = logFile, sep="\n", append = TRUE)     
  #
  # 9) additional regions of interest (only for exonic regions)
  if (filterCondition) {
    if (filterInclusion) {
      grExons <- gr_intersect(grExons, filterRegions)
    } else {
      grExons <- gr_setdiff(grExons, filterRegions)
    }
  }
  ####
  ####
  # checking for exons with very few bps (post-filtering)
  if(sum(width(grExons)) < minTargetWidth){
      cat("#", file = logFile, sep="\n", append = TRUE)
      cat(paste0("... has post-filtering exon length = ", sum(width(grExons))), sep="\n")
      cat(paste0("... has post-filtering exon length < ", minTargetWidth), file = logFile, sep="\n", append = TRUE)
      close(logFile)
      return(NULL)} 
  ####
  ####
  ########################  
  # Target & background  #
  ########################
  #-Background selection-#
  grTarget     <- grExons
  grBackground <- select_background(grIntrons, grFlanks, minBackgroundWidth, maxDistanceFromTarget)
  #
  #-Adding CpGs information-#
  if(matchingByCpGs == 1){
      grCpGs <- add_CpGs_info(grTarget, grBackground, uCpGs)
      grTarget <- grCpGs[[1]]
      grBackground <- grCpGs[[2]]}  
  #
  #-Adding mappability-information per base-pair-#
  if(matchingByMappability == 1){
  		grMapp <- add_mappability_info(grTarget, grBackground, mappScores)
  		grTarget <- grMapp[[1]]
  		grBackground <- grMapp[[2]]}
  #
  #
  # matching target and background nts
  grBackground_NoMatching <- GenomicRanges::reduce(GenomicRanges::sort(grBackground))
  cat(sprintf("tri-matching\n"))
  grBackground_tri        <- match_target_background(grTarget, grBackground, 1, matchingByMappability)
  cat(sprintf("penta-matching\n"))
  grBackground_penta      <- match_target_background(grTarget, grBackground, 0, matchingByMappability)
  #
  grDict <- list()
  grDict[["target"]]                <- GenomicRanges::reduce(GenomicRanges::sort(grTarget))
  grDict[["background_NoMatching"]] <- grBackground_NoMatching
  grDict[["background_tri"]]        <- GenomicRanges::reduce(grBackground_tri)
  grDict[["background_penta"]]      <- GenomicRanges::reduce(grBackground_penta)
  #
  cat("#", file = logFile, sep="\n", append = TRUE)
  cat(paste0("Target width = ", sum(width(grTarget))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Background width (no-matching) = ", sum(width(grBackground_NoMatching))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Background width (tri-matching) = ", sum(width(grBackground_tri))), file = logFile, sep="\n", append = TRUE)
  cat(paste0("Background width (penta-matching) = ", sum(width(grBackground_penta))), file = logFile, sep="\n", append = TRUE)
  close(logFile)  
  ###
  ###
  ###  
  if(diffInVEx_mode == 3 || diffInVEx_mode == 4){
            saveRDS(grDict, geneFile)}
  #
  return(grDict)
}
###
###
###
load_InVEx_gene_data = function(HGNC_symbol,
                                toolDirectory,                                      
                                refGenome="hg19",
                                diffInVEx_BW=50){  
  ###
  geneDB_Dir <- paste0(toolDirectory,"/references/GeneData/geneDB_",refGenome,"_",diffInVEx_BW,"Kb")
  geneFile   <- paste0(geneDB_Dir,"/",HGNC_symbol,"_diffInVEx.RDS")
  ###
  if(file.exists(geneFile)){ 
    grDict  <- readRDS(geneFile)
    return(grDict)
  }
  else{
    cat(paste0(HGNC_symbol," files do not exist, please create them first \n"))
    return(NULL)
  }
  ###
}
