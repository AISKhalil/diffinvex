###-------- loading data and libraries --------###
packages <- c("BSgenome","ensembldb","Biostrings")
invisible(quiet(lapply(packages, library, character.only = TRUE), all = T))
###
###
###
import_annotation_db = function(refGenome="hg19", annotationDb="ensemble"){
  if(annotationDb == "ensemble") 
    {
      ##reference genome + annotation  
      if(refGenome == "hg19")
        {
        suppressMessages(library(EnsDb.Hsapiens.v75))
        suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
        hs <- EnsDb.Hsapiens.v75
        hs <- addFilter(hs, SeqNameFilter(c(seq(1,22,by = 1),"X","Y"))) #to return only transcripts of interest
        hs <- addFilter(hs, TxBiotypeFilter("protein_coding"))        
        genome <- BSgenome.Hsapiens.UCSC.hg19
        }
      else if(refGeneome == "hg38") 
        {
        suppressMessages(library(EnsDb.Hsapiens.v86))
        suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
        hs <- EnsDb.Hsapiens.v86
        hs <- addFilter(hs, SeqNameFilter(c(seq(1,22,by = 1),"X","Y"))) #to return only transcripts of interest
        hs <- addFilter(hs, TxBiotypeFilter("protein_coding"))
        genome <- BSgenome.Hsapiens.UCSC.hg38
        }
        txdb <- NULL
    }
  else if(annotationDb == "UCSC")
    {
      ##reference genome + annotation  
      if(refGenome == "hg19")
        {
        suppressMessages(library(org.Hs.eg.db))
        suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
        suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
        hs <- org.Hs.eg.db
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
        genome <- BSgenome.Hsapiens.UCSC.hg19
        }
      else if(refGeneome == "hg38") 
        {
        suppressMessages(library(org.Hs.eg.db))
        suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
        suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
        hs <- org.Hs.eg.db
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        genome <- BSgenome.Hsapiens.UCSC.hg38
        }
  }
  #
  return(list(hs,genome,txdb))    
}
###
###
###
get_gene_info = function(HGNC_symbol, refGenome="hg19", annotationDb="ensemble", txPerGeneMethod=1, txGenePairs){
  libs   <- import_annotation_db(refGenome, annotationDb)
  #
  hs     <- libs[[1]]
  genome <- libs[[2]]
  txdb   <- libs[[3]]
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
      geneInfo <- data.frame(geneGR)
      geneInfo <- geneInfo[, -c(4)]
      geneInfo <- geneInfo[, c(1:4)]
      colnames(geneInfo) = c("seqnames", "start", "end","strand") # Transcript (Including UTRs)
      ##
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
      ##
      ## exons ##
      # geneExonsByTx <- exonsBy(hs, "tx", filter = TxNameFilter(geneTxs$tx_id))
      # exReduced <- GenomicRanges::reduce(unlist(geneExonsByTx))
      # seqlevelsStyle(exReduced) <- "UCSC"
      ##
      ##
      ## cds ##
      geneCdsByTx <- cdsBy(hs, "tx", filter = TxNameFilter(geneTxs$tx_id))
      cdsReduced <- GenomicRanges::reduce(unlist(geneCdsByTx))
      seqlevelsStyle(cdsReduced) <- "UCSC"
      ##
      ## 
      cdsInfo <- data.frame(cdsReduced)
      cdsInfo <- cdsInfo[,c(1:3,5)] 
      return(cdsInfo)
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
      suppressMessages(geneInfo <- genes(txdb, columns="gene_id"))
      geneInfo <- geneInfo[geneInfo$gene_id == geneNo]
      geneInfo <- data.frame(geneInfo)
      geneInfo <- geneInfo[, -c(4)]
      geneInfo <- geneInfo[, c(1:4)]
      colnames(geneInfo) = c("seqnames", "start", "end","strand") # Transcript (Including UTRs)
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
      cdsInfo <- data.frame(cdsReduced)
      cdsInfo <- cdsInfo[,c(1:3,5)] 
      return(cdsInfo)    }
}
###
###
###
Genes_in_range = function(chromosome, leftBorder, rightBorder, refGenome="hg19", annotationDb="ensemble") {
  #
  libs <- import_annotation_db(refGenome, annotationDb)
  #
  hs     <- libs[[1]]
  genome <- libs[[2]]
  txdb   <- libs[[3]]
  #
  grRanges <-  GRanges(seqnames = Rle(paste0("chr",chromosome)),
                       ranges=IRanges(start=leftBorder, end = rightBorder))
  #
  if(annotationDb == "ensemble") 
  {    
    genes <- genes(hs, filter=list(SeqNameFilter(chromosome)), columns=c("gene_name"))
    seqlevelsStyle(genes) <- "UCSC"
    #
    genesIntersected  <- gr_intersect(genes,grRanges)
    genes_HGNC_symbol <- unique(genesIntersected$gene_name)
  }
  else if(annotationDb == "UCSC")
  {
    genes <- genes(txdb, single.strand.genes.only=FALSE)
    #
    intersectIndecies <- data.frame(findOverlaps(query=grRanges, subject=genes,ignore.strand=TRUE))$subjectHits
    genesNos <- unique(data.frame(genes[intersectIndecies,])$group_name)
    hsSub <- ensembldb::select(hs,columns = c("SYMBOL", "ENSEMBL"), keys = genesNos)
    genes_HGNC_symbol <- unique(hsSub$SYMBOL)
  }
  #
  return(genes_HGNC_symbol)
}
###
###
###
get_neighbour_genes_granges = function(chromosome, start, end, HGNC_symbol, neighbors_window = 1e+06, outlierNeighborsThreshold = 0.2, onlyNames = F, providePromoters = F, refGenome="hg19", annotationDb="ensemble", txPerGeneMethod=1, txGenePairs) {
  #
  each_side_length  <- neighbors_window
  getAtLeastOneGene <- FALSE
  #
  leftmost_coordinate   <- start-each_side_length
  leftTargetGeneBorder  <- start-1
  rightmost_coordinate  <- end+each_side_length
  rightTargetGeneBorder <- end+1
  # chromosome length
  chromosome_length <- seqlengths(Hsapiens)[[paste("chr", chromosome, sep="")]]
  ##
  ##
  ##
  left_genes  <- Genes_in_range(chromosome, leftmost_coordinate, leftTargetGeneBorder, refGenome, annotationDb)
  right_genes <- Genes_in_range(chromosome, rightTargetGeneBorder, rightmost_coordinate, refGenome, annotationDb)
  ##
  ##
  ##
  # searching for more genes if no genes in the used window
  if(getAtLeastOneGene){
      if (length(left_genes) == 0) {lneighb_genes ="not exist"} else{lneighb_genes ="exist"}
      if (length(right_genes) == 0) {rneighb_genes ="not exist"} else{rneighb_genes ="exist"}
      #
      while (lneighb_genes == "not exist") {
        leftmost_coordinate = leftmost_coordinate-(each_side_length/5)
        if (leftmost_coordinate >= 0) {
          left_genes = Genes_in_range(chromosome, leftmost_coordinate, leftTargetGeneBorder, refGenome, annotationDb)
          if (length(left_genes) != 0) {
            left_genes = left_genes[length(left_genes)]
            lneighb_genes ="exist" 
          }
        }    else {
          leftmost_coordinate = 0 
          left_genes = Genes_in_range(chromosome, leftmost_coordinate, leftTargetGeneBorder, refGenome, annotationDb)
          if (length(left_genes) != 0) {
            left_genes = left_genes[length(left_genes)]
            lneighb_genes ="exist" 
          } else {
            lneighb_genes ="not exist but no more nts"
          }
        }
      }
      #
      while (rneighb_genes == "not exist") {
        rightmost_coordinate = rightmost_coordinate+(each_side_length/5)
        if (rightmost_coordinate <= chromosome_length) {
          right_genes = Genes_in_range(chromosome, rightTargetGeneBorder, rightmost_coordinate, refGenome, annotationDb)
          if (length(right_genes) != 0) {
            right_genes = right_genes[1]
            rneighb_genes ="exist" 
          }
        }    else {
          rightmost_coordinate = chromosome_length
          right_genes = Genes_in_range(chromosome, rightTargetGeneBorder, rightmost_coordinate, refGenome, annotationDb)
          if (length(right_genes) != 0) {
            right_genes = right_genes[1]
            rneighb_genes ="exist" 
          } else {
            rneighb_genes ="not exist but no more nts" 
          }
        }
      }
  }
  ##
  ##
  # getting exons & promoters of neighbor genes
  genes <- setdiff(unique(c(left_genes,right_genes)), HGNC_symbol)
  # print(genes)
  if (length(genes) != 0) {
    if (! onlyNames) {
      genes_info <- lapply(genes, get_gene_info, refGenome=refGenome, annotationDb=annotationDb, txPerGeneMethod = txPerGeneMethod, txGenePairs = txGenePairs)
      #
      cds_info <- lapply(genes_info, function(x) rbindlist(x[2]))
      cds_info <- rbindlist(cds_info)
      cds_info <- GenomicRanges::reduce(GRanges(cds_info), ignore.strand=TRUE)              
      #
      promoters_info <- lapply(genes_info, function(x) rbindlist(x[3]))
      promoters_info <- rbindlist(promoters_info)
      promoters_info <- GenomicRanges::reduce(GRanges(promoters_info), ignore.strand=TRUE)              
      #
      if(providePromoters){
        return(list(cds_info, promoters_info))
      }
      else{
        return(cds_info)
      }      
    } 
    else 
    {
      return(genes)
    }
  } 
  else 
  {
    return(NULL)
  }
}
###
###
###
# normalize each row to sum to 1. There may be a faster way.
rowNorm =function(m) {
  t( apply(m, 1, function(x) { x/sum(x) } ) )
}
##
##
euclidean = function(a, b) {
  sqrt(sum((a - b)^2))
} 
##
##
intron_sampling = function(countsAll){
  ####
  ####
  # Exlude composition in which introns have more nucloeotide than exons as we don't exclude any exonics bps
  # x <- as.vector(countsAll["exon",] < countsAll["intron",])
  # counts <- countsAll[,x]
  # countsUn <- countsAll[,!x] 
  counts <- countsAll
  ####
  ####
  maxInExRatio = 10
  stoppingCriterion = 0.01 
  maxIter = min(100*length(counts),10000)  # to prevent endless loops (in reality this can be a very high #, this works quite fast)
  iter=0
  ####
  ####  
  while ( TRUE ) {
    #
    iter=iter+1;
    #
    # setting the 'target' row to 1 (exons)
    # setting the 'offender' row to 2, so the changes will done only in introns
    target   = 1;
    offender = 2;
    #
    if(sum(counts[offender,]) <= 0){
      cat( sprintf("Cannot continue optimization - no more nts at introns \n"));
      break;
    }
    ##############################################################################
    # compute the frequencies & setting the "meanFreqs" to the exone frequencies #
    freqs=rowNorm(counts);
    meanFreqs=freqs[target,];
    diffs = meanFreqs - freqs[offender,];
    #
    # in that row, find the column which is most responsible for the difference
    # however importantly we care ONLY about the negative differences in this vector!
    # i.e. those are the cases where the offending row has HIGHER freqs
    # (meaning we can correct that by removing sites... we can't add sites!!)
    correctableCol = which.min(diffs)
    #
    if (counts[offender, correctableCol]<=0) {
      cat( sprintf("Cannot continue optimization - counts exhausted at row %d col %d\n", offender, correctableCol) );
      cat(sprintf(counts))
      cat(sprintf(freq))
      # note that I think we could continue... just mark this cell in the matrix as uncorrectable and move to another? not sure.
      break;
    }
    #
    # computing tolerance (note that tolerance is expressed via worstCol not via correctableCol -- I am not sure if that is correct/optimal)
    # 1) worst-column tolerance
    # worstCol  = which.max(abs(diffs))  # this is sometimes the same as the correctable col
    # tolerance = abs(diffs[worstCol]) 
    # 2) Eucliean-based distance
    tolerance = euclidean(freqs[offender,],freqs[target,])
    #    
    #
    ####################
    # print some stats #
    #cat( sprintf( "Iteration: %7d, Euclidean: %.3f, Tolerance: %.3f, ShortestWin: %7d, AvgWin: %7d\n",
    #              iter, euclidean(freqs[offender,],meanFreqs), abs(diffs[worstCol]), min( rowSums(counts) ), round(mean( rowSums(counts) ))  ) )
    if(iter == 1){
      initTolerance = tolerance
    }
    #
    # did we reduce the difference enough?
    if (tolerance <= stoppingCriterion) {
      cat(sprintf("Successfully completed optimization: tolerance reached.\n"));
      break;
    }
    #
    ###########################################################
    # note this adjustment (subtraction) is too conservative, #
    # but by iterating it should converge to the right value
    subtractThis = round( diffs[correctableCol] * sum(counts[offender, ]) )
    #
    if ( subtractThis == 0 ) {
      cat(sprintf("Successfully completed optimization: count reduction <=0.5.\n"));
      break;
    }
    # now simply decrease counts in the responsible column to get closer to the mean
    counts[offender, correctableCol] = counts[offender, correctableCol] + subtractThis
    #
    # number of iterations
    if (iter==maxIter) {
      cat( sprintf("Stopping optimization - maximum number of iterations reached.\n") );
      break;
    }  
  }
  ###
  ###
  cat(sprintf("tolerance: %.3f \n", tolerance))
  targetWidth   <- sum(counts[target,])
  offenderWidth <- sum(counts[offender,])
  bgTrRatio <- offenderWidth/targetWidth  
  if(bgTrRatio > maxInExRatio){
      counts[offender,] <- round(counts[offender,]*maxInExRatio/bgTrRatio)
  }
  ###
  ###
  #countsAll <- cbind(counts,countsUn)
  countsAll <- counts
  return(list(initTolerance, tolerance, countsAll))
} 
#
################################################################################
### --------------------- Base-pair representation setup ------------------- ###
###
###
###
getObjBps = function(grObj, refGenome = "hg19", isFlank = 0){
  ##
  if(refGenome == "hg19"){suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
  } else if(refGeneome == "hg38"){suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))}
  ##
  if(length(grObj) == 0){
    return (grObj)
  } 
  else{
    #
    dfObj <- data.frame(grObj)
    ###
    inSeqName <- dfObj$seqnames[1] 
    inStrand  <- dfObj$strand[1]
    objBps    <- unlist(as.vector(mapply(function(x,y){seq(x,y)}, dfObj$start, dfObj$end)))
    #
    if(isFlank == 0){
      objBpsDis <- unlist(as.vector(mapply(function(x,y){pmin(seq(x,y) - x, y -  seq(x,y))}, dfObj$start, dfObj$end)))
    } else {
      dis1 <- dfObj$end[1] - seq(dfObj$start[1], dfObj$end[1])
      dis2 <- seq(dfObj$start[2], dfObj$end[2]) - dfObj$start[2]
      objBpsDis <- unlist(as.vector(c(dis1,dis2)))            
    }
    #barplot(objBpsDis)
    objBpsDF <- data.frame(seqnames = inSeqName, start = objBps, end = objBps, strand = inStrand, distance = objBpsDis, context1 = "NNN", context2 = "NNNNN")
    grObjBps <- GRanges(objBpsDF)
    ###
    ###
    # Tri-nucleotide context
    grObjBpsContext <- grObjBps
    start(grObjBpsContext) <- start(grObjBpsContext) - 1
    end(grObjBpsContext)   <- end(grObjBpsContext) + 1
    #
    grObjBps$context1 <- getSeq(Hsapiens,grObjBpsContext)
    grObjBps$context1[(substr(grObjBps$context1, 2,2) %in% c("A","G")),] <- reverseComplement(grObjBps$context1[(substr(grObjBps$context1, 2,2) %in% c("A","G")),])
    grObjBps$context1 <- data.frame(grObjBps$context1)$grObjBps.context1
    #
    #
    # Penta-nucleotide context
    grObjBpsContext <- grObjBps
    start(grObjBpsContext) <- start(grObjBpsContext) - 2
    end(grObjBpsContext)   <- end(grObjBpsContext) + 2
    #
    grObjBps$context2 <- getSeq(Hsapiens,grObjBpsContext)
    grObjBps$context2[(substr(grObjBps$context2, 3,3) %in% c("A","G")),] <- reverseComplement(grObjBps$context2[(substr(grObjBps$context2, 3,3) %in% c("A","G")),])
    grObjBps$context2 <- data.frame(grObjBps$context2)$grObjBps.context2            
    #
    return(grObjBps)
  }
}
###
###
###
select_background = function(grIntrons, grFlanks, minBackgroundWidth, maxDistanceFromTarget){
  #
  totIntronsWidth     <- length(grIntrons)
  if(totIntronsWidth > 0){
    maxIntronBpDistance <- max(grIntrons$distance) 
  } else{
    maxIntronBpDistance <- 0
  }
  #
  if(totIntronsWidth > minBackgroundWidth){
      #
      if(maxIntronBpDistance <= maxDistanceFromTarget){
        # no-truncation
        return(grIntrons)
      }
      else{
        # remove nts to keep only sites next to exons
        distances <- c(seq(maxDistanceFromTarget,maxIntronBpDistance,100),maxIntronBpDistance)
        intronsAtDistance <- unlist(lapply(distances, function(x){length(grIntrons$distance[grIntrons$distance <= x])}))
        chosenDistance    <- distances[min(which(intronsAtDistance >= minBackgroundWidth))]
        #
        grIntrons <- grIntrons[grIntrons$distance <= chosenDistance]
        return(grIntrons)
      }
  } else {
      # get extra bps from flanks
      remainingBps       <- minBackgroundWidth - totIntronsWidth - 2
      totalFlanksWidth   <- length(grFlanks)
      maxFlankBpDistance <- max(grFlanks$distance)
      flankPenalty       <- 100000 # penalty to prioritize introns bps, compared to flanks bps
      #
      #
      if(totalFlanksWidth <= remainingBps){
        # merge grFlanks with grIntrons
        grFlanks$distance <- grFlanks$distance + flankPenalty
        grIntrons <- do.call(c,as(list(grIntrons,grFlanks), "GRangesList"))
        grIntrons <- GenomicRanges::sort(grIntrons)
        return(grIntrons)
        #
      } else{
        # remove nts to keep only sites next to exons
        distances <- c(seq(1,maxFlankBpDistance,100),maxFlankBpDistance)
        flanksAtDistance <- unlist(lapply(distances, function(x){length(grFlanks$distance[grFlanks$distance <= x])}))
        chosenDistance   <- distances[min(which(flanksAtDistance >= remainingBps))]
        #
        grFlanks <- grFlanks[grFlanks$distance <= chosenDistance]
        # merge grFlanks with grIntrons
        grFlanks$distance <- grFlanks$distance + flankPenalty
        grIntrons <- do.call(c,as(list(grIntrons,grFlanks), "GRangesList"))
        grIntrons <- GenomicRanges::sort(grIntrons)
        return(grIntrons)        
        # 
      }
  }
}    
###
###
###
add_CpGs_info = function(grExons,grIntrons,uCpGs){
        ##
        ##
        grExonsList <- list()
        #
        grExonsM <- gr_setdiff(grExons,uCpGs)
        if(length(grExonsM) != 0){
          grExonsList <- list(grExonsList,grExonsM)
          }
        #
        grExonsU <- gr_intersect(grExons,uCpGs)
        if(length(grExonsU) != 0){
          grExonsU$context1 <- paste0("u",grExonsU$context1)
          grExonsU$context2 <- paste0("u",grExonsU$context2)
          grExonsList <- list(grExonsList,grExonsU)
          }
        #
        if(length(grExonsList) != 0){
          grExonsCpG <- sort(do.call(c,as(unlist(grExonsList),"GRangesList")))    
        } else {
          grExonsCpG <- NULL
        }
        ##
        ###
        grIntronsList <- list()
        #
        grIntronsM <- gr_setdiff(grIntrons,uCpGs)
        if(length(grIntronsM) != 0){
          grIntronsList <- list(grIntronsList,grIntronsM)
          }
        #
        grIntronsU <- gr_intersect(grIntrons,uCpGs)
        if(length(grIntronsU) != 0){
          grIntronsU$context1 <- paste0("u",grIntronsU$context1)
          grIntronsU$context2 <- paste0("u",grIntronsU$context2)
          grIntronsList <- list(grIntronsList,grIntronsU)
          }
        #
        if(length(grIntronsList) != 0){
          grIntronsCpG <- sort(do.call(c,as(unlist(grIntronsList),"GRangesList")))    
        } else {
          grIntronsCpG <- NULL
        }
        ##
        ##
        return(list(grExonsCpG,grIntronsCpG))
}
###
###
###
add_mappability_info = function(grExons, grIntrons, mappTracks){
    #
    chromM <- as.vector(grExons@seqnames)[1]
    startM <- min(min(start(grExons)),min(start(grIntrons)))
    endM <- max(max(end(grExons)),max(end(grIntrons)))
    #
    mappTracks <- mappTracks[mappTracks$seqnames == chromM,]
    #
    startMapp <- max(mappTracks$start[mappTracks$start <= startM])
    endMapp   <- min(mappTracks$end[mappTracks$end >= endM])
    mappTracks <- mappTracks[mappTracks$start >= startMapp,]
    mappTracks <- mappTracks[mappTracks$end <= endMapp,]
    gc()
    ####
    grMapp  <- GRanges(mappTracks)
    grMappL <- grMapp[grMapp$Mapp == "low"]
    grMappM <- grMapp[grMapp$Mapp == "medium"]
    grMappH <- grMapp[grMapp$Mapp == "high"]
    ####
    grExonsList <- list()
    #
    grExonsL <- gr_intersect(grExons,grMappL)
    if(length(grExonsL) != 0){
      grExonsL$Mapp <- "L"
      grExonsList <- list(grExonsList,grExonsL)
      }
    # 
    grExonsM <- gr_intersect(grExons,grMappM)
    if(length(grExonsM) != 0){
      grExonsM$Mapp <- "M"
      grExonsList <- list(grExonsList,grExonsM)
      }
    #
    grExonsH <- gr_intersect(grExons,grMappH)
    if(length(grExonsH) != 0){
      grExonsH$Mapp <- "H"
      grExonsList <- list(grExonsList,grExonsH)
      }
    #
    if(length(grExonsList) != 0){
              grExonsMapp <- sort(do.call(c,as(unlist(grExonsList),"GRangesList")))    
            } else {
              grExonsMapp <- NULL
            }
    ####
    grIntronsList <- list()
    #
    grIntronsL <- gr_intersect(grIntrons,grMappL)
    if(length(grIntronsL) != 0){
      grIntronsL$Mapp <- "L"
      grIntronsList <- list(grIntronsList,grIntronsL)
      }
    #
    grIntronsM <- gr_intersect(grIntrons,grMappM)
    if(length(grIntronsM) != 0){
      grIntronsM$Mapp <- "M"
      grIntronsList <- list(grIntronsList,grIntronsM)
      }
    #
    grIntronsH <- gr_intersect(grIntrons,grMappH)
    if(length(grIntronsH) != 0){
      grIntronsH$Mapp <- "H"
      grIntronsList <- list(grIntronsList,grIntronsH)
      }
    #
    if(length(grIntronsList) != 0){
              grIntronsMapp <- sort(do.call(c,as(unlist(grIntronsList),"GRangesList")))    
            } else {
              grIntronsMapp <- NULL
            }
    ####
    return(list(grExonsMapp,grIntronsMapp))
}
###
###
###
match_target_background = function(grTarget, grBackground, tri_or_penta = 1, use_Mappability = 1){
  #
  if(use_Mappability == 0){ #Strafiction by tri-nts or penta-nts only
        if(tri_or_penta == 1){
          #
            grTarget$context <- grTarget$context1
            grBackground$context  <- grBackground$context1
        } else {
          #
            grTarget$context <- grTarget$context2
            grBackground$context  <- grBackground$context2
        }
  } else{
        if(tri_or_penta == 1){
          #
            grTarget$context <-  paste0(grTarget$context1,grTarget$Mapp)
            grBackground$context  <- paste0(grBackground$context1,grBackground$Mapp)
        } else {
          #
            grTarget$context <-  paste0(grTarget$context2,grTarget$Mapp)
            grBackground$context  <- paste0(grBackground$context2,grBackground$Mapp)
        }
        grTarget$Mapp      <- NULL
        grBackground$Mapp  <- NULL      
  }      
  #
  grTarget$context1     <- NULL
  grTarget$context2     <- NULL
  grBackground$context1 <- NULL
  grBackground$context2 <- NULL
  #
  targetComposition     <- data.frame(rbind(table(grTarget$context)))
  backgroundComposition <- data.frame(rbind(table(grBackground$context)))
  #
  targetBackground      <- rbind.fill(targetComposition, backgroundComposition)
  targetBackground[is.na(targetBackground)] <- 0
  rownames(targetBackground) <- c("exon","intron")
  #  
  #  
  newTargetBackground <- intron_sampling(targetBackground)[[3]]
  newBackgroundComposition <- newTargetBackground["intron",]
  newBackgroundComposition <- newBackgroundComposition[,newBackgroundComposition > 0]
  ##
  ##
  # select background bps next to exons based on distances
  grBackgroundsSubsets <- mapply(function(x,y){
        grBackgroundSubset <- grBackground[grBackground$context == x]
        z <- sort(grBackgroundSubset$distance, decreasing = F)[y] # distance of y minimum items
        grBackgroundSubsetSelected <- grBackgroundSubset[grBackgroundSubset$distance <= z]
        return(grBackgroundSubsetSelected)
        }, 
        colnames(newBackgroundComposition),
        newBackgroundComposition)
  #
  newGrBackground <- unlist(as(grBackgroundsSubsets, "GRangesList"))
  #
  return(newGrBackground)
}
###
###
###