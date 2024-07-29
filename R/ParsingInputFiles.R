###
"%!in%" <- function(x, y) !("%in%"(x, y))
###
###
###
gr_intersect = function(gr1,gr2){
  # subsetByOverlaps(gr1,gr2) can be used
  return(gr1[queryHits(findOverlaps(gr1, gr2, type="any", select="all",ignore.strand=TRUE)),])
}  
#
gr_setdiff = function(gr1,gr2){
  maxIndexInChr <- 300000000 #it is added to avoid a problem if there are not intersection
  return(gr1[-c(queryHits(findOverlaps(gr1, gr2, type="any", select="all",ignore.strand=TRUE)), maxIndexInChr),])
}  
###
###
###
read_GRanges <- function(inputFile) {
  if (! is.null(inputFile)) {
    regions <- fread(inputFile, header = TRUE)
    regions <- GRanges(regions)
  } 
  else {
    regions <- GRanges(NULL)
  }
  return(regions)
}
###
###
###
load_tracks = function(tool_directory,reference_genome,DiffInVEx_BW,DiffInVEx_cluster){
  #
  tracksFile <- paste0(tool_directory,"/references/GeneData/geneDB_",reference_genome,"_",DiffInVEx_BW,"Kb","_GenomicRegions.RDS")
  #
  if(file.exists(tracksFile)){
            x <- readRDS(tracksFile)
            #
            trTracks <- x[["target"]]
            if(DiffInVEx_cluster == 1){
              bgTracks <- x[["background_tri"]]
            } else if(DiffInVEx_cluster == 0){
              bgTracks <- x[["background_penta"]]
            } else{
              bgTracks <- x[["background_NoMatching"]]
            }
            #
            gTracks <- GenomicRanges::reduce(GenomicRanges::union(trTracks,bgTracks))
            #
            grDict <- list()
            grDict[["target"]]      <- trTracks
            grDict[["background"]]  <- bgTracks
            grDict[["gene"]]        <- gTracks
    } else{
            #
            grDict <- list()

    }
  #
  return(grDict)
}
###
###
###
load_sample_annotation = function(annotation_genome_wide_file){
  #
  sample_annotation <- fread(annotation_genome_wide_file)
  names(sample_annotation)[names(sample_annotation) == 'sample' | names(sample_annotation) == 'SAMPLE'] <- 'Sample'
  return(sample_annotation)
}
###
###
###
MBSs_To_SBSs = function(SNP_table){
  SBSs_table <- SNP_table %>% dplyr::filter(nchar(ref) == 1)
  MBSs_table <- SNP_table %>% dplyr::filter(nchar(ref) != 1)
  #
  if(nrow(MBSs_table) != 0){
            x <- MBSs_table
            REF <- strsplit(x$ref, split ="")
            ALT <- strsplit(x$alt, split ="")
            #
            noSBSsPerCall <- unlist(lapply(REF, function(x){length(x)}))
            CHROM  <- rep(x$chrom, noSBSsPerCall)
            Sample <- rep(x$Sample, noSBSsPerCall)
            POS    <- mapply(function(X,Y){X + seq(1,Y)-1}, X = x$pos, Y = noSBSsPerCall)
            #
            new_MBSs_table <- data.table(chrom = CHROM, pos = unlist(POS), ref = unlist(REF), alt = unlist(ALT), Sample = Sample)
            new_MBSs_table <- new_MBSs_table %>% dplyr::filter(ref != alt)
            #
            new_SNP_table <- rbind(SBSs_table, new_MBSs_table)
            new_SNP_table <- new_SNP_table[order(new_SNP_table$chrom,new_SNP_table$pos),]
  } else{
            new_SNP_table <- SBSs_table
  }
  #
  new_SNP_table <- unique(new_SNP_table)
  #
  return(new_SNP_table)
}
###
###
###
MkGRanges_SNPs = function(mutation_file, geneTrack, sample_annotation) {
  ## check if mutations file does exist
  if(is_file(mutation_file)){
      SNP_table <- fread(mutation_file, header=T, stringsAsFactors = FALSE)  
  } else {
      stop("Mutation file does not exist")
  }
  ## 
  colnames(SNP_table) <- tolower(colnames(SNP_table))
  names(SNP_table)[names(SNP_table) == 'sample' | names(SNP_table) == 'SAMPLE'] <- 'Sample'
  ##
  ##
  # keep only mutations of samples with annotations
  SNP_table <- SNP_table %>% dplyr::filter(Sample %in% sample_annotation$Sample)
  #
  # filter_out indels
  SNP_table <- SNP_table %>% dplyr::filter(nchar(ref) == nchar(alt))
  #
  # Split MBSs into SBSs
  SNP_table <- MBSs_To_SBSs(SNP_table)
  ##
  ##
  grSNPs <- GRanges(
    seqnames = Rle(SNP_table$chrom),
    ranges   = IRanges(start=SNP_table$pos, end = SNP_table$pos),
    Sample   = SNP_table$Sample,
    ref_allele = SNP_table$ref,
    mutated_to_allele = SNP_table$alt
  )
  seqlevelsStyle(grSNPs) = "UCSC" 
  #
  # Intersect with genomic regions of interest
  if(length(geneTrack) != 0){
        grSNPs <- grSNPs[queryHits(findOverlaps(grSNPs, geneTrack, type="any", select="all",ignore.strand=TRUE)),]
  }
  cat(paste0(length(grSNPs)," mutations \n"))
  ## output
  return(grSNPs)
}
###
###
###
add_MutationLoad = function(grSNPs,sample_annotation,controlled_variables_input,genomeExonicRegions,genomeIntronicRegions){
  #
  colNamesToRegress = fread(controlled_variables_input, header = F)
  colNamesToRegress = colNamesToRegress$V1    
  #
  #  
  if("MutationLoad" %in% colNamesToRegress){
        #
        trWidth <- sum(width(genomeExonicRegions))
        bgWidth <- sum(width(genomeIntronicRegions))
        #
        #
        grSNPs1 <- gr_intersect(grSNPs,genomeExonicRegions)
        grSNPs2 <- gr_intersect(grSNPs,genomeIntronicRegions)
        #
        #
        if("tumorType" %in% colNamesToRegress){
             sample_annotation2 <- sample_annotation %>% dplyr::select(c("Sample","tumorType"))
            #
            grSNPs1T <- data.frame(Sample = grSNPs1$Sample, count = 1)
            #
            grSNPs1T <- merge(grSNPs1T,sample_annotation2, by = "Sample")
            grSNPs1TT <- grSNPs1T %>% dplyr::group_by(tumorType) %>% dplyr::summarize(MutationLoad_Tr = sum(count))
            #
            grSNPs1T <- merge(grSNPs1T,grSNPs1TT,by= "tumorType")
            grSNPs1T <- grSNPs1T %>% dplyr::select(c("Sample","MutationLoad_Tr"))
            grSNPs1T <- unique(grSNPs1T)
            ##
            ##
            grSNPs2T <- data.frame(Sample = grSNPs2$Sample, count = 1)
            grSNPs2T <- merge(grSNPs2T,sample_annotation2, by = "Sample")
            grSNPs2TT <- grSNPs2T %>% dplyr::group_by(tumorType) %>% dplyr::summarize(MutationLoad_Bg = sum(count))
            #
            grSNPs2T <- merge(grSNPs2T,grSNPs2TT,by= "tumorType")
            grSNPs2T <- grSNPs2T %>% dplyr::select(c("Sample","MutationLoad_Bg"))            
            grSNPs2T <- unique(grSNPs2T)
            ##
            ##
        } else{
            grSNPs1T <- data.frame(Sample = grSNPs1$Sample, count = 1)
            grSNPs1T <- grSNPs1T %>% dplyr::group_by(Sample) %>% dplyr::summarize(MutationLoad_Tr = sum(count))
            #
            grSNPs2T <- data.frame(Sample = grSNPs2$Sample, count = 1)
            grSNPs2T <- grSNPs2T %>% dplyr::group_by(Sample) %>% dplyr::summarize(MutationLoad_Bg = sum(count))
        }
        #
        # uncomment if single (exonic/intronic) mutation rate per study
        #grSNPs1T$MutationLoad_Tr <- sum(grSNPs1T$MutationLoad_Tr)        
        #grSNPs2T$MutationLoad_Bg <- sum(grSNPs2T$MutationLoad_Bg)
        #
        grSNPs1T$MutationLoad_Tr <- grSNPs1T$MutationLoad_Tr/trWidth*1000 # mutation rate (per KB)
        grSNPs2T$MutationLoad_Bg <- grSNPs2T$MutationLoad_Bg/bgWidth*1000 # mutation rate (per KB)
        #
        #
        mutationLoad_df <- merge(grSNPs1T,grSNPs2T, by = "Sample")
        sample_annotation <- merge(sample_annotation, mutationLoad_df, by = "Sample")
  }
  #
  return(sample_annotation)
}
###
###
###
read_genomic_tracks = function(reference_genome, 
                               toolDirectory, 
                               filtering_tracks, 
                               filter_inclusion){
  ##
  ##  auxiliary files
  if(reference_genome == "hg19"){
      MicrosatellitesTracks  <- paste0(toolDirectory,"/references/GenomicTracks/Repeats/hg19_repeats_more6nts.bed")
      CTCF_Tracks            <- paste0(toolDirectory,"/references/GenomicTracks/CTCF_regions/CTCF_Cohesion_peaks_RenLab_Motifs_hg19.bed")
      BlackListedTracks      <- paste0(toolDirectory,"/references/GenomicTracks/BlackListedRegions/hg19_consensusBlacklist.v3.bed")
      mappabilityTracks      <- paste0(toolDirectory,"/references/GenomicTracks/MappTracks/uMap/UMAP_hg19_all_mappability_100bpTrack.csv")
      mappabilityScores      <- paste0(toolDirectory,"/references/GenomicTracks/MappTracks/uMap/UMAP_hg19_high_medium_low_mappability_track.csv")
      backgroundHotspotsFile <- paste0(toolDirectory,"/references/GenomicTracks/Background_hotspots/hg19_background_hotspots_min3muts.csv")
      uCpGsFile              <- paste0(toolDirectory,"/references/GenomicTracks/Methylation/hg19_uMethylated_CpGs_solid_tumors.csv")
      promoterMotifFile      <- paste0(toolDirectory,"/references/GenomicTracks/PromotersMotifs/hg19_promoter_NYTTCCG_regions.csv")
      shmOnFile              <- paste0(toolDirectory,"/references/GenomicTracks/SHM/hg19_shmOnTargets_PROCESSED_CRG75.bed") 
      shmOffFile             <- paste0(toolDirectory,"/references/GenomicTracks/SHM/hg19_shmOffTargets_PROCESSED_CRG75.bed")
      all_cds_file           <- paste0(toolDirectory,"/references/GeneData/geneDB_hg19_ensemble_cdss.RDS")

  } else if(reference_genome == "hg38"){ #To-be-updated
      MicrosatellitesTracks  <- paste0(toolDirectory,"/references/GenomicTracks/Repeats/hg19_repeats_more6nts.bed")
      CTCF_Tracks            <- paste0(toolDirectory,"/references/GenomicTracks/CTCF_regions/CTCF_Cohesion_peaks_RenLab_Motifs_hg19.bed")
      BlackListedTracks      <- paste0(toolDirectory,"/references/GenomicTracks/BlackListedRegions/hg38_consensusBlacklist.v3.bed")
      mappabilityTracks      <- paste0(toolDirectory,"/references/GenomicTracks/MappTracks/uMap/UMAP_hg19_all_mappability_100bpTrack.csv")
      mappabilityScores      <- paste0(toolDirectory,"/references/GenomicTracks/MappTracks/uMap/UMAP_hg19_high_medium_low_mappability_track.csv")
      backgroundHotspotsFile <- paste0(toolDirectory,"/references/GenomicTracks/Background_hotspots/hg19_background_hotspots_min3muts.csv") 
      uCpGsFile              <- paste0(toolDirectory,"/references/GenomicTracks/Methylation/hg19_uMethylated_CpGs_solid_tumors.csv")
      promoterMotifFile      <- paste0(toolDirectory,"/references/GenomicTracks/PromotersMotifs/hg19_promoter_NYTTCCG_regions.csv")
      shmOnFile              <- paste0(toolDirectory,"/references/GenomicTracks/SHM/hg19_shmOnTargets_PROCESSED_CRG75.bed") 
      shmOffFile             <- paste0(toolDirectory,"/references/GenomicTracks/SHM/hg19_shmOffTargets_PROCESSED_CRG75.bed")
      all_cds_file           <- paste0(toolDirectory,"/references/GeneData/geneDB_hg38_ensemble_cdss.RDS")
  } else {
      stop("Error: Reference genome should be hg19 or hg38")
  }
  ##
  ##
  filterRegions   <- read_GRanges(filtering_tracks)
  filterCondition <- !is.null(filtering_tracks)
  ##  
  Microsatellites       <- import(MicrosatellitesTracks, format = "bed")
  Microsatellites$name  <- NULL
  Microsatellites       <- keepStandardChromosomes(Microsatellites,pruning.mode = "tidy")
  start(Microsatellites) = start(Microsatellites) - 2
  end(Microsatellites) = end(Microsatellites) + 2  
  ##
  CTCF_regions      <- import(CTCF_Tracks, format = "bed")
  CTCF_regions$name <- NULL  
  CTCF_regions$score <- NULL    
  ##
  BlackListed_regions <- import(BlackListedTracks, format = "bed")
  BlackListed_regions$name <- NULL  
  BlackListed_regions$score <- NULL     
  ##
  mappRegions <- read_GRanges(mappabilityTracks)
  mappScores  <- fread(mappabilityScores, header = TRUE)
  ##
  backgroundHotspots <- read_GRanges(backgroundHotspotsFile)      
  ##
  uCpGs <- read_GRanges(uCpGsFile)
  ##
  promoterMotifs <- read_GRanges(promoterMotifFile)
  ##
  shmOnRegions <- read_GRanges(shmOnFile)
  shmOffRegions <- read_GRanges(shmOffFile)
  shmRegions <- GenomicRanges::union(shmOnRegions,shmOffRegions)
  ##
  ## All coding-regions to be excluded from introns and flanks
  if(file.exists(all_cds_file)){
      allTargetRegions <- readRDS(all_cds_file)
  } else{
      quiet(source(paste0(toolDirectory,"/R/extract_all_cds_regions.R")), all = T)
      #
      allTargetRegions <- extract_all_cds_regions(toolDirectory, reference_genome, "ensemble")
      start(allTargetRegions) <- start(allTargetRegions) - 5
      end(allTargetRegions)   <- end(allTargetRegions) + 5
      allTargetRegions <- GenomicRanges::sort(GenomicRanges::reduce(allTargetRegions))
      #
      saveRDS(allTargetRegions,all_cds_file)
  }
  ##
  ##
  genomic_tracks <- list(mappRegions,mappScores,filterRegions,filterCondition,filter_inclusion,BlackListed_regions,Microsatellites,CTCF_regions,backgroundHotspots,uCpGs,promoterMotifs,shmRegions,allTargetRegions)
  return(genomic_tracks)
}
###
###
###
load_gene_annotation = function(sample_annotation,
                                sample_annotation_perGene_file,
                                HGNC_symbol){
  #
  sample_annotation_perGene <- fread(sample_annotation_perGene_file)
  names(sample_annotation_perGene)[names(sample_annotation_perGene) == 'sample' | names(sample_annotation_perGene) == 'SAMPLE'] <- 'Sample'
  names(sample_annotation_perGene)[names(sample_annotation_perGene) == 'gene' | names(sample_annotation_perGene) == 'GENE'] <- 'Gene'
  #
  sample_annotation_perGene <- sample_annotation_perGene[sample_annotation_perGene$Sample %in% sample_annotation$Sample,]
  #
  sample_annotation_perGene <- sample_annotation_perGene[Gene==HGNC_symbol][, -"Gene"]
  #
  sample_annotation <- merge(sample_annotation, sample_annotation_perGene, by="Sample")
  #
  return(sample_annotation)
}
###
###
###
sample_wise_filtering = function(grSNPs, annotations, grExons){
  # any sample (patient) with more than/equal to 3 exonic mutations will be filtered out
  threshold = 3 
  #
  x <- grSNPs[queryHits(findOverlaps(grSNPs, grExons, type="any", select="all",ignore.strand=TRUE)),]
  excludedSamples <- names(which(table(x$Sample) >= threshold))
  #
  grSNPs <- grSNPs[! grSNPs$Sample %in% excludedSamples]
  annotations <- annotations[!annotations$Sample %in% excludedSamples,]
  #
  return(list(grSNPs,annotations))
}
###
###
###
writing_mutations = function(directory, HGNC_symbol, grSNPs, grExons, grIntrons, Method){
    #
    if (!dir.exists(directory)){
        dir.create(directory)
    }
    #
    exonsWidth   <- sum(width(grExons))
    intronsWidth <- sum(width(grIntrons))
    #
    SNPsExonsIntersect   <- gr_intersect(grSNPs,grExons)
    SNPsIntronsIntersect <- gr_intersect(grSNPs,grIntrons)
    #
    SNPsInExons <- (sum(width(SNPsExonsIntersect)) > 0)
    SNPsInIntrons <- (sum(width(SNPsIntronsIntersect)) > 0)
    #
    if(SNPsInExons && SNPsInIntrons){
        dfExons <- as.data.frame(SNPsExonsIntersect)
        dfIntrons <- as.data.frame(SNPsIntronsIntersect)
        # 
        dfExons <- dfExons[c("seqnames","start", "end", "ref_allele", "mutated_to_allele", "Sample")]
        dfExons$DiffInVEx_gene  <- HGNC_symbol
        dfExons$DiffInVEx_class <- "target"
        dfExons$DiffInVEx_width <- exonsWidth
        #
        dfIntrons <- dfIntrons[c("seqnames","start", "end", "ref_allele", "mutated_to_allele", "Sample")]
        dfIntrons$DiffInVEx_gene  <- HGNC_symbol
        dfIntrons$DiffInVEx_class <- "background"
        dfIntrons$DiffInVEx_width <- intronsWidth
        #
        dfGene <- rbind(dfExons, dfIntrons)
        #
    } else if(SNPsInExons){
        dfExons <- as.data.frame(SNPsExonsIntersect)
        #
        dfExons <- dfExons[c("seqnames","start", "end", "ref_allele", "mutated_to_allele", "Sample")]
        dfExons$DiffInVEx_gene  <- HGNC_symbol
        dfExons$DiffInVEx_class <- "target"
        dfExons$DiffInVEx_width <- exonsWidth
        #        
        dfGene <- dfExons
        #     
    } else if(SNPsInIntrons){
        dfIntrons <- as.data.frame(SNPsIntronsIntersect)
        #
        dfIntrons <- dfIntrons[c("seqnames","start", "end", "ref_allele", "mutated_to_allele", "Sample")]
        dfIntrons$DiffInVEx_gene  <- HGNC_symbol
        dfIntrons$DiffInVEx_class <- "background"
        dfIntrons$DiffInVEx_width <- intronsWidth
        #        
        dfGene <- dfIntrons
        #
    }
    #
    if(SNPsInExons || SNPsInIntrons){
          filename = paste(directory, "/",HGNC_symbol, "_MutationProfile_", Method, ".tsv", sep = "")
          tryCatch({fwrite(dfGene, filename, sep="\t", row.names=F, quote=T) },
                error=function(err){  NULL  },  
                finally = {}
          )
    }
    #
    return(SNPsInExons || SNPsInIntrons)
}
