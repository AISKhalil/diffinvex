#' Compute Intersection of Two GRanges Objects
#'
#' This function calculates the intersection of two genomic ranges (`GRanges`) objects. It identifies 
#' overlaps between the two sets of ranges and returns a subset of the first `GRanges` object 
#' that overlaps with the second.
#'
#' @param gr1 A `GRanges` object representing the first set of genomic ranges.
#' @param gr2 A `GRanges` object representing the second set of genomic ranges.
#'
#' @return A `GRanges` object containing the ranges from `gr1` that overlap with any range in `gr2`.
#'
#' @import GenomicRanges
#' @export
gr_intersect <- function(gr1,gr2){
  # subsetByOverlaps(gr1,gr2) can be used
  return(gr1[queryHits(findOverlaps(gr1, gr2, type="any", select="all",ignore.strand=TRUE)),])
}



#' Compute Set Difference of Two GRanges Objects
#'
#' This function calculates the set difference between two genomic ranges (`GRanges`) objects. 
#' It returns a subset of the first `GRanges` object that does not overlap with any range in the second.
#'
#' @param gr1 A `GRanges` object representing the first set of genomic ranges.
#' @param gr2 A `GRanges` object representing the second set of genomic ranges.
#'
#' @return A `GRanges` object containing the ranges from `gr1` that do not overlap with any range in `gr2`.
#'
#' @import GenomicRanges
#' @export
gr_setdiff <- function(gr1,gr2){
  maxIndexInChr <- 300000000 #it is added to avoid a problem if there are not intersection
  return(gr1[-c(queryHits(findOverlaps(gr1, gr2, type="any", select="all",ignore.strand=TRUE)), maxIndexInChr),])
}  



#' Read Genomic Ranges from File
#'
#' This function reads genomic range data from a file and converts it into a `GRanges` object. 
#' If the input file is not provided or is `NULL`, it returns an empty `GRanges` object.
#'
#' @param inputFile A string representing the path to the input file containing genomic range data. 
#' The file should have a header.
#'
#' @return A `GRanges` object containing the genomic ranges read from the file. If the file is `NULL`, 
#' returns an empty `GRanges` object.
#'
#' @import data.table
#' @import GenomicRanges
#' @export
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



#' Load Genomic Tracks
#'
#' This function loads genomic tracks from an RDS file based on the specified tool directory, reference genome, 
#' bandwidth, and clustering option. It reads the genomic regions data and prepares a list containing target, 
#' background, and union of these tracks.
#'
#' @param tool_directory A string representing the directory path where the tool's reference data is located.
#' @param reference_genome A string specifying the reference genome version (e.g., "hg38").
#' @param DiffInVEx_BW A numeric value indicating the bandwidth in kilobases used to select the genomic regions.
#' @param DiffInVEx_cluster An integer (0, 1, or other) specifying the type of clustering to use for background tracks: 
#' 1 for "background_tri", 0 for "background_penta", and any other value for "background_NoMatching".
#'
#' @return A list containing three elements: 
#' - `target`: a `GRanges` object representing the target genomic tracks.
#' - `background`: a `GRanges` object representing the background genomic tracks based on the clustering option.
#' - `gene`: a `GRanges` object representing the reduced union of target and background tracks.
#' 
#' If the specified RDS file does not exist, the function returns an empty list.
#'
#' @import GenomicRanges
#' @export
load_tracks <- function(tool_directory,reference_genome,DiffInVEx_BW,DiffInVEx_cluster){
      #
      tracksFile <- paste0(tool_directory,"/data/GeneData/geneDB_",reference_genome,"_",DiffInVEx_BW,"Kb","_GenomicRegions.RDS")
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



#' Load Sample Annotation
#'
#' This function reads a genome-wide annotation file and standardizes the column name for sample identifiers.
#' It ensures that the sample column is consistently named 'Sample'.
#'
#' @param annotation_genome_wide_file A string representing the file path to the genome-wide annotation data.
#' The file should be in a format readable by `fread`.
#'
#' @return A data table containing the sample annotation data, with the column for sample identifiers 
#' standardized to 'Sample'.
#'
#' @import data.table
#' @export
load_sample_annotation <- function(annotation_genome_wide_file){
  #
  sample_annotation <- fread(annotation_genome_wide_file)
  names(sample_annotation)[names(sample_annotation) == 'sample' | names(sample_annotation) == 'SAMPLE'] <- 'Sample'
  return(sample_annotation)
}



#' Convert Multi-base Substitutions to Single-base Substitutions
#'
#' This function processes a table of SNPs (Single Nucleotide Polymorphisms) to separate 
#' multi-base substitutions (MBSs) from single-base substitutions (SBSs) and converts 
#' MBSs into individual SBSs.
#'
#' @param SNP_table A data frame or data table containing SNP data. It must have columns 
#' named `chrom`, `pos`, `ref`, `alt`, and `Sample`, where `ref` and `alt` are strings 
#' representing reference and alternate alleles, respectively.
#'
#' @return A data table containing single-base substitutions, including both originally 
#' single-base substitutions and those derived from multi-base substitutions. The result 
#' is sorted by chromosome and position, and duplicates are removed.
#'
#' @import dplyr
#' @import data.table
#' @export
MBSs_To_SBSs <- function(SNP_table){
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



#' Convert SNP Data to GRanges Object
#'
#' This function reads a file containing SNP mutation data, filters and processes the data, and 
#' converts it into a `GRanges` object. The resulting `GRanges` object represents SNP mutations 
#' with additional annotations, and optionally intersects with provided genomic regions of interest.
#'
#' @param mutation_file A string representing the path to the SNP mutation file. The file should 
#' contain columns for chromosome (`chrom`), position (`pos`), reference allele (`ref`), 
#' alternate allele (`alt`), and sample identifier (`Sample`).
#' 
#' @param geneTrack A `GRanges` object representing the genomic regions of interest. 
#' If empty, no intersection is performed.
#' 
#' @param sample_annotation A data frame or data table containing sample annotations. It must 
#' include a column named `Sample` to filter the SNPs based on the provided sample annotations.
#'
#' @return A `GRanges` object containing the filtered and processed SNP mutations. The object 
#' includes information on sample identifiers, reference alleles, and mutated alleles. The ranges 
#' are intersected with the provided `geneTrack` if it is not empty.
#'
#' @import dplyr
#' @import data.table
#' @import GenomicRanges
#' @export
MkGRanges_SNPs <- function(mutation_file, geneTrack, sample_annotation) {
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



#' Add Mutation Load Information to Sample Annotation
#'
#' This function calculates the mutation load (number of mutations per kilobase) for both 
#' exonic and intronic regions and integrates this information into the sample annotation data.
#' The mutation load is calculated based on the provided genomic regions and is added to the 
#' sample annotation data frame.
#'
#' @param grSNPs A `GRanges` object containing SNP mutation data. It includes sample identifiers 
#' and mutation positions.
#' 
#' @param sample_annotation A data frame or data table containing sample annotations. It should 
#' have a column named `Sample`, and optionally a column named `tumorType` for tumor type annotations.
#'
#' @param controlled_variables_input A string representing the file path to a file containing 
#' controlled variables. This file should list variables for which mutation load information is to 
#' be computed. The presence of `MutationLoad` in this file will trigger the mutation load calculation.
#'
#' @param genomeExonicRegions A `GRanges` object representing exonic regions.
#'
#' @param genomeIntronicRegions A `GRanges` object representing intronic regions.
#'
#' @return A data frame or data table containing the updated sample annotations with added mutation 
#' load columns (`MutationLoad_Tr` for exonic regions and `MutationLoad_Bg` for intronic regions). 
#' Mutation load is expressed as the number of mutations per kilobase of the genomic region.
#'
#' @import dplyr
#' @import data.table
#' @import GenomicRanges
#' @export
add_MutationLoad <- function(grSNPs,sample_annotation,controlled_variables_input,genomeExonicRegions,genomeIntronicRegions){
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


#' Load Genomic Tracks for Analysis
#'
#' This function loads various genomic tracks based on the reference genome specified. It includes
#' tracks related to genomic repeats, CTCF binding sites, blacklist regions, mappability, background
#' hotspots, methylation sites, promoter motifs, and somatic hypermutation (SHM) regions. Additionally,
#' it retrieves and processes coding regions to exclude from introns and flanks.
#'
#' @param reference_genome A string indicating the reference genome. Supported values are `"hg19"` 
#' and `"hg38"`. This determines which set of genomic track files to load.
#'
#' @param toolDirectory A string representing the directory path where genomic track files are located.
#'
#' @param filtering_tracks A string representing the file path to a file containing regions for filtering. 
#' The file should be in a format compatible with the `read_GRanges` function.
#'
#' @param filter_inclusion A boolean value indicating whether to include or exclude the regions specified 
#' in `filtering_tracks` from the analysis.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item `mappRegions`: A `GRanges` object with mappability regions.
#'   \item `mappScores`: A data frame with mappability scores.
#'   \item `filterRegions`: A `GRanges` object with regions specified for filtering (if applicable).
#'   \item `filterCondition`: A boolean indicating whether filtering regions are used.
#'   \item `filter_inclusion`: The value of the `filter_inclusion` parameter.
#'   \item `BlackListed_regions`: A `GRanges` object with blacklist regions.
#'   \item `Microsatellites`: A `GRanges` object with microsatellite regions.
#'   \item `CTCF_regions`: A `GRanges` object with CTCF binding sites.
#'   \item `backgroundHotspots`: A `GRanges` object with background hotspot regions.
#'   \item `uCpGs`: A `GRanges` object with unmethylated CpG sites.
#'   \item `promoterMotifs`: A `GRanges` object with promoter motifs.
#'   \item `shmRegions`: A `GRanges` object with combined SHM on and off targets.
#'   \item `allTargetRegions`: A `GRanges` object with all coding regions, processed and saved.
#' }
#'
#' @importFrom GenomicRanges import GRanges keepStandardChromosomes start end union reduce sort
#' @importFrom data.table fread
#' @importFrom methods is
#' @export
read_genomic_tracks <- function(reference_genome, 
                                toolDirectory, 
                                filtering_tracks, 
                                filter_inclusion){
  ##
  ##  auxiliary files
  if(reference_genome == "hg19"){
      MicrosatellitesTracks  <- paste0(toolDirectory,"/data/GenomicTracks/Repeats/hg19_repeats_more6nts.bed")
      CTCF_Tracks            <- paste0(toolDirectory,"/data/GenomicTracks/CTCF_regions/CTCF_Cohesion_peaks_RenLab_Motifs_hg19.bed")
      BlackListedTracks      <- paste0(toolDirectory,"/data/GenomicTracks/BlackListedRegions/hg19_consensusBlacklist.v3.bed")
      mappabilityTracks      <- paste0(toolDirectory,"/data/GenomicTracks/MappTracks/uMap/UMAP_hg19_all_mappability_100bpTrack.csv")
      mappabilityScores      <- paste0(toolDirectory,"/data/GenomicTracks/MappTracks/uMap/UMAP_hg19_high_medium_low_mappability_track.csv")
      backgroundHotspotsFile <- paste0(toolDirectory,"/data/GenomicTracks/Background_hotspots/hg19_background_hotspots_min3muts.csv")
      uCpGsFile              <- paste0(toolDirectory,"/data/GenomicTracks/Methylation/hg19_uMethylated_CpGs_solid_tumors.csv")
      promoterMotifFile      <- paste0(toolDirectory,"/data/GenomicTracks/PromotersMotifs/hg19_promoter_NYTTCCG_regions.csv")
      shmOnFile              <- paste0(toolDirectory,"/data/GenomicTracks/SHM/hg19_shmOnTargets_PROCESSED_CRG75.bed") 
      shmOffFile             <- paste0(toolDirectory,"/data/GenomicTracks/SHM/hg19_shmOffTargets_PROCESSED_CRG75.bed")
      all_cds_file           <- paste0(toolDirectory,"/data/GeneData/geneDB_hg19_ensemble_cdss.RDS")

  } else if(reference_genome == "hg38"){ #To-be-updated
      MicrosatellitesTracks  <- paste0(toolDirectory,"/data/GenomicTracks/Repeats/hg19_repeats_more6nts.bed")
      CTCF_Tracks            <- paste0(toolDirectory,"/data/GenomicTracks/CTCF_regions/CTCF_Cohesion_peaks_RenLab_Motifs_hg19.bed")
      BlackListedTracks      <- paste0(toolDirectory,"/data/GenomicTracks/BlackListedRegions/hg38_consensusBlacklist.v3.bed")
      mappabilityTracks      <- paste0(toolDirectory,"/data/GenomicTracks/MappTracks/uMap/UMAP_hg19_all_mappability_100bpTrack.csv")
      mappabilityScores      <- paste0(toolDirectory,"/data/GenomicTracks/MappTracks/uMap/UMAP_hg19_high_medium_low_mappability_track.csv")
      backgroundHotspotsFile <- paste0(toolDirectory,"/data/GenomicTracks/Background_hotspots/hg19_background_hotspots_min3muts.csv") 
      uCpGsFile              <- paste0(toolDirectory,"/data/GenomicTracks/Methylation/hg19_uMethylated_CpGs_solid_tumors.csv")
      promoterMotifFile      <- paste0(toolDirectory,"/data/GenomicTracks/PromotersMotifs/hg19_promoter_NYTTCCG_regions.csv")
      shmOnFile              <- paste0(toolDirectory,"/data/GenomicTracks/SHM/hg19_shmOnTargets_PROCESSED_CRG75.bed") 
      shmOffFile             <- paste0(toolDirectory,"/data/GenomicTracks/SHM/hg19_shmOffTargets_PROCESSED_CRG75.bed")
      all_cds_file           <- paste0(toolDirectory,"/data/GeneData/geneDB_hg38_ensemble_cdss.RDS")
  } else {
      stop("Error: Reference genome should be hg19 or hg38")
  }
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



#' Integrate Gene Annotations into Sample Metadata
#'
#' This function merges gene-specific annotations into a sample metadata table based on a specified gene symbol.
#' It reads gene annotations from a file and integrates them into an existing sample annotation table, 
#' retaining only samples that have annotations for the specified gene.
#'
#' @param sample_annotation A data frame containing sample metadata. Must include a column named `"Sample"`.
#'
#' @param sample_annotation_perGene_file A string representing the file path to a file containing gene annotations. 
#' The file should have columns for sample identifiers and gene names. Columns can be named `"sample"` or `"SAMPLE"` 
#' for sample identifiers and `"gene"` or `"GENE"` for gene names.
#'
#' @param HGNC_symbol A string representing the HGNC symbol of the gene for which annotations should be retrieved.
#'
#' @return A data frame with the original sample metadata, augmented with gene-specific annotations for the given HGNC symbol. 
#' Only samples present in both the original metadata and the gene annotations will be included.
#'
#' @importFrom data.table fread
#' @export
load_gene_annotation <- function(sample_annotation,
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



#' Filter Samples Based on Exonic Mutation Count
#'
#' This function filters out samples from a set of SNP data and annotations if they have 3 or more mutations located in exonic regions.
#' It identifies and excludes samples based on the number of exonic mutations and updates both the SNP data and the annotations accordingly.
#'
#' @param grSNPs A `GRanges` object containing SNP data with a column named `"Sample"` indicating sample identifiers.
#'
#' @param annotations A data frame with sample metadata, including a column named `"Sample"` that matches the identifiers in `grSNPs`.
#'
#' @param grExons A `GRanges` object representing exonic regions to be used for filtering SNPs.
#'
#' @return A list containing:
#' \describe{
#'   \item{grSNPs}{A `GRanges` object with SNP data, excluding samples with 3 or more exonic mutations.}
#'   \item{annotations}{A data frame with sample metadata, excluding samples with 3 or more exonic mutations.}
#' }
#'
#' @importFrom GenomicRanges findOverlaps queryHits
#' @importFrom dplyr %>%
#' @export
sample_wise_filtering <- function(grSNPs, annotations, grExons){
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



#' Write Mutation Data to File Based on Genomic Regions
#'
#' This function filters and writes mutation data to a file based on their overlap with exonic and intronic regions. 
#' The data is categorized and saved into separate files depending on whether mutations fall within exonic or intronic regions.
#'
#' @param directory A character string specifying the directory where the output files will be saved. If the directory does not exist, it will be created.
#' 
#' @param HGNC_symbol A character string representing the HGNC symbol for the gene associated with the mutations.
#' 
#' @param grSNPs A `GRanges` object containing SNP data with columns `seqnames`, `start`, `end`, `ref_allele`, `mutated_to_allele`, `Sample`.
#' 
#' @param grExons A `GRanges` object representing exonic regions for filtering SNPs.
#' 
#' @param grIntrons A `GRanges` object representing intronic regions for filtering SNPs.
#' 
#' @param Method A character string specifying the method or type of analysis for naming the output file.
#'
#' @return A logical value indicating whether any SNPs were present in either exonic or intronic regions. Returns `TRUE` if SNPs were found in either region, otherwise `FALSE`.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Checks for overlaps between SNPs and exonic or intronic regions.
#'   \item Creates data frames for SNPs in exonic and intronic regions, including additional columns for gene symbol, class, and region width.
#'   \item Combines the data frames and writes them to a tab-separated values (TSV) file in the specified directory.
#' }
#'
#' @importFrom data.table fwrite
#' @importFrom GenomicRanges gr_intersect
#' @export
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