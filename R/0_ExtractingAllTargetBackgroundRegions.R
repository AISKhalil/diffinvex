suppressMessages(library("ddpcr"))
suppressMessages(library("optparse"))
#
packages <- c("arm","biomaRt","rtracklayer","GenomicRanges","GenomicFeatures", "dplyr", "reshape", "data.table","tidyverse","readxl","qvalue","Matrix","gplots","ggplot2","ggrepel","ggpubr")
invisible(quiet(lapply(packages, library, character.only = TRUE), all = T)) # invisible for suppressing output of the lapply, suppresMessages for libraries
#
#
#####################
### input parsing ###
#####################
option_list <- list(
  make_option(c("-r", "--reference_genome"), type="character", action="store", default="hg19",
              help="reference genome"),  
  make_option(c("-b", "--DiffInVEx_BW"), type="integer", action="store", default=20,
              help="background width (in Kb)"),
  make_option(c("-t", "--tool_directory"), type="character", action="store", default=F,
              help="Full path for the diffInVEx directory with R and auxiliary files")
)
###
###
parser <- OptionParser(option_list=option_list)
optInfo <- parse_args(parser)
#
# Default values
if (! "DiffInVEx_BW" %in% names(optInfo)) {
  optInfo <- c(optInfo, list(DiffInVEx_BW = 50))
}
#
if (! "reference_genome" %in% names(optInfo)) {
  optInfo <- c(optInfo, list(reference_genome = "hg19"))
}
#
BW              <- optInfo$DiffInVEx_BW
refGenome       <- optInfo$reference_genome
toolDirectory   <- optInfo$tool_directory
#
#
# toolDirectory <- "/home/akhalil/GDS/Studies/Tools/DiffInVEx_star3"
# BW <- 50
# refGenome <- "hg19"
#
####################
####################
load_InVEx_gene_data = function(HGNC_symbol,
                                toolDirectory,                                      
                                refGenome="hg19",
                                diffInVEx_BW=50){  
        geneDB_Dir <- paste0(toolDirectory,"/references/GeneData/geneDB_",refGenome,"_",diffInVEx_BW,"Kb")
        geneFile   <- paste0(geneDB_Dir,"/",HGNC_symbol,"_diffInVEx.RDS")
        if(file.exists(geneFile)){ 
                grDict  <- readRDS(geneFile)
                return(grDict)
        }
        else{
                cat(paste0(HGNC_symbol," files do not exist, please create them first \n"))
                return(NULL)
        }
}
#
# get exon/intron information
get_info <- function(geneList,toolDirectory,region, BW = 50){
        #
        geneFiles <- paste0(toolDirectory,"/references/GeneData/geneDB_hg19_",BW,"Kb")
        geneFiles <- list.files(geneFiles)
        geneFiles <- gsub("_diffInVEx.RDS","",geneFiles)
        geneList <- intersect(geneList,geneFiles)
        #
        exonInfo <- lapply(geneList, function(geneID){
              print(geneID)
              geneRDS  <- load_InVEx_gene_data(geneID, toolDirectory, refGenome="hg19", diffInVEx_BW=BW)
              #
              geneExon <- geneRDS[[region]]
              names(geneExon) <- NULL
              geneExon <- GenomicRanges::reduce(geneExon)      
              geneExon <- as.data.frame(geneExon)
              geneExon <- geneExon[,c("seqnames","start","end")]
              return(geneExon)})
        exonInfoDF <- do.call(rbind, exonInfo)
        exonInfoGR <- GRanges(exonInfoDF)
        exonInfoGR <- GenomicRanges::reduce(exonInfoGR)
        #
        return(exonInfoGR)
}
#
which_gene <- function(geneList, toolDirectory, region, BW = 50, iTracks){
        #
        geneFiles <- paste0(toolDirectory,"/references/GeneData/geneDB_hg19_",BW,"Kb")
        geneFiles <- list.files(geneFiles)
        geneFiles <- gsub("_diffInVEx.RDS","",geneFiles)
        geneList  <- intersect(geneList,geneFiles)
        #
        exonInfo <- lapply(geneList, function(geneID){
              geneRDS  <- load_InVEx_gene_data(geneID, toolDirectory, refGenome="hg19", diffInVEx_BW=BW)
              #
              geneExon <- geneRDS[[region]]
              iExon <- GenomicRanges::intersect(geneExon,iTracks, ignore.strand = TRUE)
              if(length(iExon) != 0){
                  print(geneID)
                  print(sum(width(iExon)))
                  print(sum(width(geneExon)))
                  return(geneID)
              } else{
                  return(NULL)
              }
              })
        #
        exonInfo <- do.call(rbind,exonInfo)
        return(exonInfo)
}
#####
#####
#####
#
#
#
geneDir   <- paste0(toolDirectory,"/references/GeneData/geneDB_hg19_",BW,"Kb")
geneFiles <- list.files(geneDir)
geneList  <- gsub("_diffInVEx.RDS","",geneFiles)
#
genomeExonicRegions    <- get_info(geneList,toolDirectory,"target")
genomeIntronicRegions_tri <- get_info(geneList,toolDirectory,"background_tri")
genomeIntronicRegions_penta  <- get_info(geneList,toolDirectory,"background_penta")
genomeIntronicRegions_NoMatching <- get_info(geneList,toolDirectory,"background_NoMatching")
#
#
totalGenomicRegions <- list()
totalGenomicRegions$target <- genomeExonicRegions
totalGenomicRegions$background_tri  <- genomeIntronicRegions_tri
totalGenomicRegions$background_penta   <- genomeIntronicRegions_penta
totalGenomicRegions$background_NoMatching <- genomeIntronicRegions_NoMatching
#
#
regionFile <- paste0(toolDirectory,"/references/GeneData/geneDB_hg19_",BW,"Kb","_GenomicRegions.RDS")
saveRDS(totalGenomicRegions,regionFile)
#
#
#
ExInIntersect <- GenomicRanges::intersect(genomeExonicRegions, genomeIntronicRegions_tri)
print(length(ExInIntersect))