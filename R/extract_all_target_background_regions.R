###  loading data abd libraries ##
suppressMessages(library("ddpcr"))
suppressMessages(library("optparse"))
#
packages <- c("arm","biomaRt","rtracklayer","GenomicRanges","GenomicFeatures", "dplyr", "reshape2", "data.table","tidyverse")
#
invisible(quiet(lapply(packages, library, character.only = TRUE), all = T)) # invisible for suppressing output of the lapply, suppresMessages for libraries
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
###
###
###
#
quiet(source(paste0(toolDirectory,"/R/get_InVEx_gene_data.R")), all = T)
#
#
geneDir   <- paste0(toolDirectory,"/references/GeneData/geneDB_hg19_",BW,"Kb")
geneFiles <- list.files(geneDir)
geneList  <- gsub("_diffInVEx.RDS","",geneFiles)
#
genomeExonicRegions    <- get_info(geneList, toolDirectory, "target", BW)
genomeIntronicRegions_tri <- get_info(geneList, toolDirectory, "background_tri", BW)
genomeIntronicRegions_penta  <- get_info(geneList, toolDirectory, "background_penta", BW)
genomeIntronicRegions_NoMatching <- get_info(geneList, toolDirectory, "background_NoMatching", BW)
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