#' Load Libraries and Process Genomic Data
#'
#' This R script loads required libraries, parses input arguments, and processes genomic data. It initializes 
#' necessary packages, handles input options for reference genome and background width, and performs operations 
#' to extract and save genomic regions data.
#'
#' @section Libraries:
#' The script loads the following libraries:
#' \itemize{
#'   \item `ddpcr`: For digital PCR analysis.
#'   \item `optparse`: For command-line option parsing.
#' }
#' Additionally, the script loads a set of additional packages:
#' \itemize{
#'   \item `arm`, `biomaRt`, `rtracklayer`, `GenomicRanges`, `GenomicFeatures`: For genomic data manipulation and analysis.
#'   \item `dplyr`, `reshape2`, `data.table`, `tidyverse`: For data manipulation and processing.
#' }
#'
#' @section Input Parsing:
#' The script accepts the following command-line options:
#' \itemize{
#'   \item `-r` or `--reference_genome`: Character string specifying the reference genome (default: "hg19").
#'   \item `-b` or `--DiffInVEx_BW`: Integer specifying the background width in kilobases (default: 20).
#'   \item `-t` or `--tool_directory`: Character string specifying the full path to the directory containing R and auxiliary files.
#' }
#'
#' @section Data Processing:
#' The script performs the following steps:
#' \itemize{
#'   \item Loads the necessary R libraries and suppresses messages.
#'   \item Parses input arguments to get the reference genome and background width.
#'   \item Sources an R script for retrieving InVEx gene data.
#'   \item Constructs file paths for gene data based on the provided background width.
#'   \item Retrieves genomic regions data for exonic and intronic regions using the `get_info` function.
#'   \item Compiles the retrieved data into a list and saves it to an RDS file.
#' }
#'
#' @section Outputs:
#' - `totalGenomicRegions`: A list containing genomic regions data:
#'   - `target`: Exonic regions.
#'   - `background_tri`: Intronic regions for background type "tri".
#'   - `background_penta`: Intronic regions for background type "penta".
#'   - `background_NoMatching`: Intronic regions for background type "NoMatching".
#' - The data is saved to a file named `geneDB_hg19_<BW>Kb_GenomicRegions.RDS` in the specified tool directory.
#'
#' @export
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
geneDir   <- paste0(toolDirectory,"/data/GeneData/geneDB_hg19_",BW,"Kb")
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
regionFile <- paste0(toolDirectory,"/data/GeneData/geneDB_hg19_",BW,"Kb","_GenomicRegions.RDS")
saveRDS(totalGenomicRegions,regionFile)