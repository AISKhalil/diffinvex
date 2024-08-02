#' Run DiffInVEx Analysis with Custom Parameters
#'
#' This R script performs the DiffInVEx (Differential Insertion/Deletion Variants and Exonic Variants) analysis based on 
#' user-provided parameters. It loads necessary libraries, parses command-line arguments, and executes the analysis functions. 
#' This script includes functionality to handle custom gene lists, mutation files, and annotation data.
#'
#' @section Libraries:
#' The script loads the following libraries:
#' \itemize{
#'   \item `optparse`: For parsing command-line options.
#'   \item `ddpcr`: Required for performing the analysis.
#' }
#'
#' @section Command-line Arguments:
#' The script accepts the following command-line arguments:
#' \itemize{
#'   \item `-g` or `--gene_file`: Path to the gene list file.
#'   \item `-m` or `--mutation_file`: Path to the mutation file.
#'   \item `-a` or `--annotation_file`: Path to the annotation file.
#'   \item `-v` or `--variable_file`: Path to the variable file.
#'   \item `-s` or `--sample_annotation_perGene_file`: Path to the gene-wise annotation file.
#'   \item `-r` or `--reference_genome`: Reference genome (default is "hg19").
#'   \item `-n` or `--no_cores`: Number of CPU cores to use (default is 1).
#'   \item `-b` or `--DiffInVEx_BW`: Background width in Kb (default is 50).
#'   \item `-c` or `--DiffInVEx_cluster`: Cluster number for background matching (default is 1).
#'   \item `-d` or `--DiffInVEx_mode`: Mode for gene data processing (default is 1).
#'   \item `-t` or `--tool_directory`: Full path to the directory containing DiffInVEx R scripts and auxiliary files.
#'   \item `-o` or `--output_directory`: Directory where the output files will be saved.
#' }
#'
#' @section Script Execution:
#' \itemize{
#'   \item The script starts by creating a directory for figures if it does not exist.
#'   \item It then sources various R scripts from the `tool_directory`, including functions for computing selection coefficients, handling gene data, and generating regression tables.
#'   \item If the specified gene file exists, it is read and used; otherwise, a default gene list is used.
#'   \item The main function `diffinvex_coefficients` is called with the parsed arguments to perform the analysis and generate output files.
#' }
#'
#' @examples
#' # To run this script from the command line:
#' Rscript this_script.R --gene_file /path/to/gene_file.txt --mutation_file /path/to/mutation_file.txt --annotation_file /path/to/annotation_file.txt --variable_file /path/to/variable_file.txt --sample_annotation_perGene_file /path/to/sample_annotation_perGene_file.txt --reference_genome hg19 --no_cores 4 --DiffInVEx_BW 50 --DiffInVEx_cluster 1 --DiffInVEx_mode 1 --tool_directory /path/to/tool_directory --output_directory /path/to/output_directory
#'
#' @export
suppressMessages(library("optparse"))
suppressMessages(library("ddpcr"))
#
#####################
### input parsing ###
#####################
option_list <- list(
  make_option(c("-g", "--gene_file"), type="character", action="store", default=NULL,
              help="gene list"),
  make_option(c("-m", "--mutation_file"), type="character", action="store", default=NULL,
              help="mutation file"),
  make_option(c("-a", "--annotation_file"), type="character", action="store", default=NULL,
              help="annotation file"),
  make_option(c("-v", "--variable_file"), type="character", action="store", default=NULL,
              help="variable file"),
  make_option(c("-s", "--sample_annotation_perGene_file"), type="character", action="store", default=NULL,
              help="gene-wise annotation file"),  
  make_option(c("-r", "--reference_genome"), type="character", action="store", default="hg19",
              help="reference genome"),  
  make_option(c("-n", "--no_cores"), type="integer", action="store", default=1,
              help="number of CPUs"),
  make_option(c("-b", "--DiffInVEx_BW"), type="integer", action="store", default=50,
              help="background width (in Kb)"),
  make_option(c("-c", "--DiffInVEx_cluster"), type="integer", action="store", default=1,
              help="cluster number: 0 (penta-matching), 1(tri-matching), 96 (tri no-matching)"),
  make_option(c("-d", "--DiffInVEx_mode"), type="integer", action="store", default=1,
              help="mode: 1) use offline gene database (default), 2) compute & use gene information online (don't save), 3) compute,save & use gene information, 4) compute & save gene information"),  
  make_option(c("-t", "--tool_directory"), type="character", action="store", default=F,
              help="Full path for the diffInVEx directory with R and auxiliary files"),   
  make_option(c("-o", "--output_directory"), type="character", action="store", default=NULL,
              help="output directory")
)
#
#
parser <- OptionParser(option_list=option_list)
optInfo <- parse_args(parser)
#
#------------- default values ------------#
if (! "reference_genome" %in% names(optInfo)) {
  optInfo <- c(optInfo, list(reference_genome = "hg19"))
}
#
if (! "DiffInVEx_BW" %in% names(optInfo)) {
  optInfo <- c(optInfo, list(DiffInVEx_BW = 50))
}
#
if (! "no_cores" %in% names(optInfo)) {
  optInfo <- c(optInfo, list(no_cores = 1))
}
#
if (! "DiffInVEx_cluster" %in% names(optInfo)) {
  optInfo <- c(optInfo, list(DiffInVEx_cluster = 0))
}
#
if (! "DiffInVEx_mode" %in% names(optInfo)) {
  optInfo <- c(optInfo, list(DiffInVEx_mode = 1))
}
#
if (! "reference_genome" %in% names(optInfo)) {
  optInfo <- c(optInfo, list(reference_genome = "hg19"))
}
#
gene_file         <- optInfo$gene_file
mutation_file     <- optInfo$mutation_file
annotation_file   <- optInfo$annotation_file
sample_annotation_perGene_file <- optInfo$sample_annotation_perGene_file
variable_file     <- optInfo$variable_file
no_cores          <- optInfo$no_cores
reference_genome  <- optInfo$reference_genome
DiffInVEx_BW      <- optInfo$DiffInVEx_BW
DiffInVEx_cluster <- optInfo$DiffInVEx_cluster
DiffInVEx_mode    <- optInfo$DiffInVEx_mode
tool_directory    <- optInfo$tool_directory
output_directory  <- optInfo$output_directory
# 
#######################
## DiffInVEx Module ###
#######################
# 
quiet(source(paste0(tool_directory,"/R/diffinvex_coefficients.R")), all = T)
quiet(source(paste0(tool_directory,"/R/get_InVEx_selection_coefficients.R")), all = T)
quiet(source(paste0(tool_directory,"/R/get_InVEx_gene_data.R")), all = T)
quiet(source(paste0(tool_directory,"/R/get_InVEx_regression_table.R")), all = T)
quiet(source(paste0(tool_directory,"/R/ParsingInputFiles.R")), all = T)
quiet(source(paste0(tool_directory,"/R/MutationsToTable.R")), all = T) 
quiet(source(paste0(tool_directory,"/R/RegressionsAndWriting.R")), all = T)
quiet(source(paste0(tool_directory,"/R/glm.nb.safe_fast.R")), all = T)
#
#-------#
#-Genes-#
#-------#
if(is_file(gene_file)){
    gene_list  <- read.table(file = gene_file, col.names = "Gene", sep = "")$Gene
} else{
    print("gene file does not exist")
    print("DiffInVEx will use the default TREG gene list instead")
    gene_file  <- paste0(tool_directory,"/data/GeneData/geneInfo/TREG_18215Genes_goldenProteinCoding.txt")
    gene_list <- read.table(file = gene_file, col.names = "Gene", sep = "")$Gene
}
#
############################
## DiffInvex coefficients ##
############################
#
diffinvex_coefficients(gene_list, mutation_file, annotation_file, variable_file, tool_directory, output_directory, 
                      no_cores, reference_genome, DiffInVEx_BW, DiffInVEx_cluster, DiffInVEx_mode, sample_annotation_perGene_file)