# Functions #
#############
cat("diffInVEx pipeline ...\n")
#
# loading data abd libraries #
suppressMessages(library("ddpcr"))
suppressMessages(library("optparse"))
packages <- c("arm","biomaRt","rtracklayer","GenomicRanges","GenomicFeatures", "plyr", "dplyr", "reshape", "data.table","tidyverse","readxl")
invisible(quiet(lapply(packages, library, character.only = TRUE), all = T)) # invisible for suppressing output of the lapply, suppresMessages for libraries
#
#####################
### input parsing ###
#####################
option_list <- list(
  make_option(c("-g", "--gene_symbol"), type="character", action="store", default="TP53",
              help="HGNC symbol for a gene"),
  make_option(c("-c", "--cluster_number"), type="integer", action="store", default=0,
              help="mutation representation:- 0: pentancleotide relative proportion matching, 1: trinucleotide relative proportion matching, others: MS96 clusterization"),
  make_option(c("-m", "--mutations_data"), type="character", action="store",
              help="Mutations file (csv) for input"),
  make_option(c("-a", "--genome_wide_sample_annotation_file"), type="character", action="store",
              help="A file, columns: sample, condition1, condition2 to control for. if conditions more than 1 there will be an interaction term included"),
  make_option(c("-b", "--sample_annotation_perGene_file"), type="character", action="store", default=NULL,
              help="A file, columns: gene, condition1, condition2 to control for. if conditions more than 1 there will be an interaction term included"),
  make_option(c("-v", "--variables_to_control_file"), type="character", action="store", default=F,
              help="File with variables (one per line) to supply to a regression"),
  make_option(c("-f", "--filtering_tracks"), type="character", action="store", default=NULL,
              help="filter mask to be applied on mutations, RDS file for a bed format"),
  make_option(c("-i", "--filter_inclusion"), type="character", action="store_true", default=F,
              help="True if only mutations withing the mask should be taken into account, False if otherwise (default)"),
  make_option(c("-o", "--output_directory"), type="character", action="store", default=F,
              help="Full path for the directory where output files should be stored"),  
  make_option(c("-t", "--tool_directory"), type="character", action="store", default=F,
              help="Full path for the diffInVEx directory with R and auxiliary files"),  
  make_option(c("-r", "--reference_genome"), type="character", action="store_true", default="hg19",
              help="Reference genome)"),
  make_option(c("-n", "--neighbors_window"), type="double", action="store", default=25000,
              help="length of genomic window for neighbors search"),
  make_option(c("-d", "--diffInVEx_mode"), type="integer", action="store", default=1,
              help="DiffInVEx mode: 1) use offline gene database (default), 2) compute & use gene information online (don't save), 3) compute,save & use gene information, 4) compute & save gene information"),
  make_option(c("-w", "--diffInVEx_BW"), type="integer", action="store", default=20,
              help="DiffInVEx background width (in Kb)")
)
###
###
###
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)
#
if (opt$filtering_tracks=="NULL") {opt$filtering_tracks = NULL}
if (! "filtering_tracks" %in% names(opt)) {
  opt <- c(opt, list(filtering_tracks = NULL))
}
#
if (opt$sample_annotation_perGene_file=="NULL") {opt$sample_annotation_perGene_file = NULL}
if (! "sample_annotation_perGene_file" %in% names(opt)) {
  opt <- c(opt, list(sample_annotation_perGene_file = NULL))
}
#
print(opt)
###
###
###
# loading DiffInVEx
quiet(source(paste0(opt$tool_directory,"/R/get_InVEx_selection_coefficients.R")), all = T)
quiet(source(paste0(opt$tool_directory,"/R/get_InVEx_gene_data.R")), all = T)
quiet(source(paste0(opt$tool_directory,"/R/get_InVEx_regression_table.R")), all = T)
quiet(source(paste0(opt$tool_directory,"/R/ParsingInputFiles.R")), all = T)
quiet(source(paste0(opt$tool_directory,"/R/MutationsToTable.R")), all = T) 
quiet(source(paste0(opt$tool_directory,"/R/RegressionsAndWriting.R")), all = T)
quiet(source(paste0(opt$tool_directory,"/R/glm.nb.safe_fast.R")), all = T)
###
###
###
# Computing coefficients
get_InVEx_selection_coefficients(HGNC_symbol = opt$gene_symbol,
                                clusterNumber = opt$cluster_number,
                                mutations = opt$mutations_data, 
                                annotation_genome_wide_file = opt$genome_wide_sample_annotation_file, 
                                sample_annotation_perGene_file = opt$sample_annotation_perGene_file, 
                                controlled_variables_input = opt$variables_to_control_file,
                                filterRegionsFile = opt$filtering_tracks, 
                                filterInclusion   = opt$filter_inclusion, 
                                outputDirectory = opt$output_directory,
                                toolDirectory = opt$tool_directory,
                                refGenome  = opt$reference_genome,
                                neighbors_window = opt$neighbors_window,
                                diffInVEx_mode = opt$diffInVEx_mode,
                                diffInVEx_BW = opt$diffInVEx_BW
                                )