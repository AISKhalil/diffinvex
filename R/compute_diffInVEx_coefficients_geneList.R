###  loading data abd libraries
suppressMessages(library("ddpcr"))
suppressMessages(library("optparse"))
#
packages <- c("arm","biomaRt","rtracklayer","GenomicRanges","GenomicFeatures", "dplyr", "reshape", "data.table","tidyverse","readxl")
invisible(quiet(lapply(packages, library, character.only = TRUE), all = T)) # invisible for suppressing output of the lapply, suppresMessages for libraries
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
  make_option(c("-b", "--DiffInVEx_BW"), type="integer", action="store", default=20,
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
reference_genome  <- optInfo$reference_genome
DiffInVEx_BW      <- optInfo$DiffInVEx_BW
DiffInVEx_cluster <- optInfo$DiffInVEx_cluster
DiffInVEx_mode    <- optInfo$DiffInVEx_mode
tool_directory    <- optInfo$tool_directory
output_directory  <- optInfo$output_directory
# 
print(gene_file)
print(mutation_file)
print(annotation_file)
print(sample_annotation_perGene_file)
print(variable_file)
print(reference_genome)
print(DiffInVEx_BW)
print(DiffInVEx_cluster)
print(DiffInVEx_mode)
print(tool_directory)
print(output_directory)
#
#######################
## DiffInVEx Module ###
#######################
# 
quiet(source(paste0(tool_directory,"/R/get_InVEx_selection_coefficients.R")), all = T)
quiet(source(paste0(tool_directory,"/R/get_InVEx_gene_data.R")), all = T)
quiet(source(paste0(tool_directory,"/R/get_InVEx_regression_table.R")), all = T)
quiet(source(paste0(tool_directory,"/R/ParsingInputFiles.R")), all = T)
quiet(source(paste0(tool_directory,"/R/MutationsToTable.R")), all = T) 
quiet(source(paste0(tool_directory,"/R/RegressionsAndWriting.R")), all = T)
quiet(source(paste0(tool_directory,"/R/glm.nb.safe_fast.R")), all = T)
#
#########################
##### Data loading  #####
#########################
print("Data loading ...")
start_time = Sys.time()
#-------#
#-Genes-#
#-------#
if(is_file(gene_file)){
    genesNames   <- read.table(file = gene_file, col.names = "Gene", sep = "")$Gene
    noGenes      <- length(genesNames)
} else{
    print("gene file does not exist")
    print("DiffInVEx will use the default TREG gene list instead")
    gene_file <- paste0(tool_directory,"/references/GeneData/geneInfo/TREG_18215Genes_goldenProteinCoding.txt")
    genesNames   <- read.table(file = gene_file, col.names = "Gene", sep = "")$Gene
    noGenes      <- length(genesNames)    
}
#
#--------------------------#
# target/background tracks #
#--------------------------#
if(DiffInVEx_mode != 4){
    trBg_tracks <- load_tracks(tool_directory,reference_genome,DiffInVEx_BW,DiffInVEx_cluster)
    trTracks <- trBg_tracks[["target"]]
    bgTracks <- trBg_tracks[["background"]]
    gTracks  <- trBg_tracks[["gene"]]
  } else{
    trTracks <- NULL
    bgTracks <- NULL
    gTracks  <- NULL
  }
#
#------#
#-SNVs-#
#------#
if(DiffInVEx_mode != 4){
    sample_annotation  <- load_sample_annotation(annotation_genome_wide_file = annotation_file)
    mutationData       <- MkGRanges_SNPs(mutation_file, gTracks, sample_annotation)
    sample_annotation  <- add_MutationLoad(mutationData,sample_annotation,variable_file,trTracks,bgTracks)
} else{
    sample_annotation <- NULL
    mutationData <- NULL
}
#
#----------------#
#-Genomic Tracks-#
#----------------#
filtering_tracks <- NULL
filter_inclusion <- TRUE
#
# read genomic tracks if exonic/intronic regions will be computed
if(DiffInVEx_mode != 1){
    quiet(source(paste0(tool_directory,"/R/GettingGenomicFeatures.R")), all = T)
    genomic_tracks <- read_genomic_tracks(reference_genome, tool_directory, filtering_tracks, filter_inclusion)
} else{
    genomic_tracks <- NULL
}
#
end_time = Sys.time()
end_time - start_time
#
#############
# variables #
#############
opt <- list()
opt$mutations_data                     <- mutationData
opt$sample_annotation                  <- sample_annotation
opt$variables_to_control_file          <- variable_file
opt$sample_annotation_perGene_file     <- sample_annotation_perGene_file
opt$genomic_tracks                     <- genomic_tracks
#
opt$tool_directory   <- tool_directory
opt$output_directory <- output_directory
#
opt$cluster_number   <- DiffInVEx_cluster
opt$reference_genome <- reference_genome
opt$diffInVEx_mode   <- DiffInVEx_mode
opt$diffInVEx_BW     <- DiffInVEx_BW
###
###
if (!dir.exists(opt$output_directory)){
  dir.create(opt$output_directory)
}
######################################
### compute selection coefficients ###
######################################
print("processing ...")
start_time = Sys.time()
##
##-Parallel code-##
require(parallel)
noCores <- detectCores()
noCoresToUse <- max(as.integer(noCores - 2),1)
#
mclapply(genesNames, 
         function(gene) get_InVEx_selection_coefficients(HGNC_symbol = gene,
                                                        clusterNumber = opt$cluster_number,
                                                        mutations = opt$mutations_data, 
                                                        sample_annotation = opt$sample_annotation, 
                                                        sample_annotation_perGene_file = opt$sample_annotation_perGene_file, 
                                                        controlled_variables_input = opt$variables_to_control_file,
                                                        genomic_tracks = opt$genomic_tracks,  
                                                        outputDirectory = opt$output_directory,
                                                        toolDirectory = opt$tool_directory,
                                                        refGenome  = opt$reference_genome,
                                                        diffInVEx_mode = opt$diffInVEx_mode,
                                                        diffInVEx_BW = opt$diffInVEx_BW),
          mc.cores = noCoresToUse)
##
##
end_time = Sys.time()
end_time - start_time
###
###
#genesNames <- c("AR","KRAS","NRAS","TP53","KDM6A")
#for(gene in genesNames){
#                       get_InVEx_selection_coefficients(HGNC_symbol = gene,
#                                                        clusterNumber = opt$cluster_number,
#                                                        mutations = opt$mutations_data, 
#                                                        sample_annotation = opt$sample_annotation, 
#                                                        sample_annotation_perGene_file = opt$sample_annotation_perGene_file, 
#                                                        controlled_variables_input = opt$variables_to_control_file,
#                                                        genomic_tracks = opt$genomic_tracks,  
#                                                        outputDirectory = opt$output_directory,
#                                                        toolDirectory = opt$tool_directory,
#                                                        refGenome  = opt$reference_genome,
#                                                        diffInVEx_mode = opt$diffInVEx_mode,
#                                                        diffInVEx_BW = opt$diffInVEx_BW)
# }
#
