suppressMessages(library("ddpcr"))
suppressMessages(library("optparse"))



#' Perform DiffInVEx Analysis with Specified Parameters
#'
#' This function performs a DiffInVEx (Differential Insertion/Deletion Variants and Exonic Variants) analysis.
#' It sources necessary modules and scripts, computes DiffInVEx regression coefficients, and identifies putative drivers
#' based on the specified parameters. The function integrates various scripts and utilities to handle genomic data
#' and perform analysis.
#'
#' @param gene_list Character vector. List of gene symbols to be analyzed.
#' @param mutation_file Character string. Path to the file containing mutation data.
#' @param annotation_file Character string. Path to the file containing gene annotation data.
#' @param variable_file Character string. Path to the file containing variable data.
#' @param tool_directory Character string. Path to the directory containing DiffInVEx R scripts and auxiliary files.
#' @param output_directory Character string. Directory where the output files will be saved.
#' @param no_cores Integer. Number of CPU cores to use (default is 1).
#' @param reference_genome Character string. Reference genome to be used (default is "hg19").
#' @param diffinvex_matching Integer. Parameter for choosing background matching 0: penta-matching (default), 1: tri-matching.
#' @param diffinvex_BW Integer. Background width in kilobases (default is 50).
#' @param diffinvex_mode Integer. Mode for gene data processing (default is 1).
#' @param regression_var Character string. Name of the regression variable (default is "isTarget1").
#' @param regression_type Character string. Type of regression model to be used (default is "bayes.poisson").
#' @param sample_annotation_perGene_file Character string. Path to the file containing gene-wise sample annotations (default is NULL).
#'
#' @section Process:
#' The function performs the following steps:
#' \itemize{
#'   \item Sources necessary R scripts from the specified `tool_directory` to load required modules and functions.
#'   \item Computes DiffInVEx regression coefficients using the `diffinvex_coefficients` function.
#'   \item Identifies putative drivers using the `diffinvex_drivers` function.
#' }
#'
#' @return None. The function performs analysis and saves results to the specified `output_directory`.
#'
#' @export
diffinvex <- function(gene_list,
                     mutation_file,
                     annotation_file,
                     variable_file,
                     tool_directory,
                     output_directory,
                     no_cores = 1,
                     reference_genome = "hg19",
                     diffinvex_matching = 0,
                     diffinvex_BW = 50,
                     diffinvex_mode = 1,
                     regression_var = "isTarget1",
                     regression_type = "bayes.poisson",
                     sample_annotation_perGene_file = NULL){
            #
            # diffinvex modules
            quiet(source(paste0(tool_directory,"/R/diffinvex_coefficients.R")), all = T)
            quiet(source(paste0(tool_directory,"/R/diffinvex_drivers.R")), all = T)
            quiet(source(paste0(tool_directory,"/R/get_InVEx_selection_coefficients.R")), all = T)
            quiet(source(paste0(tool_directory,"/R/get_InVEx_gene_data.R")), all = T)
            quiet(source(paste0(tool_directory,"/R/get_InVEx_regression_table.R")), all = T)
            quiet(source(paste0(tool_directory,"/R/ParsingInputFiles.R")), all = T)
            quiet(source(paste0(tool_directory,"/R/MutationsToTable.R")), all = T) 
            quiet(source(paste0(tool_directory,"/R/RegressionsAndWriting.R")), all = T)
            quiet(source(paste0(tool_directory,"/R/glm.nb.safe_fast.R")), all = T)            
            #
            # compute diffinvex regression coefficients
            diffinvex_coefficients(gene_list, mutation_file, annotation_file, variable_file, tool_directory, output_directory, 
                                no_cores, reference_genome, diffinvex_BW, diffinvex_matching, diffinvex_mode, sample_annotation_perGene_file)
            #
            # get diffinvex putative drivers
            diffinvex_drivers(output_directory,regression_var,regression_type)
            #
}