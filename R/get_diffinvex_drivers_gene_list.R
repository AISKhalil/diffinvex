#' Run DiffInVEx Analysis and Generate Figures
#'
#' This R script executes the DiffInVEx (Differential Insertion/Deletion Variants and Exonic Variants) analysis and generates associated figures. It loads necessary libraries, parses command-line arguments, and runs the `diffinvex_drivers` function from a specified script.
#'
#' @section Libraries:
#' The script loads the following libraries:
#' \itemize{
#'   \item `optparse`: For parsing command-line options.
#'   \item `ddpcr`: A library required for the analysis.
#' }
#'
#' @section Command-line Arguments:
#' The script accepts the following command-line arguments:
#' \itemize{
#'   \item `-o` or `--output_directory`: The directory where the DiffInVEx output files will be stored.
#'   \item `-v` or `--regression_var`: The name of the regression variable used in the analysis.
#'   \item `-r` or `--regression_type`: The type of regression to be used in the analysis.
#'   \item `-t` or `--tool_directory`: The directory containing the DiffInVEx R script and auxiliary files.
#' }
#'
#' @section Script Execution:
#' \itemize{
#'   \item The script first creates a directory for figures, if it does not already exist.
#'   \item It then sources the `diffinvex_drivers.R` script located in the specified `tool_directory`.
#'   \item The `diffinvex_drivers` function is called with the parsed command-line arguments to perform the analysis and generate figures.
#' }
#'
#' @examples
#' # To run this script from the command line:
#' Rscript this_script.R --output_directory /path/to/output --regression_var var_name --regression_type type_name --tool_directory /path/to/tool_directory
#'
#' @export
suppressMessages(library("optparse"))
suppressMessages(library("ddpcr"))
##
## input arguments
option_list <- list(make_option(c("-o", "--output_directory"), type="character", action="store", help="DiffInVEx output directory"),
                    make_option(c("-v", "--regression_var"),   type="character", action="store", help="Regression variable"),
                    make_option(c("-r", "--regression_type"),   type="character", action="store", help="Regression type"),                    
                    make_option(c("-t", "--tool_directory"),   type="character", action="store", default=F, help="Full path for the diffInVEx directory with R and auxiliary files"))
##
parser <- OptionParser(option_list=option_list)
optInfo <- parse_args(parser)
#
output_directory <- optInfo$output_directory
regression_var   <- optInfo$regression_var
tool_directory   <- optInfo$tool_directory
regression_type  <- optInfo$regression_type 
##
##
figures_directory <- paste0(output_directory,"_figures")
if (!dir.exists(figures_directory)){
  dir.create(figures_directory)
}
##
##
quiet(source(paste0(tool_directory,"/R/diffinvex_drivers.R")), all = T)
diffinvex_drivers(output_directory,regression_var,regression_type)