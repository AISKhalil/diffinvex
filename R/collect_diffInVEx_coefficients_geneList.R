###  loading data and libraries
#BiocManager::install("qvalue")
suppressMessages(library("ddpcr"))
suppressMessages(library(optparse))
suppressMessages(library(qvalue))
suppressMessages(library("dplyr"))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
####################################
####       input argments      #####
####################################
##
option_list <- list(make_option(c("-o", "--output_directory"), type="character", action="store", help="DiffInVEx output directory"),
                    make_option(c("-v", "--regression_var"),   type="character", action="store", help="Regression variable"),
                    make_option(c("-r", "--regression_type"),   type="character", action="store", help="Regression type"),                    
                    make_option(c("-t", "--tool_directory"),   type="character", action="store", default=F, help="Full path for the diffInVEx directory with R and auxiliary files"))
##
parser <- OptionParser(option_list=option_list)
optInfo <- parse_args(parser)
#
output_directory <- optInfo$output_directory # "output_50Kb_tri_topAbundantTumors"
regression_var   <- optInfo$regression_var   # "isTarget1"
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
source(paste0(tool_directory,"/R/DiffInVEx_coefficients_geneList.R"))
DiffInVEx_candidate_drivers(output_directory,regression_var,regression_type)