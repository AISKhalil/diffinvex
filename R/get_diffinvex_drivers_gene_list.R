###  loading data abd libraries ###
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