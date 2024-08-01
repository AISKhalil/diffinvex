#' Compute InVEx Selection Coefficients
#'
#' This function calculates the selection coefficients for InVEx (Insertion/Deletion Variants and Exonic Variants) analysis. It retrieves gene-specific mutation data, performs regression analysis, and outputs selection coefficients.
#'
#' @param HGNC_symbol Character string. The HGNC symbol of the gene of interest.
#' @param clusterNumber Numeric. The cluster number used for background selection. Possible values are 0 (penta-matching), 1 (tri-matching), and other values for no-matching.
#' @param mutations `GRanges` object. Genomic ranges representing mutation data.
#' @param sample_annotation `data.frame`. A data frame containing sample annotations.
#' @param sample_annotation_perGene_file Character string. Path to a file with additional sample annotations specific to the gene. Default is `NULL`.
#' @param controlled_variables_input Character string. Path to a file containing controlled variables for regression analysis.
#' @param genomic_tracks `GRanges` object or other data structure. Genomic tracks used for retrieving gene data when `diffInVEx_mode` is not 1.
#' @param outputDirectory Character string. The path to the directory where results will be saved.
#' @param toolDirectory Character string. The path to the directory containing tool references.
#' @param refGenome Character string. Reference genome assembly. Default is "hg19".
#' @param diffInVEx_mode Numeric. Mode for retrieving InVEx gene data. Default is 1. Other values are handled by `get_InVEx_gene_data`.
#' @param diffInVEx_BW Numeric. Bin width for differential InVEx analysis. Default is 50.
#' @param neighbors_window Numeric. Window size for neighboring regions. Default is 25000.
#'
#' @return NULL. The function does not return a value but writes the selection coefficients to files in the specified output directory.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Retrieves gene-specific data based on the specified mode.
#'   \item Selects background data based on the `clusterNumber` parameter.
#'   \item Constructs a regression table using the `get_InVEx_regression_table` function.
#'   \item Computes selection coefficients using specified regression methods.
#'   \item Writes the results to files in the specified output directory.
#' }
#'
#' @export
get_InVEx_selection_coefficients <- function(HGNC_symbol,
                                            clusterNumber,
                                            mutations,
                                            sample_annotation,
                                            sample_annotation_perGene_file,
                                            controlled_variables_input,
                                            genomic_tracks,
                                            outputDirectory,
                                            toolDirectory,
                                            refGenome,
                                            diffInVEx_mode=1,
                                            diffInVEx_BW=50,
                                            neighbors_window=25000
                                            ){
    ##
    # Getting gene information (target vs background) #
    cat(paste0(HGNC_symbol,"...DiffInVEx gene data \n"))
    #
    if(diffInVEx_mode == 1){
                  geneData <- load_InVEx_gene_data(HGNC_symbol = HGNC_symbol,  
                                                 toolDirectory = toolDirectory,
                                                 refGenome     = refGenome,
                                                 diffInVEx_BW  = diffInVEx_BW)}
    else{
                  geneData <- get_InVEx_gene_data(HGNC_symbol     = HGNC_symbol, 
                                                genomic_tracks    = genomic_tracks,
                                                toolDirectory     = toolDirectory,
                                                refGenome         = refGenome,
                                                neighbors_window  = neighbors_window,
                                                outputDirectory   = outputDirectory,
                                                diffInVEx_mode    = diffInVEx_mode,
                                                diffInVEx_BW      = diffInVEx_BW)
    }
    ##
    ##
    # selct background based on the clusterNumber
    if(is.null(geneData) | (diffInVEx_mode == 4)){   
            return(NULL)
    } else{
            grExons   <- geneData[["target"]]
            #
            if(clusterNumber == 0){
                    # penta-matching
                    grIntrons <- geneData[["background_penta"]]
            } else if(clusterNumber == 1){
                    # tri-matching
                    grIntrons <- geneData[["background_tri"]]                
            } else{
                    # No-matching (tri-nucleotide)
                    grIntrons <- geneData[["background_NoMatching"]]                            
            }
            #
    }
    ##
    ##
    # getting regression table
    cat(paste0(HGNC_symbol,"...DiffInVEx mutation table \n"))
    #
    #
    TableToRegress <- get_InVEx_regression_table(HGNC_symbol = HGNC_symbol,
                                                 grExons = grExons,
                                                 grIntrons = grIntrons, 
                                                 clusterNumber = clusterNumber,      
                                                 mutations = mutations,
                                                 sample_annotation = sample_annotation,                                                
                                                 sample_annotation_perGene_file = sample_annotation_perGene_file,
                                                 controlled_variables_input = controlled_variables_input,
                                                 toolDirectory = toolDirectory,
                                                 outputDirectory = outputDirectory,
                                                 refGenome = refGenome)
    #
    if (is.null(TableToRegress)){   
        return(NULL)
    }
    ##
    ##
    # computing selection coefficients
    cat(paste0(HGNC_symbol,"...DiffInVEx selection coefficients \n"))
    # 
    regressions = c(# "bayes.negative.binomial", "negative.binomial",
                    # "zeroinfl.negative.binomial", "loglin.poisson", "loglin.bayes.poisson","poisson",
                    "bayes.poisson")
    #
    Regression_output = compute_regress_table(TableToRegress,
                                              clusterNumber = clusterNumber,
                                              regression_types = regressions,
                                              controlled_variables_input = controlled_variables_input,
                                              HGNC_symbol = HGNC_symbol)
    ##
    ##
    # outputing the coefficents     
    #writing_table(TableToRegress, outputDirectory, HGNC_symbol, Method = "diffInVEx")
    writing_files(Regression_output, outputDirectory, HGNC_symbol, Method = "diffInVEx")
}
###
###
###