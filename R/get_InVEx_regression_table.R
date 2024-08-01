#' Generate Regression Table for InVEx Analysis
#'
#' This function constructs a regression table for InVEx (Insertion/Deletion Variants and Exonic Variants) analysis. It processes mutation data, filters samples, and prepares a regression table for both target gene exons and neighboring introns. The output is suitable for statistical modeling and analysis.
#'
#' @param HGNC_symbol Character string. The HGNC symbol of the gene of interest.
#' @param grExons `GRanges` object. Genomic ranges representing exonic regions of the target gene.
#' @param grIntrons `GRanges` object. Genomic ranges representing intronic regions of the target gene.
#' @param clusterNumber Numeric. The cluster number to use for filtering samples. Default is 1.
#' @param mutations `GRanges` object. Genomic ranges representing mutation data.
#' @param sample_annotation `data.frame`. A data frame containing sample annotations.
#' @param sample_annotation_perGene_file Character string. Path to a file containing additional sample annotations specific to the gene. Default is `NULL`.
#' @param controlled_variables_input Character string. Path to a file containing controlled variables for the regression analysis.
#' @param toolDirectory Character string. The path to the directory containing tool references.
#' @param outputDirectory Character string. The path to the directory where results will be saved.
#' @param refGenome Character string. Reference genome assembly. Default is "hg19".
#'
#' @return A `data.frame` containing the regression table for mutations in exons and introns. Returns `NULL` if no mutations are found within the specified regions.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Filters the mutation data to include only those located on the same chromosome as the exonic regions.
#'   \item Applies sample-wise filtering based on the mutation data and sample annotations.
#'   \item Loads additional sample annotations if a specific file for the gene is provided.
#'   \item Writes mutation data to a file for further analysis.
#'   \item Reads controlled variables from a specified input file.
#'   \item Prepares regression tables for both exonic and intronic regions using the `Mutations_toRegress` function.
#'   \item Combines these tables, updates them with relevant variables, and prepares the final regression table.
#'   \item Returns the final regression table or `NULL` if no relevant data is present.
#' }
#'
#' @export
get_InVEx_regression_table <- function(HGNC_symbol,
                                      grExons,
                                      grIntrons,
                                      clusterNumber=1,
                                      mutations,
                                      sample_annotation,
                                      sample_annotation_perGene_file,
                                      controlled_variables_input,
                                      toolDirectory,
                                      outputDirectory,
                                      refGenome="hg19"){
  ###########################
  ####  loading mutations ###
  ###########################
  #
  chr <- as.character(grExons@seqnames[1])
  grSNPs <- mutations[mutations@seqnames == chr]
  #
  temp   <- sample_wise_filtering(grSNPs, sample_annotation, grExons)
  grSNPs <- temp[[1]]
  sample_annotation <- temp[[2]] 
  #
  if(! is.null(sample_annotation_perGene_file)){
      sample_annotation <- load_gene_annotation(sample_annotation, sample_annotation_perGene_file, HGNC_symbol)
  }
  #
  SNPsInGene <- writing_mutations(outputDirectory, HGNC_symbol, grSNPs, grExons, grIntrons, Method = "diffInVEx")
  if(!SNPsInGene){
        return(NULL)
  }
  #################################
  ####  making regression table ###
  #################################
  #
  colNamesToRegress = fread(controlled_variables_input, header = F)
  colNamesToRegress = colNamesToRegress$V1  
  #  
  if("MutationLoad" %in% colNamesToRegress){
            sample_annotation1 <- sample_annotation %>% dplyr::select(setdiff(colnames(sample_annotation),"MutationLoad_Bg"))
            names(sample_annotation1)[names(sample_annotation1) == "MutationLoad_Tr"] <- "MutationLoad"
            #
            sample_annotation2 <- sample_annotation %>% dplyr::select(setdiff(colnames(sample_annotation),"MutationLoad_Tr"))
            names(sample_annotation2)[names(sample_annotation2) == "MutationLoad_Bg"] <- "MutationLoad"
  } else {
            sample_annotation1 <- sample_annotation
            sample_annotation2 <- sample_annotation
  }
  #
  #
  TargetGeneTABLEforGLM <- Mutations_toRegress(grObject_within = grExons, 
                                               grMutations = grSNPs,
                                               HGNC_symbol = HGNC_symbol,
                                               annotation  = sample_annotation1,
                                               clusterNumber = clusterNumber)
  #
  NeighboursTABLEforGLM <- Mutations_toRegress(grObject_within = grIntrons, 
                                              grMutations = grSNPs,
                                              HGNC_symbol = HGNC_symbol,
                                              annotation  = sample_annotation2,
                                              clusterNumber=clusterNumber)
  #
  # Adding isTarget variable (grExons: 1, grIntrons: 0)
  TargetGeneTABLEforGLM$isTarget = factor(1)
  NeighboursTABLEforGLM$isTarget = factor(0)
  # binding the exon and intron tables
  TableToRegress = rbind(TargetGeneTABLEforGLM,NeighboursTABLEforGLM)
  TableToRegress$isTarget = relevel(TableToRegress$isTarget, "0")
  #
  # Regression table update (nts at risk)
  TableToRegress[, ntAtRisk := Counts*noSamples]
  TableToRegress$Context    <- NULL
  TableToRegress$Counts     <- NULL
  TableToRegress$noSamples  <- NULL
  ###
  ###
  ###
  if (nrow(TableToRegress) == 0) {
    TableToRegress = NULL
  }
  #
  return(TableToRegress)
}