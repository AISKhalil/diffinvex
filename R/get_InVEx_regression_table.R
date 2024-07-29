get_InVEx_regression_table = function(HGNC_symbol,
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
  # temp   <- sample_wise_filtering(grSNPs, sample_annotation, grExons)
  # grSNPs <- temp[[1]]
  # sample_annotation <- temp[[2]] 
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