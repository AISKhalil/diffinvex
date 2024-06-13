get_InVEx_selection_coefficients = function(HGNC_symbol, 
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
    }   # else{
        # suppressMessages(TableToRegressAggregrated <- TableToRegress %>% group_by_at(setdiff(colnames(TableToRegress), c("Mutation","MutationNumber","ntAtRisk"))) %>%
        #                                                dplyr::summarize(MutationNumber = sum(MutationNumber)))
        #  print(TableToRegressAggregrated)
        # }
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
