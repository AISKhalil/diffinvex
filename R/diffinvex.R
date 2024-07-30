###  loading data abd libraries ##
suppressMessages(library("ddpcr"))
suppressMessages(library("optparse"))
#
diffinvex = function(gene_list,
                     mutation_file,
                     annotation_file,
                     variable_file,
                     tool_directory,
                     output_directory,
                     no_cores= 1,
                     reference_genome = "hg19",
                     DiffInVEx_BW = 50,
                     DiffInVEx_cluster = 0,
                     DiffInVEx_mode = 1,
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
                                no_cores, reference_genome, DiffInVEx_BW, DiffInVEx_cluster, DiffInVEx_mode, sample_annotation_perGene_file)
            #
            # get diffinvex putative drivers
            diffinvex_drivers(output_directory,regression_var,regression_type)
            #
}