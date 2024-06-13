DiffInVEx_Path = "/home/akhalil/GDS/Tools/DiffInVEx"
DiffInVEx_topModule <- paste0(DiffInVEx_Path,"/R/DiffInVEx_coefficients_geneList.R")
source(DiffInVEx_topModule)
#
geneFile = paste0(DiffInVEx_Path,"/references/GeneData/geneInfo/TREG_18215Genes_goldenProteinCoding.txt")
###
###
mutFile = "input/POG570_topAdundantTumors_SNVs.csv"
annFile = "input/POG570_topAdundantTumors_annotationFile.csv"
varFile = "input/POG570_variablesControlled_NoInteraction.txt"
#
output_directory = "output_POG570_example_noInteraction2"
###
###
DiffInVEx_BW = 50
DiffInVEx_cluster = 1
DiffInVEx_mode = 1
tool_directory = DiffInVEx_Path
#
regression_var = "isTarget1"
regression_type= "bayes.poisson"
###
###
DiffInVEx_selection_coefficients(gene_file          = geneFile,
                                  mutation_file     = mutFile,
                                  annotation_file   = annFile,
                                  variable_file     = varFile,
                                  DiffInVEx_BW      = DiffInVEx_BW,
                                  DiffInVEx_cluster = DiffInVEx_cluster,
                                  DiffInVEx_mode    = DiffInVEx_mode,
                                  tool_directory    = tool_directory,
                                  output_directory  = output_directory)
#
DiffInVEx_candidate_drivers(output_directory  = output_directory,
                            regression_var    = regression_var,
                            reg_type = regression_type)