###  loading data abd libraries ##
suppressMessages(library("ddpcr"))
suppressMessages(library("optparse"))
#
packages <- c("arm","biomaRt","rtracklayer","GenomicRanges","GenomicFeatures", "dplyr", "reshape2", "data.table","tidyverse",
               "readxl","parallel","qvalue","Matrix","ggplot2",
               "MASS","stringr","fastglm","pscl")
#
invisible(quiet(lapply(packages, library, character.only = TRUE), all = T)) # invisible for suppressing output of the lapply, suppresMessages for libraries
###
###
###-------------------------------------###
diffinvex_coefficients = function(gene_list,
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
                                   sample_annotation_perGene_file = NULL){
        #--------------------------#
        # target/background tracks #
        #--------------------------#
        print("Data loading ...")
        start_time = Sys.time()
        #
        if(DiffInVEx_mode != 4){
            # use pre-computed target and background tracks
            trBg_tracks <- load_tracks(tool_directory,reference_genome,DiffInVEx_BW,DiffInVEx_cluster)
            trTracks <- trBg_tracks[["target"]]
            bgTracks <- trBg_tracks[["background"]]
            gTracks  <- trBg_tracks[["gene"]]
          } else{
            trTracks <- NULL
            bgTracks <- NULL
            gTracks  <- NULL
          }
        #
        #------#
        #-SNVs-#
        #------#
        if(DiffInVEx_mode != 4){
            sample_annotation  <- load_sample_annotation(annotation_genome_wide_file = annotation_file)
            mutationData       <- MkGRanges_SNPs(mutation_file, gTracks, sample_annotation)
            sample_annotation  <- add_MutationLoad(mutationData,sample_annotation,variable_file,trTracks,bgTracks)
        } else{
            sample_annotation <- NULL
            mutationData <- NULL
        }
        #
        #----------------#
        #-Genomic Tracks-#
        #----------------#
        filtering_tracks <- NULL
        filter_inclusion <- TRUE
        #
        # read genomic tracks if exonic/intronic regions will be computed
        if(DiffInVEx_mode != 1){
            quiet(source(paste0(tool_directory,"/R/GettingGenomicFeatures.R")), all = T)
            genomic_tracks <- read_genomic_tracks(reference_genome, tool_directory, filtering_tracks, filter_inclusion)
        } else{
            genomic_tracks <- NULL
        }
        #
        end_time = Sys.time()
        end_time - start_time
        #
        #############
        # variables #
        #############
        opt <- list()
        opt$mutations_data                     <- mutationData
        opt$sample_annotation                  <- sample_annotation
        opt$variables_to_control_file          <- variable_file
        opt$sample_annotation_perGene_file     <- sample_annotation_perGene_file
        opt$genomic_tracks                     <- genomic_tracks
        #
        opt$tool_directory   <- tool_directory
        opt$output_directory <- output_directory
        #
        opt$cluster_number   <- DiffInVEx_cluster
        opt$reference_genome <- reference_genome
        opt$diffInVEx_mode   <- DiffInVEx_mode
        opt$diffInVEx_BW     <- DiffInVEx_BW
        #
        if (!dir.exists(opt$output_directory)){
          dir.create(opt$output_directory)
        }
        #
        ######################################
        ### compute selection coefficients ###
        ######################################
        print("processing ...")
        start_time = Sys.time()
        ##
        ##-Parallel code-##
        #
        mclapply(gene_list, 
                 function(gene) get_InVEx_selection_coefficients(HGNC_symbol = gene,
                                                                clusterNumber = opt$cluster_number,
                                                                mutations = opt$mutations_data, 
                                                                sample_annotation = opt$sample_annotation, 
                                                                sample_annotation_perGene_file = opt$sample_annotation_perGene_file, 
                                                                controlled_variables_input = opt$variables_to_control_file,
                                                                genomic_tracks = opt$genomic_tracks,  
                                                                outputDirectory = opt$output_directory,
                                                                toolDirectory = opt$tool_directory,
                                                                refGenome  = opt$reference_genome,
                                                                diffInVEx_mode = opt$diffInVEx_mode,
                                                                diffInVEx_BW = opt$diffInVEx_BW),
                  mc.cores = no_cores)
        #
        end_time = Sys.time()
        end_time - start_time
}