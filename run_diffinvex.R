# DiffInvex needs 4 input files:
# 	- mutation_file: input .csv file with mutation data (single base substitutions and multi base substitutions) as follows:
# 							chrom,pos,ref,alt,Sample
# 							1,781002,G,A,PD1456			
# 	- annotation_file: input .csv file with annotations per sample (one row per sample) as follows:
# 							Sample,tumor_type,isTreated
# 							PD1456,NSCLC,1
# 	- variable_file: input .txt file with DiffInvex regression explantory variables 					
# 					 For a regression forumla #mutations ~ isTarget + isTreated + isTarget:isTreated + tumor_type,
# 					 each row should contain one regression term as follows:
# 					 		isTarget
# 					 		isTreated
# 					 		isTarget:isTreated
# 					 		tumor_type
# 	- gene_list: vector with gene names: c("KRAS","TP53","BRAF")

# Additionally, you need to provide the path to diffInvex and output directories:
# 	- output_directory: path for saving output files of the diffinvex method
# 	- diffinvex_directory: path to the diffInvex folder

# Finally, some DiffInvex parameters:
# 	- reference_genome: reference genome that is used for DiffInVEx (default "hg19")
# 	- no_cores: number of CPUs to be used. So, DiffInvex can be applied for many genes in parallel

# You can run this script on a data sample in the "input" directory inside the diffinvex_directory.
# Please modify this script based on your data
#
#
# 1) source diffinvex top modules
diffinvex_directory  <- "../diffinvex"
diffinvex_top_module <- paste0(diffinvex_directory,"/R/diffinvex.R")
source(diffinvex_top_module)
#
# 2) define your input files
mutation_file   <- paste0(diffinvex_directory,"/tests/pog570_SBSs.csv")
annotation_file <- paste0(diffinvex_directory,"/tests/pog570_annotations.csv")
variable_file   <- paste0(diffinvex_directory,"/tests/pog570_variables.txt") 
#
gene_file  <- paste0(diffinvex_directory,"/tests/cgc_tcga_genes.txt")
gene_list  <- read.table(file = gene_file, col.names = "Gene", sep = "")$Gene
#
# 3) set the output directory
output_directory <- paste0(diffinvex_directory,"/tests/pog570_diffinvex_results")
#
# 4) define DiffInvex parameters
reference_genome <- "hg19"
no_cores <- 4
#
# 5) run diffinvex
diffinvex(gene_list,
          mutation_file,
          annotation_file,
          variable_file,
          diffinvex_directory,
          output_directory,
          no_cores,
          reference_genome)