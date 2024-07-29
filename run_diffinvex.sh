: '
This bash script for running DiffInvex methods to identify selection and conditional selection in cancer.
DiffInvex needs 4 input files:
	- mutation_file: input .csv file with mutation data (single base substitutions and multi base substitutions) as follows:
							chrom,pos,ref,alt,Sample
							1,781002,G,A,PD1456			
	- annotation_file: input .csv file with annotations per sample (one row per sample) as follows:
							Sample,tumor_type,isTreated
							PD1456,NSCLC,1,0
	- variable_file: input .txt file with DiffInvex regression explantory variables 					
					 For a regression forumla #mutations ~ isTarget + isTreated + isTarget:isTreated + tumor_type,
					 each row should contain one regression term as follows:
					 		isTarget
					 		isTreated
					 		isTarget:isTreated
					 		tumor_type
	- gene_file: input .txt file with genes to be tested as follows
					"KRAS"					 		
					"TP53"
					"BRAF"

Additionally, you need to provide the path to diffInvex and output directories:
	- output_directory: path for saving output files of the diffinvex method
	- diffinvex_directory: path to the diffInvex folder

Finally, some DiffInvex parameters:
	- reference_genome: reference genome that is used for DiffInVEx (default "hg19")
	- no_cores: number of CPUs to be used. So, DiffInvex can be applied for many genes in parallel

You can test the diffinvex installation by running the given script on a pog570 data sample in the "example" directory inside the diffinvex_directory.
Please modify the script based on your data.
'

# 1) add DiffInvex directory to the environment variable $PATH
diffinvex_directory=../diffinvex
diffinvex_R="$diffinvex_directory"/R
export PATH="$PATH:$diffinvex_directory:$diffinvex_R"

# 2) define your input files
data_dir="./example"
mutation_file="$data_dir"/pog570_SBSs.csv
annotation_file="$data_dir"/pog570_annotations.csv
variable_file="$data_dir"/pog570_variables.txt
gene_file="$data_dir"/cgc_tcga_genes.txt

# 3) set the output directory
output_directory=./example/pog570_diffinvex_results

# 4) define DiffInvex parameters
reference_genome=hg19
no_cores=4

# 5) run the DiffInvex framework
bash diffinvex.sh $mutation_file $annotation_file $variable_file $gene_file $output_directory $diffinvex_directory $reference_genome $no_cores