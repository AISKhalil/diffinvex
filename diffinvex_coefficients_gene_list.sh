mutation_file=${1}
annotation_file=${2}
variable_file=${3}
gene_file=${4}
output_directory=${5}
tool_directory=${6}
reference_genome=${7}
no_cores=${8}
DiffInVEx_BW=${9}
DiffInVEx_cluster=${10}
DiffInVEx_mode=${11}
#
#-- running DiffInvex for computing selection coefficients --#
DiffInVEx_topModule1="$tool_directory"/R/compute_diffInVEx_coefficients_gene_list.R
#
Rscript $DiffInVEx_topModule1 \
--gene_file $gene_file \
--mutation_file $mutation_file \
--annotation_file $annotation_file \
--variable_file $variable_file \
--reference_genome $reference_genome \
--no_cores $no_cores \
--DiffInVEx_BW $DiffInVEx_BW \
--DiffInVEx_cluster $DiffInVEx_cluster \
--DiffInVEx_mode $DiffInVEx_mode \
--tool_directory $tool_directory \
--output_directory $output_directory \
#
#
DiffInVEx_topModule2="$tool_directory"/R/extract_all_target_background_regions.R
#
if [ $DiffInVEx_mode == 4 ]; then
	Rscript $DiffInVEx_topModule2 --DiffInVEx_BW $DiffInVEx_BW --reference_genome $reference_genome --tool_directory $tool_directory
fi

