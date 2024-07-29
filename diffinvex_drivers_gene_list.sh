tool_directory=${1}
output_directory=${2}
regression_var=${3}
regression_type=${4}
#
#
DiffInVEx_topModule="$tool_directory"/R/get_diffinvex_drivers_gene_list.R
#
Rscript $DiffInVEx_topModule \
--output_directory $output_directory \
--regression_var $regression_var \
--tool_directory $tool_directory \
--regression_type $regression_type \