#!/bin/bash
#
# SBATCH --job-name=dCdF
# Specify number of tasks, in general --ntasks should always be equal to 1 (unless you use MPI)
# SBATCH --ntasks=1
# Specify CPUs per task, if your application requires n threads, then set --cpus-per-tasks=n 
# (unless MPI is being used)
# Specify the memory required per node (1024MB by default)
# SBATCH --cpus-per-task=4
# module load R/3.5.0

gene_file=${1}
mutation_file=${2}
annotation_file=${3}
variable_file=${4}
DiffInVEx_BW=${5}
DiffInVEx_cluster=${6}
DiffInVEx_mode=${7}
tool_directory=${8}
output_directory=${9}
#
echo $gene_file $mutation_file $annotation_file $variable_file $DiffInVEx_BW $DiffInVEx_cluster $DiffInVEx_mode $tool_directory $output_directory
#
DiffInVEx_topModule1=$tool_directory 
DiffInVEx_topModule1+="/R/compute_diffInVEx_coefficients_geneList.R"
#
Rscript $DiffInVEx_topModule1 \
--gene_file $gene_file \
--mutation_file $mutation_file \
--annotation_file $annotation_file \
--variable_file $variable_file \
--DiffInVEx_BW $DiffInVEx_BW \
--DiffInVEx_cluster $DiffInVEx_cluster \
--DiffInVEx_mode $DiffInVEx_mode \
--tool_directory $tool_directory \
--output_directory $output_directory \
#
#
DiffInVEx_topModule2=$tool_directory 
DiffInVEx_topModule2+="/R/0_ExtractingAllTargetBackgroundRegions.R"
#
if [ $DiffInVEx_mode == 4 ]; then
	Rscript $DiffInVEx_topModule2 --DiffInVEx_BW $DiffInVEx_BW --tool_directory $tool_directory
fi

