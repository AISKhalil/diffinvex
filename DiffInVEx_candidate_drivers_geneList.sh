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

tool_directory=${1}
output_directory=${2}
regression_var=${3}
regression_type=${4}
#
echo $tool_directory $output_directory $regression_var
#
#
DiffInVEx_topModule2=$tool_directory 
DiffInVEx_topModule2+="/R/collect_diffInVEx_coefficients_geneList.R"
#
Rscript $DiffInVEx_topModule2 \
--output_directory $output_directory \
--regression_var $regression_var \
--tool_directory $tool_directory \
--regression_type $regression_type \