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

gene_symbol=${1}
cluster_number=${2}
mutations_data=${3}
genome_wide_sample_annotation_file=${4}
sample_annotation_perGene_file=${5}
variables_to_control_file=${6}
filtering_tracks=${7}
filter_inclusion=${8}
output_directory=${9}
tool_directory=${10}
reference_genome=${11}
neighbors_window=${12}
diffInVEx_mode=${13}
diffInVEx_BW=${14}
#
echo $gene_symbol $cluster_number $mutations_data $genome_wide_sample_annotation_file $sample_annotation_perGene_file $variables_to_control_file $output_directory $filtering_tracks $filter_inclusion $diffInVEx_mode $diffInVEx_BW
#
DiffInVEx_topModule=$tool_directory 
DiffInVEx_topModule+="/R/compute_diffInVEx_coefficients.R"
#
Rscript $DiffInVEx_topModule \
--gene_symbol $gene_symbol \
--cluster_number $cluster_number \
--mutations_data $mutations_data \
--genome_wide_sample_annotation_file $genome_wide_sample_annotation_file \
--sample_annotation_perGene_file $sample_annotation_perGene_file \
--variables_to_control_file $variables_to_control_file \
--filtering_tracks $filtering_tracks \
--filter_inclusion $filter_inclusion \
--output_directory $output_directory \
--tool_directory $tool_directory \
--reference_genome $reference_genome \
--neighbors_window $neighbors_window \
--diffInVEx_mode $diffInVEx_mode \
--diffInVEx_BW $diffInVEx_BW \