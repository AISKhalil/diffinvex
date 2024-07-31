# The diffinvex script for computing selection coefficients and driver genes.
# User need to first activate the diffinvex_env conda environment.
# Then, run this provided scripts using the following 8 arguments.
mutFile=${1}
annFile=${2}
varFile=${3}
geneFile=${4}
outDir=${5}
diffinvexDir=${6}
refGenome=${7}
noCores=${8}
#
echo $mutFile $annFile $varFile $geneFile $outDir $diffinvexDir
#
#-- creating output directory --#
output_directory=$outDir
output_directory2="$output_directory"_figures
mkdir -p $output_directory
mkdir -p $output_directory2
#
logFile1="$output_directory2"/output_diffinvex_coefficients.log
logFile2="$output_directory2"/output_diffinvex_drivers.log
#
#-- DiffInvex parameters --#
diffinvex_BW=50
diffinvex_cluster=0
diffinvex_mode=1
diffinvex_directory=$diffinvexDir
reference_genome=$refGenome
no_cores=$noCores
regression_var=isTarget1
regression_type=bayes.poisson
#
#-- computing selection coefficients per gene --#
bash diffinvex_coefficients_gene_list.sh $mutFile $annFile $varFile $geneFile $output_directory $diffinvex_directory $reference_genome $no_cores $diffinvex_BW $diffinvex_cluster $diffinvex_mode > $logFile1
#
#-- collecting DiffInvex results --#
bash diffinvex_drivers_gene_list.sh $diffinvex_directory $output_directory $regression_var $regression_type > $logFile2
# 
#-- deleting DiffInvex temporary files --#
# find $output_directory -maxdepth 1 -type f -name "*.tsv" -delete