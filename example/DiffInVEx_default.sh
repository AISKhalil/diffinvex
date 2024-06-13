mutFile=${1}
annFile=${2}
varFile=${3}
outDir=${4}
DiffInVExDir=${5}
##
##
echo $mutFile $annFile $varFile $outDir $DiffInVExDir
##
## DiffInVEx directory
DiffInVEx_Path=$DiffInVExDir
DiffInVEx_R="$DiffInVEx_Path"/R
export PATH="$PATH:$DiffInVEx_Path:$DiffInVEx_R"
## Protein-coding mutations
geneFile="$DiffInVEx_Path"/references/GeneData/geneInfo/TREG_18215Genes_goldenProteinCoding.txt
##
## Output directory
output_directory=$outDir
output_directory2="$output_directory"_figures
echo $output_directory
mkdir -p $output_directory
mkdir -p $output_directory2
#
logFile1="$output_directory2"/output_DiffInVEx_coefficients.log
logFile2="$output_directory2"/output_DiffInVEx_drivers.log
##
## DiffInVEx parameters
DiffInVEx_BW=50
DiffInVEx_cluster=1
DiffInVEx_mode=1
tool_directory=$DiffInVEx_Path
##
##
bash DiffInVEx_coefficients_geneList.sh $geneFile $mutFile $annFile $varFile $DiffInVEx_BW $DiffInVEx_cluster $DiffInVEx_mode $tool_directory $output_directory > $logFile1
#
regression_var=isTarget1
regression_type=bayes.poisson
bash DiffInVEx_candidate_drivers_geneList.sh $tool_directory $output_directory $regression_var $regression_type > $logFile2
#
find $output_directory -maxdepth 1 -type f -name "*.tsv" -delete