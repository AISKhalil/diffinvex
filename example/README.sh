Here, we provided two scripts (DiffInVEx_run.sh and DiffInVEx_run.R) for running DiffInVEx computational method 
using its default parameters.
#
First, please install the R libraries in rPackages_installation.R
#
Then, you can use these scripts on your own data by providing three main files:
mutFile: input file with SNVs in the following format: chrom,pos,ref,alt,Sample (one mutation per line).
annFile: - annotation of samples. For example, it can contain tumor type and treated/untreated for each sample (one sample per line)
         in the following format: Sample,tumorType,isTreated.
         - it should include at least one annotation per sample (even if all samples have the same tumor type or cohort).
varFile: list of variables to control in the regression in addition to "isTarget" and "Mutation".
