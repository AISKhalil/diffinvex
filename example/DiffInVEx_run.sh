DiffInVExDir=/home/akhalil/GDS/Tools/DiffInVEx
###
###
# mutFile: input file with SNVs as follows
# 		chrom,pos,ref,alt,Sample
# 		1,781002,G,A,CPCT02010003T
#
# annFile: samples annotation as follows
# 		Sample,tumorType,isMetastatic,isTreated
# 		CPCT02010035T,NSCLC,1,0
#
# varFile: variables to control for as follows
# 		isTarget
#		isTreated
#		isTarget:isTreated
#		tumorType
#		isMetastatic
###
###
mutFile=input/POG570_topAdundantTumors_SNVs.csv #mutation file
annFile=input/POG570_topAdundantTumors_annotationFile.csv
varFile=input/POG570_variablesControlled_NoInteraction.txt
#
outputDirectory=output_POG570_example_noInteraction 

bash DiffInVEx_default.sh $mutFile $annFile $varFile $outputDirectory $DiffInVExDir
###
###
