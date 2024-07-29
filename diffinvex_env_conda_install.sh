# 1) create conda environment 
# this installs r-data.table, r-dplyr, r-ggplot2, r-reshape2, r-tidyverse, r-readxl, r-parallelly, r-matrix, r-ggplot2, r-mass and r-stringr
conda create -n diffinvex_env r-base r-essentials
#
# 2) activate the conda environment
conda activate diffinvex_env
#
# 3) install R packages using conda install
conda install r-optparse r-ddpcr r-arm r-tidyverse r-fastglm r-pscl
conda install bioconda::bioconductor-biomart bioconda::bioconductor-rtracklayer bioconda::bioconductor-genomicranges bioconda::bioconductor-genomicfeatures bioconda::bioconductor-qvalue
#
# 4) export conda enviroment
conda env export > diffinvex_env.yml