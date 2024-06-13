### Different methods for installing libraries  ###
###################################################
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
#
# sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("devtools")
library(devtools)
#
## Examples 
# install.packages("tidyverse")
# BiocManager::install("pcaMethods")
# install_github("velocyto-team/velocyto.R")
# remotes::install_github("mojaveazure/seurat-disk")
# remotes::install_github('satijalab/seurat-wrappers')
######################################################
# diffInVEx libraries
install.packages(c("ddpcr","optparse","arm","dplyr","plyr","reshape","data.table","tidyverse","readxl","stringr","MASS","fastglm","pscl"))
BiocManager::install(c("GenomicRanges","GenomicFeatures", "BSgenome","ensembldb","Biostrings","biomaRt","rtracklayer"))
#
## ref-genome + annotation database (hg19)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
##
# ref-genome + annotation database (hg38)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
