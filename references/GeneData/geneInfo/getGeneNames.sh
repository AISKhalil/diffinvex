# Genes information
https://www.genenames.org/download/custom/
#
# TREGT website, that ranks the transcripts per gene using GTEx expression data (detailed analyses in Tung et. al 2020).
# Briefly, they ranked the 145,571 transcripts for the 19,591 human protein-coding genes based on 1) Expression level,
# 2) Transcript or CDS length if the expression levels are identical for two or more transcripts.
# Notably, rank 1 transcripts of the 19,591 protein-coding genes are 17,472 protein-coding transcripts and 2,119 processed-transcript or retained-intron. 
# They suggested that these noncoding isoforms in protein-coding genes may possess additional modulation functions.
# 
# We create a list of the highest-rank protein-coding transcripts for each gene.
#
# Also, some genes such as AGPAT1 has many ensemble-id, is found on a haplotypic region. 
# Haplotypes are regions of the genome which have two or more versions, which we find in full in different individuals. These may have the same genes in a different order, or even different genes.
# We also selected highest-transcripted enemble-id (version) of each gene.
https://www.nature.com/articles/s41598-020-73081-5
https://tregt.ibms.sinica.edu.tw/index.php