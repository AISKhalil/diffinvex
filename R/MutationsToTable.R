#' Compute One-Nucleotide Context Composition
#'
#' This function calculates the composition of nucleotide contexts (A, T, G, C) within the given genomic regions. It provides counts of each nucleotide in the specified regions.
#'
#' @param grObject A `GRanges` object representing genomic regions of interest. The sequence for these regions is extracted and analyzed.
#'
#' @return A `data.frame` with two columns:
#' \itemize{
#'   \item \code{Context}: The nucleotide context ("A", "T", "G", "C").
#'   \item \code{Counts}: The total counts of each nucleotide context within the given genomic regions.
#' }
#' 
#' @details
#' The function computes the nucleotide frequencies within the regions defined in \code{grObject}. Only nucleotides "T" and "C" are reported in the final output. The counts are based on the sum of sense strand nucleotides.
#'
#' @importFrom Biostrings letterFrequency getSeq
#' @importFrom GenomeInfoDb Hsapiens
#' @export
CompositionOnenucl0ntContext <- function(grObject) {
  Exon_nucleotides= data.frame("Sense_counts"=colSums(as.data.frame(letterFrequency(getSeq(Hsapiens,grObject), letters = c("A","T", "G","C")))))
  Composition1nt = data.frame(Context=c("A","T", "G","C"), Counts=c(rep(sum(Exon_nucleotides[c("A","T"),]), 2),rep(sum(Exon_nucleotides[c("G","C"),]), 2)))
  Composition1nt = Composition1nt[! (Composition1nt$Context %in% c("A","G")),]
  row.names(Composition1nt) = NULL
  return(Composition1nt)
}



#' Compute tri-Nucleotide Context Composition
#'
#' This function calculates the composition of trinucleotide contexts within the given genomic regions. It provides counts of each trinucleotide, considering both sense and antisense strands.
#'
#' @param grObject A `GRanges` object representing genomic regions of interest. The sequence for these regions is extracted and analyzed.
#'
#' @return A `data.frame` with two columns:
#' \itemize{
#'   \item \code{Context}: The trinucleotide context.
#'   \item \code{Counts}: The total counts of each trinucleotide context within the given genomic regions.
#' }
#'
#' @details
#' The function computes the trinucleotide frequencies within the regions defined in \code{grObject}. Only trinucleotides where the middle nucleotide is not "A" or "G" are included in the final output. Counts are aggregated from both sense and antisense strands.
#'
#' @importFrom Biostrings trinucleotideFrequency reverseComplement getSeq
#' @importFrom GenomeInfoDb Hsapiens
#' @export
CompositionTrinucl1ntContext <- function(grObject) {
  start(grObject) = start(grObject)-1
  end(grObject) = end(grObject)+1
  Exon_trinucleotides= data.frame("Sense_counts"=colSums(as.data.frame(trinucleotideFrequency(getSeq(Hsapiens,grObject)))))
  Exon_trinucleotides$Antisence_counts = colSums(as.data.frame(trinucleotideFrequency(reverseComplement(getSeq(Hsapiens,grObject)))))
  Exon_trinucleotides$Counts= as.numeric(Exon_trinucleotides$Sense_counts) +
    as.numeric(Exon_trinucleotides$Antisence_counts)
  Exon_trinucleotides$Context = rownames(Exon_trinucleotides)
  Exon_trinucleotides = Exon_trinucleotides[! (substr(Exon_trinucleotides$Context, 2,2) %in% c("A","G")),]
  row.names(Exon_trinucleotides) = NULL
  Composition3nt = data.frame(Context=Exon_trinucleotides$Context,Counts=Exon_trinucleotides$Counts)
  return(Composition3nt)
}



#' Compute penta-nucleotide Context Composition
#'
#' This function calculates the composition of pentanucleotide contexts within the given genomic regions. It provides counts of each pentanucleotide, considering both sense and antisense strands.
#'
#' @param grObject A `GRanges` object representing genomic regions of interest. The sequence for these regions is extracted and analyzed.
#'
#' @return A `data.frame` with two columns:
#' \itemize{
#'   \item \code{Context}: The pentanucleotide context.
#'   \item \code{Counts}: The total counts of each pentanucleotide context within the given genomic regions.
#' }
#'
#' @details
#' The function computes the pentanucleotide frequencies within the regions defined in \code{grObject}. Only pentanucleotides where the third nucleotide is not "A" or "G" are included in the final output. Counts are aggregated from both sense and antisense strands.
#'
#' @importFrom Biostrings oligonucleotideFrequency reverseComplement getSeq
#' @importFrom GenomeInfoDb Hsapiens
#' @export
CompositionPentanucl2ntContext <- function(grObject) {
  start(grObject) = start(grObject)-2
  end(grObject) = end(grObject)+2
  Exon_trinucleotides= data.frame("Sense_counts"=colSums(as.data.frame(oligonucleotideFrequency(getSeq(Hsapiens,grObject), width=5))))
  Exon_trinucleotides$Antisence_counts = colSums(as.data.frame(oligonucleotideFrequency(reverseComplement(getSeq(Hsapiens,grObject)), width=5)))
  Exon_trinucleotides$Counts= as.numeric(Exon_trinucleotides$Sense_counts) +
    as.numeric(Exon_trinucleotides$Antisence_counts)
  Exon_trinucleotides$Context = rownames(Exon_trinucleotides)
  Exon_trinucleotides = Exon_trinucleotides[! (substr(Exon_trinucleotides$Context, 3,3) %in% c("A","G")),]
  row.names(Exon_trinucleotides) = NULL
  Composition3nt = data.frame(Context=Exon_trinucleotides$Context,Counts=Exon_trinucleotides$Counts)
  return(Composition3nt)
}



#' Prepare Mutation Data for Regression Analysis
#'
#' This function prepares a dataset for regression analysis by calculating the number of mutations in various contexts (trinucleotide or pentanucleotide) 
#' and stratifying them based on sample cohorts and mutation patterns.
#'
#' @param grObject_within A `GRanges` object representing genomic regions of interest for context extraction.
#' @param grMutations A `GRanges` object containing mutation data.
#' @param HGNC_symbol An optional character string representing the HGNC symbol for a gene. If provided, it is used in the analysis.
#' @param annotation A `data.table` containing sample annotations with columns including at least "Sample" and cohort information.
#' @param clusterNumber An integer specifying the type of nucleotide context to use:
#' \itemize{
#'   \item \code{0}: Pentanucleotide context + matching.
#'   \item \code{1}: Trinucleotide context + matching.
#'   \item \code{96}: Trinucleotide context + stratification.
#'   \item \code{1536}: Pentanucleotide context + stratification.
#' }
#'
#' @return A `data.table` with the following columns:
#' \itemize{
#'   \item \code{Mutation}: Mutation type or pentaMatch or triMatch in case of matching.
#'   \item \code{Context}: The nucleotide context associated with the mutation.
#'   \item \code{Counts}: Total counts of each nucleotide context.
#'   \item \code{noSamples}: Number of samples in each cohort.
#'   \item \code{MutationNumber}: Number of mutations for each group.
#' }
#'
#' @details
#' The function determines possible mutation patterns based on the `clusterNumber` and calculates the mutation counts within specified genomic regions. It handles both trinucleotide and pentanucleotide contexts, including special cases for different clustering schemes. The result is a dataset ready for further regression analysis, including normalization for sample size and mutation counts.
#'
#' @importFrom Biostrings getSeq reverseComplement DNAStringSet
#' @importFrom GenomeInfoDb Hsapiens
#' @importFrom data.table data.table .SD .SDcols :=
#' @importFrom plyr count
#' @importFrom dplyr slice
#' @export
Mutations_toRegress <- function(grObject_within,
                                grMutations,
                                HGNC_symbol = NULL,
                                annotation = NULL,
                                clusterNumber = 0) {

      # ================================================
      # 1. Counting the number of samples in each cohort
      # based on the different classes of samples 
      # (e.g. treatment, cancer type)
      # ================================================
      sample_cohorts <- setdiff(colnames(annotation), "Sample")
      strat_features <- c("Mutation", sample_cohorts)
      # 
      SampleNumber_by_cancer_type <- plyr::count(annotation[, -"Sample", with = F])
      colnames(SampleNumber_by_cancer_type)[colnames(SampleNumber_by_cancer_type) == "freq"] <- "SampleNumber_by_cancer_type"
      #
      # Possible mutation patterns for trinucleotide or pentanucleotide matching 
      nucleotide_alphabet <- c("A", "C", "T", "G")
      pentanucleotide_combination <- lapply(seq_len(5), function(i) nucleotide_alphabet) # list(nucleotide_alphabet, nucleotide_alphabet, nucleotide_alphabet, nucleotide_alphabet, nucleotide_alphabet)
      #
      if (clusterNumber == 96) {
        # trinucleotide representation
        PossibleMutations <- apply(do.call(expand.grid, pentanucleotide_combination[1:3]), 1, paste, collapse = "")
        PossibleMutations <- apply(expand.grid(PossibleMutations, nucleotide_alphabet), 1, paste, collapse = ">")
        PossibleMutations <- data.frame("Mutation" = PossibleMutations[substr(PossibleMutations, 2, 2) != substr(PossibleMutations, 5, 5) & !(substr(PossibleMutations, 2, 2) %in% c("A", "G"))])
        #
        Composition_nucleotides <- CompositionTrinucl1ntContext(grObject_within)
        ##
      } else if (clusterNumber == 1536) {
        # pentanucleotide representation
        PossibleMutations <- apply(do.call(expand.grid, pentanucleotide_combination[1:5]), 1, paste, collapse = "")
        PossibleMutations <- apply(expand.grid(PossibleMutations, nucleotide_alphabet), 1, paste, collapse = ">")
        PossibleMutations <- data.frame("Mutation" = PossibleMutations[substr(PossibleMutations, 3, 3) != substr(PossibleMutations, 7, 7) & !(substr(PossibleMutations, 3, 3) %in% c("A", "G"))])
        #
        Composition_nucleotides <- CompositionPentanucl2ntContext(grObject_within)
        ##    
      } else if (clusterNumber == 1){
        # trinucleotide matching
        PossibleMutations <- data.frame("Mutation" = "triMatch")
        Composition_nucleotides <- data.frame(Context = "tri",Counts= sum(width(grObject_within)))
        ##
      } else if (clusterNumber == 0){
        # pentanucleotide representation
        PossibleMutations <- data.frame("Mutation" = "pentaMatch")
        Composition_nucleotides <- data.frame(Context = "penta",Counts= sum(width(grObject_within)))
        ##
      } 
      #
      # possible combinations of cohorts and mutation patterns
      annotation_file_cohort_combinations <- annotation[, -"Sample"] %>% unique()
      #
      noPossibleMutations <- dim(PossibleMutations)[1]
      noPossibleAnnotations <- dim(annotation_file_cohort_combinations)[1]
      annotation_file_cohort_combinations_dt <- annotation_file_cohort_combinations %>% dplyr::slice(rep(1:n(), each = noPossibleMutations))
      PossibleMutations_dt <- PossibleMutations %>% dplyr::slice(rep(1:n(), noPossibleAnnotations))
      PossibleMutations_by_cohort <- cbind(PossibleMutations_dt, annotation_file_cohort_combinations_dt)
      #
      # ===========================================
      # 3. Collecting mutations within the grObject
      # and creating a table with stratifications
      # ===========================================
      #
      Mutations <- findOverlaps(query = grMutations, subject = grObject_within, type = "within", ignore.strand = TRUE)
      if (length(Mutations) > 0) {
            Mutations.df <- data.table::data.table(as.data.frame(grMutations[queryHits(Mutations), ]))
            if (clusterNumber == 96) {
                  Mutations.df$Context <- c(as.character(getSeq(
                    Hsapiens, Mutations.df$seqnames, Mutations.df$start - 1,
                    Mutations.df$end + 1
                  )))
                  Mutations.df$Context <- as.character(paste(substr(Mutations.df$Context, 1, 1),
                                                            Mutations.df$ref_allele,
                                                            substr(Mutations.df$Context, 3, 3),
                                                            sep = ""
                  ))
                  Mutations.df$Mutation <- as.character(paste(Mutations.df$Context,
                                                              Mutations.df$mutated_to_allele,
                                                              sep = ">"
                  ))
                  ##
                  # reverseComplement of A or G mutations
                  cond_AG <- Mutations.df$ref_allele %in% c("A", "G")
                  if (nrow(Mutations.df[cond_AG, ]) > 0) {
                        Mutations.df[cond_AG, ]$Context <- as.character(reverseComplement(DNAStringSet(Mutations.df[cond_AG, ]$Context)))          
                        Mutations.df[cond_AG, ]$Mutation <- as.character(paste(Mutations.df[cond_AG, ]$Context,
                                                                      chartr("ATGC", "TACG", Mutations.df[cond_AG, ]$mutated_to_allele),
                                                                      sep = ">"))
                                                          }
            ###
            } else if(clusterNumber == 1536){
                  Mutations.df$Context <- c(as.character(getSeq(
                    Hsapiens, Mutations.df$seqnames, Mutations.df$start - 2,
                    Mutations.df$end + 2
                  )))
                  Mutations.df$Context <- as.character(paste(substr(Mutations.df$Context, 1, 2),
                                                            Mutations.df$ref_allele,
                                                            substr(Mutations.df$Context, 4, 5),
                                                            sep = ""
                  ))
                  Mutations.df$Mutation <- as.character(paste(Mutations.df$Context,
                                                              Mutations.df$mutated_to_allele,
                                                              sep = ">"
                  ))
                  # reverseComplement of A or G mutations
                  cond_AG <- Mutations.df$ref_allele %in% c("A", "G")
                  if (nrow(Mutations.df[cond_AG, ]) > 0) {
                        Mutations.df[cond_AG, ]$Context <- as.character(reverseComplement(DNAStringSet(Mutations.df[cond_AG, ]$Context)))          
                        Mutations.df[cond_AG, ]$Mutation <- as.character(paste(Mutations.df[cond_AG, ]$Context,
                                                                      chartr("ATGC", "TACG", Mutations.df[cond_AG, ]$mutated_to_allele),
                                                                      sep = ">"))
                                                          }              
            ###
            } else if(clusterNumber == 1){
                  Mutations.df$Context <- "tri"
                  Mutations.df$Mutation <- "triMatch"
            ###
            } else if(clusterNumber == 0){
                  Mutations.df$Context <- "penta"
                  Mutations.df$Mutation <- "pentaMatch"
            }
            ###
            # reduce number of columns with unneccessary now data:
            Mutations.df = Mutations.df[, -c("mutated_to_allele", "ref_allele", "seqnames", "start", "end", "width", "strand")]
            ##
            ##
            # Merging mutations with annotations
            Mutations.df = merge(Mutations.df, annotation, by = "Sample")
            ##
            ##
            # Count mutation of same type + cohort
            MutationNumber_percohort_andMutType <- plyr::count(Mutations.df[, ..strat_features]) %>% data.table()
            colnames(MutationNumber_percohort_andMutType)[colnames(MutationNumber_percohort_andMutType) == "freq"] <- "MutationNumber"
            #
            Mutations.df = merge(Mutations.df, MutationNumber_percohort_andMutType, by = strat_features)
            Mutations.df = merge(Mutations.df, PossibleMutations_by_cohort, by = strat_features, all.y = T)
            Mutations.df[is.na(Mutations.df$MutationNumber), MutationNumber := 0]
      } else {
            ## Mutations.df with no mutations
            Mutations.df = PossibleMutations_by_cohort
            Mutations.df = data.table::data.table(Mutations.df)
            # Mutations.df = merge(Mutations.df, SampleNumber_by_cancer_type, by = setdiff(strat_features, "Mutation"))
            Mutations.df$MutationNumber <- 0
      }
      #
      # Computing the nucleotide at risk
      Mutations.df = merge(Mutations.df, SampleNumber_by_cancer_type, by = setdiff(strat_features, "Mutation"))
      # Fill Context of combinations with no mutations
      if(clusterNumber == 1){
          Mutations.df[Mutations.df$MutationNumber == 0, Context := "tri"]  
      } else if(clusterNumber == 0){
          Mutations.df[Mutations.df$MutationNumber == 0, Context := "penta"]      
      } else{
          Mutations.df[Mutations.df$MutationNumber == 0, Context := substr(as.character(Mutation), 1, nchar(as.character(Mutation)) - 2)]  
      }  
      #
      Mutations.df = merge(Mutations.df, Composition_nucleotides[, c("Context", "Counts")], by = "Context")
      #
      ### For handling numeric features
      nonNumericFeatures <- dplyr::intersect(names(which(sapply(Mutations.df, class) != "numeric")), strat_features)
      Mutations.df[, (nonNumericFeatures) := lapply(.SD, as.factor), .SDcols = nonNumericFeatures]
      #
      Mutations.df = Mutations.df[, c(strat_features, "Context","Counts","SampleNumber_by_cancer_type", "MutationNumber"), with = F]
      colnames(Mutations.df)[colnames(Mutations.df) == "SampleNumber_by_cancer_type"] <- "noSamples"
      Mutations.df = unique(Mutations.df) # in case of mutational same rows can be left
      #
      return(Mutations.df)
}