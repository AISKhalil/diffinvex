suppressMessages(library("ddpcr"))
suppressMessages(library("optparse"))
#
packages <- c("arm","biomaRt","rtracklayer","GenomicRanges","GenomicFeatures", "dplyr", "reshape2", "data.table","tidyverse",
               "readxl","qvalue","Matrix","ggplot2")
#
invisible(quiet(lapply(packages, library, character.only = TRUE), all = T)) # invisible for suppressing output of the lapply, suppresMessages for libraries



#' Identify DiffInVEx Putative Drivers and Generate Results
#'
#' This function analyzes the results of a DiffInVEx (Differential Insertion/Deletion Variants and Exonic Variants) analysis 
#' to identify putative drivers. It generates various outputs, including mutation profiles, regression coefficients, 
#' and figures such as QQ-plots of p-values. The function supports single or multiple cohorts and applies multiple 
#' filtering criteria based on the provided parameters.
#'
#' @param output_directory Character string. Path to the directory containing the output files from the DiffInVEx analysis.
#' @param regression_var Character string. Name of the regression variable to be used for analysis (default is "isTarget1").
#' @param reg_type Character string. Type of regression model used (default is "bayes.poisson").
#' @param fdr Numeric. False discovery rate threshold for significance (default is 0.1).
#' @param genesSelected Character vector. List of selected genes for filtering. If `NULL`, no gene filtering is applied (default is `NULL`).
#' @param mutationsTh Integer. Threshold for the number of mutations to filter genes (default is 10).
#'
#' @section Process:
#' The function performs the following steps:
#' \itemize{
#'   \item Creates a directory for figures if it does not already exist.
#'   \item Outputs mutation profiles for each gene analyzed, saving the results to a TSV file.
#'   \item Outputs regression coefficients for the specified regression variable and type, applying FDR correction.
#'   \item Generates QQ-plots of p-values for the regression results and saves them as PNG files.
#' }
#'
#' @return None. The function generates and saves various results and figures to the specified `output_directory`.
#'
#' @export
diffinvex_drivers <- function(output_directory,
                                       regression_var = "isTarget1",
                                       reg_type = "bayes.poisson",
                                       fdr = 0.1,
                                       genesSelected = NULL,
                                       mutationsTh = 10){
    ##------------------ output figures ------------------##
    figures_directory <- paste0(output_directory,"_figures")
    if (!dir.exists(figures_directory)){
      dir.create(figures_directory)
    }
    ##
    ##----------- output mutation profile -----------##
    outFiles     <- list.files(output_directory)
    suffix       <- paste0("_MutationProfile_diffInVEx.tsv")
    regOutFiles  <- outFiles[grep(suffix, outFiles)]
    genesTested  <- sub(suffix,"",regOutFiles)
    if(length(genesTested) > 0){
                annResults <- lapply(genesTested, function(x){
                                  fileName <- paste0(output_directory,"/",x,suffix)
                                  geneResults <- fread(fileName, header = TRUE)
                                  geneResults$Sample <- as.character(geneResults$Sample)
                                  return(geneResults)
                              })
                #
                annResultsAll <- bind_rows(annResults)
                #
                noMutations <- dim(annResultsAll)[1]
                annResultsAll$DiffInVEx_ID <- seq(1,noMutations,1)
                #
                fwrite(annResultsAll, file = paste0(figures_directory,"/cohort_mutation_profile.tsv"), sep="\t", row.names=F, quote=T)
    }    
    ##
    ##----------- output coefficients -----------##
    outFiles        <- list.files(output_directory)
    regression_type <- paste0("_",reg_type,"_diffInVEx.csv")
    regOutFiles     <- outFiles[grep(regression_type, outFiles)]
    genesTested     <- sub(regression_type,"",regOutFiles)
    #
    if(length(genesSelected) != 0){
        if(genesSelected[1] == "genesWtargetMut"){
            x <- fread(file = paste0(figures_directory,"/cohort_mutation_profile.tsv"))
            x <- x %>% dplyr::filter(DiffInVEx_class == "target")
            x <- x %>% group_by_at(c("DiffInVEx_gene")) %>% dplyr::summarize(MutationNumber = length(Sample))
            #
            x <- x %>% dplyr::filter(MutationNumber > mutationsTh) %>% dplyr::select(DiffInVEx_gene)
            genesTested <- intersect(genesTested,x$DiffInVEx_gene)       
        }
        else {
            genesTested <- intersect(genesTested,genesSelected)
        }
    }
    ##
    ## single or multiple cohorts / regression variables ##
    fileName     <- paste0(output_directory,"/",genesTested[1],regression_type)
    geneResults  <- read.csv(fileName, header = TRUE)
    singleCohort <- !(c("ID") %in% colnames(geneResults))
    regression_variables <- grep(regression_var, geneResults$coefName, value = TRUE)
    ##
    ## output for each regression variables ##
    for(reg_var in regression_variables){
            ##
            ##
            if(singleCohort){
                          ## results ##
                          regResults <- lapply(genesTested, function(x){
                                            fileName <- paste0(output_directory,"/",x,regression_type)
                                            geneResults <- read.csv(fileName)
                                            geneResults <- geneResults[geneResults$coefName == reg_var & geneResults$converged == TRUE,]
                                            geneResults <- geneResults[c("coef","pVal","Std.error","HGNC")]
                                        })
                          ##
                          ## FDR ##
                          tumorRegResults    <- bind_rows(lapply(regResults, function(x){x[,c("coef","pVal","Std.error")]}))
                          convergedGeneNames <- unlist(lapply(regResults, function(x){x[,c("HGNC")]}))
                          #
                          rownames(tumorRegResults) <- convergedGeneNames
                          tumorRegResults$gene <- convergedGeneNames
                          #
                          #
                          tumorRegResults$pVal_adj <- p.adjust(tumorRegResults$pVal, method = "fdr", n = length(tumorRegResults$pVal)) 
                          tumorRegResults$qVal <- qvalue(tumorRegResults$pVal)$qvalues
                          tumorRegResults <- tumorRegResults[order(tumorRegResults$qVal),]
                          #
                          write.csv(tumorRegResults[c("gene","coef","pVal","Std.error","pVal_adj","qVal")], paste0(figures_directory,"/cohort_",reg_var,"_",reg_type,"_all.csv"), row.names = FALSE)
                          ##
                          ##
                          ## QQ-plot (p-value) ##
                          #
                          pValues <- tumorRegResults$pVal
                          genes   <- tumorRegResults$gene 
                          #
                          genes <- genes[order(pValues)]
                          pValues <- sort(pValues)
                          ##
                          ##
                          figName2 <- paste0(figures_directory,"/cohort_",reg_var,"_",reg_type,"_pValues_QQplot.png")
                          ##
                          print(paste0("λ = ",inflation(pValues)))
                          ##
                          gg_qqplot(pValues) +
                          theme_bw(base_size = 10) +
                          theme(axis.ticks = element_line(linewidth = 0.4), panel.grid = element_blank()) +
                          annotate(
                          geom = "text",
                          x = -Inf,
                          y = Inf,
                          hjust = -0.15,
                          vjust = 1 + 0.15 * 3,
                          label = sprintf("λ = %.2f", inflation(pValues)),
                          size = 5) +
                          coord_cartesian(xlim = c(0,4), expand=0) + 
                          coord_cartesian(ylim = c(0,15), expand=0)      
                          ggsave(figName2, width = 10, height = 10, units = "cm")
            } else {
            ## results ##
                          regResults <- lapply(genesTested, function(x){
                                            fileName <- paste0(output_directory,"/",x,regression_type)
                                            geneResults <- read.csv(fileName)
                                            geneResults <- geneResults[geneResults$coefName == reg_var & geneResults$converged == TRUE,]
                                            geneResults <- geneResults[c("coef","pVal","Std.error","HGNC","ID")]
                                        })
                          ## FDR ##
                          IDs <- (unique(bind_rows(lapply(regResults, function(x){x[c("ID")]}))))$ID
                          print(IDs)
                          tumorRegResults <- list()
                          #
                          for(i in IDs){
                            print(i)
                            #
                            y <- lapply(regResults, function(x){x[x$ID==i,c("coef","pVal","Std.error")]})
                            tumorRegResults[[i]] <- bind_rows(y)
                            convergedGeneNames <- unlist(lapply(regResults, function(x){x[x$ID==i,c("HGNC")]}))
                            #
                            rownames(tumorRegResults[[i]]) <- convergedGeneNames
                            tumorRegResults[[i]]$gene <- convergedGeneNames
                            #
                            tumorRegResults[[i]]$pVal_adj <- p.adjust(tumorRegResults[[i]]$pVal, method = "fdr", n = length(tumorRegResults[[i]]$pVal)) 
                            tumorRegResults[[i]]$qVal <- qvalue(tumorRegResults[[i]]$pVal)$qvalues
                            tumorRegResults[[i]] <- tumorRegResults[[i]][order(tumorRegResults[[i]]$qVal),]
        		                #
                            write.csv(tumorRegResults[[i]][c("gene","coef","pVal","Std.error","pVal_adj","qVal")], paste0(figures_directory,"/",i,"_",reg_var,"_",reg_type,"_all.csv"), row.names = FALSE)
                            #
                          }
                          ##
                          ## QQ-plot (p-value) ##
                          tumorRegResults <- list()
                          #
                          for(i in IDs){
                            tumorRegResults[[i]] <- read.csv(paste0(figures_directory,"/",i,"_",reg_var,"_",reg_type,"_all.csv"))
                            #
                            pValues <- tumorRegResults[[i]]$pVal
                            genes   <- tumorRegResults[[i]]$gene 
                            #
                            genes <- genes[order(pValues)]
                            pValues <- sort(pValues)
                            ##
                            ##
                            figName2 <- paste0(figures_directory,"/",i,"_pValues_QQplot_",reg_var,"_",reg_type,".png")
                            #
                            gg_qqplot(pValues) +
                            theme_bw(base_size = 10) +
                            theme(axis.ticks = element_line(linewidth = 0.4), panel.grid = element_blank()) +
                            annotate(
                            geom = "text",
                            x = -Inf,
                            y = Inf,
                            hjust = -0.15,
                            vjust = 1 + 0.15 * 3,
                            label = sprintf("λ = %.2f", inflation(pValues)),
                            size = 5) +
                            coord_cartesian(xlim = c(0,4), expand=0) + 
                            coord_cartesian(ylim = c(0,15), expand=0)      
                            ggsave(figName2, width = 10, height = 10, units = "cm")
                          }
              }
    }
}



#' Generate a QQ-Plot for p-values
#'
#' This function creates a QQ-plot to visualize the distribution of p-values. The QQ-plot compares the observed 
#' p-values against the expected distribution under the null hypothesis, with confidence intervals provided for 
#' the expected values.
#'
#' @param ps Numeric vector. The p-values to be plotted.
#' @param ci Numeric. The confidence interval level to be used for the confidence bands (default is 0.95).
#'
#' @return A ggplot2 object representing the QQ-plot of the p-values.
#'
#' @details
#' The function calculates the observed and expected -log10(p-values) and plots them against each other. 
#' It also includes confidence intervals for the expected -log10(p-values) using the beta distribution.
#'
#' The plot includes:
#' \itemize{
#'   \item A ribbon representing the confidence interval for the expected -log10(p-values).
#'   \item Points representing the observed -log10(p-values).
#'   \item A diagonal line indicating where the observed p-values would fall if they were uniformly distributed.
#' }
#'
#' @import ggplot2
#'
#' @examples
#' # Example usage:
#' p_values <- runif(1000) # Generate random p-values
#' gg_qqplot(p_values)
#'
#' @export
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
  geom_ribbon(
    mapping = aes(x = expected, ymin = clower, ymax = cupper),
    alpha = 0.1
  ) +
  geom_point(aes(expected, observed), shape = 1, size = 2) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
  xlab(log10Pe)  +
  ylab(log10Po)
}



#' Calculate the Inflation Factor for p-values
#'
#' This function computes the inflation factor (λ) for a set of p-values. The inflation factor is used to assess 
#' the extent of deviation from the null hypothesis in a dataset, often in the context of genomic studies. It 
#' helps to detect whether there is an excess of significant p-values compared to what would be expected under 
#' the null hypothesis.
#'
#' @param ps Numeric vector. The p-values for which the inflation factor is to be calculated.
#'
#' @return Numeric. The inflation factor (λ), which is the median of the chi-squared statistics divided by 
#' the chi-squared value at the 0.5 quantile.
#'
#' @export
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}
