###  loading data abd libraries ##
suppressMessages(library("ddpcr"))
suppressMessages(library("optparse"))
#
packages <- c("arm","biomaRt","rtracklayer","GenomicRanges","GenomicFeatures", "dplyr", "reshape2", "data.table","tidyverse",
               "readxl","qvalue","Matrix","ggplot2")
#
invisible(quiet(lapply(packages, library, character.only = TRUE), all = T)) # invisible for suppressing output of the lapply, suppresMessages for libraries
###
###
###
diffinvex_drivers = function(output_directory,
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
                fwrite(annResultsAll, file = paste0(figures_directory,"/cohort_MutationProfile_diffInVEx.tsv"), sep="\t", row.names=F, quote=T)
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
            x <- fread(file = paste0(figures_directory,"/cohort_MutationProfile_diffInVEx.tsv"))
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
    #
    ## single or multiple cohorts
    fileName     <- paste0(output_directory,"/",genesTested[1],regression_type)
    geneResults  <- read.csv(fileName, header = TRUE)
    singleCohort <- !(c("ID") %in% colnames(geneResults))
    ##
    ##
    ##
    if(singleCohort){
                  ## results ##
                  regResults <- lapply(genesTested, function(x){
                                    fileName <- paste0(output_directory,"/",x,regression_type)
                                    geneResults <- read.csv(fileName)
                                    geneResults <- geneResults[geneResults$coefName == regression_var & geneResults$converged == TRUE,]
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
                  write.csv(tumorRegResults[c("gene","coef","pVal","Std.error","pVal_adj","qVal")], paste0(figures_directory,"/cohort_",regression_var,"_",reg_type,"_all.csv"), row.names = FALSE)
                  #
                  #
                  sigGenes <- rownames(tumorRegResults[tumorRegResults$qVal < fdr,])
                  write.table(sigGenes, file = paste0(figures_directory,"/cohort_",regression_var,"_",reg_type,"_sigGenes_qval_",fdr,"percent.csv"), col.names = FALSE, row.names = FALSE, sep = "")
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
                  figName2 <- paste0(figures_directory,"/cohort_",regression_var,"_",reg_type,"_DiffInVEx_pValues_QQplot_lambda.png")
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
                                    geneResults <- geneResults[geneResults$coefName == regression_var & geneResults$converged == TRUE,]
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
                    write.csv(tumorRegResults[[i]][c("gene","coef","pVal","Std.error","pVal_adj","qVal")], paste0(figures_directory,"/",i,"_",regression_var,"_",reg_type,"_all.csv"), row.names = FALSE)
                    #
                  }
                  ##
                  ## QQ-plot (p-value) ##
                  tumorRegResults <- list()
                  #
                  for(i in IDs){
                    tumorRegResults[[i]] <- read.csv(paste0(figures_directory,"/",i,"_",regression_var,"_",reg_type,"_all.csv"))
                    #
                    pValues <- tumorRegResults[[i]]$pVal
                    genes   <- tumorRegResults[[i]]$gene 
                    #
                    genes <- genes[order(pValues)]
                    pValues <- sort(pValues)
                    ##
                    ##
                    figName1 <- paste0(figures_directory,"/",i,"_pValues_qqplot_DiffInVEx_",regression_var,"_",reg_type,".png")
                    testPvalues(pValues, figName1)
                    ##
                    ##
                    figName2 <- paste0(figures_directory,"/",i,"_pValues_QQplot_lambda_DiffInVEx_",regression_var,"_",reg_type,".png")
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
############################################
############################################
###-------- plotting p-values -----------###
#
testPvalues = function(pvalues, qqPlotFile){
  #
  #png(qqPlotFile, width = 800, height = 600)
  #qqunif(pvalues, logscale = TRUE, col = col)
  #dev.off()
  #
  n =length(pvalues)
  y = runif(n,min = 0, max = 1)
  png(qqPlotFile, width = 600, height = 600, res = 100)
  attach(USJudgeRatings)
  q <- qqplot(-log10(y),-log10(pvalues), main="qq-plot (uniform distribution)", 
         xlab = "-log10(Expected p-value)",
         ylim = c(0,15),
         xlim = c(0,4.5),
         ylab = "-log10(Observed p-value)",
         lwd = 0.5) 
  q + theme(axis.text=element_text(size=20),
            axis.title=element_text(size=14,face="bold")) 
  abline(coef = c(0,1), lwd=2)
  detach(USJudgeRatings)
  dev.off()
  #
  return(ks.test(pvalues,"punif",0,1))
}
#
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
#
inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}
