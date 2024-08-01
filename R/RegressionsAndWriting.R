#' Summarize Regression Table
#'
#' This function processes a regression table by summarizing mutation data based on the specified mutation number (MSnumber).
#' It groups and aggregates data based on the matching method we are applying.
#'
#' @param TableToRegress A data frame or data table containing regression data, including mutation-related columns.
#' @param MSnumber An integer indicating the the DiffInvex matching approach.
#' MSnumber = 0: penta-nucleotide matching between target and background regions.
#' MSnumber = 1: tri-nucleotide matching between target and background regions.
#'
#' @return A processed data table with aggregated mutation data. 
#'
#' @import dplyr
#' @import data.table
#' @export
summarize_regress_table <- function(TableToRegress, MSnumber) {
  if (MSnumber==0 | MSnumber==1) { 
      #
      if(MSnumber== 1) {
        TableToRegress[, Mutation := "TriRelMatching"]
      } else {
        TableToRegress[, Mutation := "PentaRelMatching"]
      }
      #
      # Aggregate over Mutation types
      TableToRegress = TableToRegress %>%
        dplyr::group_by_at(vars(-ntAtRisk, -MutationNumber, -Mutation)) %>%
        dplyr::mutate(ntAtRisk = sum(ntAtRisk), MutationNumber = sum(MutationNumber),) %>% unique() %>% data.table()
      #
      TableToRegress$ln_ntAtRisk = as.numeric(log(TableToRegress$ntAtRisk))
      TableToRegress$Mutation = factor(TableToRegress$Mutation)
      return(TableToRegress)
      ##
  } else {
      #
      TableToRegress$ln_ntAtRisk = as.numeric(log(TableToRegress$ntAtRisk))
      return(TableToRegress)
  }
}



#' Perform Regression on a Data Table
#'
#' This function applies generalized linear model (GLM) regressions on a specified data table. 
#' The table is first preprocessed by creating unique identifiers for groups of interest, 
#' then summarized if necessary, and finally, regressions are run on either grouped or 
#' the entire dataset based on the number of grouping variables.
#'
#' @param TableToRegress A data frame or data table containing the data to be regressed.
#' @param clusterNumber An integer indicating the mutation type clustering; 
#' 0: penta-nucelotide matching and 1: tri-nucelotide matching.
#' @param regression_types A character vector specifying the types of regression models to be applied 
#' (e.g., "bayes.possion", "poisson").
#' @param colNamesToRegress A character vector specifying the names of columns to be included in the regression.
#'
#' @return A list of regression results, with each element corresponding to a different regression type. 
#' Each result is a data table of regression coefficients and statistics.
#'
#' @import dplyr
#' @import data.table
#' @export
regress_table = function(TableToRegress,
                         clusterNumber,
                         regression_types,
                         colNamesToRegress) {
  ##
  ##
  strat_features_forMSdescreasing = colnames(TableToRegress) %>% setdiff(c("isTarget", "ntAtRisk", "MutationNumber", "Mutation"))
  TableToRegress[, ID := Reduce(function(...) paste(..., sep = "__"), .SD[, mget(strat_features_forMSdescreasing)])]
  ##
  ##
  # reduce mutation table if clusterNumber %in% c(0,1)
  TableToRegress_clustered = rbindlist(lapply(unique(TableToRegress$ID), function(id) {
    return(summarize_regress_table(TableToRegress[ID==id, -"ID"],clusterNumber))
  }))
  TableToRegress[, ID := NULL]
  TableToRegress_clustered = data.table::data.table(TableToRegress_clustered)
  colnames_perGroup =  strat_features_forMSdescreasing %>% setdiff(colNamesToRegress)
  ##
  ##
  ##
  # Appling GLM regression 
  # 1) No more variables, apply regression once
  if (length(colnames_perGroup) == 0) {
    Gene_coeffs = lapply(regression_types, function(x) {
                  glm_output = glm.safe(data = TableToRegress_clustered,
                            family = x,
                            targetVar = "MutationNumber",
                            offsetVar = "ntAtRisk",
                            featureVars = colNamesToRegress,
                            interVars =  NULL,
                            providePval = TRUE,
                            contrasts = NULL,
                            provideConfint = FALSE,
                            baseIterationNumber = 25)
                  glm_output = glm_output[!(str_detect(glm_output$coefName, "Mutation")),]
                  return(glm_output)
                  })
  } else {
  ##
  ##
  # 2) if there are more variables than ones that you are controlling for, than the regression will run saparately for each level of variable that is not in this file.  
    setDT(TableToRegress_clustered)[, ID := Reduce(function(...) paste(..., sep = "__"), .SD[, mget(colnames_perGroup)])]
    unique_IDS = unique(TableToRegress_clustered$ID)
    #    
    N = length(unique_IDS)
    Gene_coeffs <- lapply(1:N, function(i) {
                  #cat(paste0("Working on group of ", unique_IDS[i], " (", i, "/", N, ")\n"))
                  TableToRegress_ID = TableToRegress_clustered[ID == unique_IDS[i], ]
                  outputs = lapply(regression_types, function(x) {
                              glm_output = glm.safe(data = TableToRegress_ID,
                                                  family = x,
                                                  targetVar = "MutationNumber",
                                                  offsetVar = "ntAtRisk",
                                                  featureVars = colNamesToRegress,
                                                  interVars =  NULL,
                                                  providePval = TRUE,
                                                  contrasts = NULL,
                                                  provideConfint = FALSE,
                                                  baseIterationNumber = 25)
                              glm_output = glm_output[!(str_detect(glm_output$coefName, "Mutation")),]  
                              glm_output =  glm_output[,ID := unique_IDS[i]]
                              return(glm_output)})
                  names(outputs) = regression_types
                  return(outputs)
                  })
    #
    Gene_coeffs = lapply(regression_types, function(regression_type) {
                  rbindlist(lapply(Gene_coeffs, function(x) x[[regression_type]]))
                  })
  }
  #
  names(Gene_coeffs) = regression_types
  return(Gene_coeffs)
}



#' Perform Regression on the input table
#'
#' This function reads a file containing controlled variables and prepares a data table for regression analysis.
#' It applies GLM regression across specified variables and returns summarized regression results.
#'
#' @param TableToRegress A data frame or data table containing data for regression.
#' @param clusterNumber An integer specifying the type of matching (penta-nucelotide:0, tri-nucleotide:1).
#' @param regression_types A character vector of regression types to be applied.
#' @param controlled_variables_input A file path to a list of variables to control for in the regression.
#' @param HGNC_symbol A character string representing the gene symbol for which regression is being computed.
#'
#' @return A list of summarized regression results, each element corresponding to a different regression type.
#'
#' @import dplyr
#' @import data.table
#' @export
compute_regress_table <- function(TableToRegress,
                                 clusterNumber,
                                 regression_types,
                                 controlled_variables_input,
                                 HGNC_symbol){
  # Reading variables to control for
  colNamesToRegress = fread(controlled_variables_input, header = F)
  colNamesToRegress = colNamesToRegress$V1  
  #
  # Excluding mutation class ("Mutation") from variables to control for in we apply trinucleotide/pentanucleotide matching 
  if (clusterNumber==0 | clusterNumber==1) {
    colNamesToRegress = colNamesToRegress %>% setdiff("Mutation")
  }
  #
  # cat("variable to control for: \n")
  # print(colNamesToRegress)
  #
  for (columnNameToRegress in colNamesToRegress) {
    if (!(str_detect(columnNameToRegress, ":"))) {
      if  ("0" %in% levels(TableToRegress[[columnNameToRegress]]) ) {
        TableToRegress[[columnNameToRegress]] = relevel(TableToRegress[[columnNameToRegress]], "0")
    }}
  }
  #
  # Applying the GLM regression iteratively
  ntimes <- 1
  regression_output = lapply(1:ntimes, function(seed){
                      set.seed(seed)
                      tryCatch({regress_table(TableToRegress, clusterNumber, regression_types, colNamesToRegress)},
                      error=function(err){  NULL  },  
                      finally = {})
                      })
  # 
  regression_output = lapply(regression_types, function(regression_type) {
                     rbindlist(lapply(regression_output, function(x) x[[regression_type]]))
                     })
  names(regression_output) = regression_types
  #
  # summarize the regression table
  regression_output_summarized = lapply(regression_types, function(regression) {
                                 single_regression_output = regression_output[[regression]]
                                 if (!(is.null(single_regression_output))) {
                                  # for  having one summarized result instead of many
                                    if ("coefName" %in% colnames(single_regression_output)) {
                                        if (! "ID" %in% colnames(single_regression_output)) {
                                              single_regression_output = single_regression_output %>% 
                                                                                      data.frame() %>%
                                                                                      dplyr::group_by(coefName, converged) %>% 
                                                                                      dplyr::mutate(coef=median(as.numeric(coef), na.rm=T),
                                                                                            pVal = exp(median(log(pVal), na.rm=T)),
                                                                                            upperTailPval = exp(median(log(as.numeric(upperTailPval)), na.rm=T)),
                                                                                            lowerTailPval = exp(median(log(as.numeric(lowerTailPval)), na.rm=T)),
                                                                                            Std.error = median(as.numeric(Std.error), na.rm=T)) %>%
                                                                                      unique() %>%
                                                                                      data.table()
                                        } else {
                                              single_regression_output = single_regression_output %>% 
                                                                                      data.frame() %>%
                                                                                      dplyr::group_by(coefName, converged, ID) %>% 
                                                                                      dplyr::mutate(coef=median(as.numeric(coef), na.rm=T),
                                                                                            pVal = exp(median(log(pVal), na.rm=T)),
                                                                                            upperTailPval = exp(median(log(as.numeric(upperTailPval)), na.rm=T)),
                                                                                            lowerTailPval = exp(median(log(as.numeric(lowerTailPval)), na.rm=T)),
                                                                                            Std.error = median(as.numeric(Std.error), na.rm=T)) %>%
                                                                                      unique() %>%
                                                                                      data.table()
                                        }
                                    }
                                }
                                #
                                single_regression_output$HGNC = HGNC_symbol
                                return(single_regression_output)
                                })
  #
  names(regression_output_summarized) = regression_types
  return(regression_output_summarized)
}



#' Write Summarized Regression Table to File
#'
#' This function writes a summarized version of the regression table to a specified directory.
#'
#' @param TableToRegress A data frame or data table containing regression data to be summarized and saved.
#' @param directory A string specifying the directory where the output file will be saved.
#' @param HGNC_symbol A character string representing the gene symbol to be included in the output.
#' @param Method A string indicating the method used for regression, included in the output filename.
#'
#' @return None. The function writes the output directly to a file.
#'
#' @import dplyr
#' @import data.table
#' @export
writing_table <- function(TableToRegress, directory, HGNC_symbol, Method){
    #
    if (!dir.exists(directory)){
        dir.create(directory)
    }
    #
    suppressMessages(TableToRegress <- TableToRegress %>% group_by_at(setdiff(colnames(TableToRegress), c("MutationNumber","ntAtRisk"))) %>%
                                         dplyr::summarize(MutationNumber = sum(MutationNumber), ntAtRisk = sum(ntAtRisk)))    
    #  
    TableToRegress$HGNC = HGNC_symbol
    filename = paste(directory, "/",HGNC_symbol, "_TableToRegress_", Method, ".txt", sep = "")
    tryCatch({fwrite(TableToRegress, filename, sep=",", row.names=F, quote=T) },
          error=function(err){  NULL  },  
          finally = {}
          )
}



#' Write Regression Coefficients to Files
#'
#' This function writes regression coefficients for different regression types to separate CSV files 
#' in a specified directory. Each file includes the HGNC symbol and is named based on the regression 
#' type and method used.
#'
#' @param Gene_coeffs A list of data tables containing regression coefficients for each regression type.
#' @param directory A string specifying the directory where the output files will be saved.
#' @param HGNC_symbol A character string representing the gene symbol to be included in each output file.
#' @param Method A string indicating the method used for regression, included in the output filenames.
#'
#' @return None. The function writes the output directly to files.
#'
#' @import data.table
#' @export
writing_files <- function(Gene_coeffs, directory, HGNC_symbol, Method){
    #
    if (!dir.exists(directory)){
        dir.create(directory)
    }
    #
    regression_types = names(Gene_coeffs)
    temp = lapply(regression_types, function(regression) {
      if (!(is.null(Gene_coeffs[[regression]]))) {
          Gene_coeffs[[regression]]$HGNC = HGNC_symbol
          filename = paste(directory, "/",HGNC_symbol, "_", regression, "_", Method, ".csv", sep = "")
          write.csv(Gene_coeffs[[regression]],filename, row.names = F)
      }
    })
}