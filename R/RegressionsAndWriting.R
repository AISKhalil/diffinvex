summarize_regress_table = function(TableToRegress, MSnumber) {
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
###
###
###
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
  #print(TableToRegress_clustered)
  #print(colNamesToRegress)
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
###
###
###
compute_regress_table = function(TableToRegress,
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
  ##
  ##
  # Applying the GLM regression iteratively
  ntimes <- 1
  regression_output = lapply(1:ntimes, function(seed){
                      set.seed(seed)
                      tryCatch({regress_table(TableToRegress, clusterNumber, regression_types, colNamesToRegress)},
                      error=function(err){  NULL  },  
                      finally = {})
                      })
  # print(regression_output)
  # 
  regression_output = lapply(regression_types, function(regression_type) {
                     rbindlist(lapply(regression_output, function(x) x[[regression_type]]))
                     })
  names(regression_output) = regression_types
  ##
  ##
  # summarize the regression table
  regression_output_summarized = lapply(regression_types, function(regression) {
                                 single_regression_output = regression_output[[regression]]
                                 if (!(is.null(single_regression_output))) {
                                  # for  having one summarized result instead of many
                                    if ("coefName" %in% colnames(single_regression_output)) {
                                        if ("ID" %!in% colnames(single_regression_output)) {
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
##
##
writing_table = function(TableToRegress, directory, HGNC_symbol, Method){
    #
    if (!dir.exists(directory)){
        dir.create(directory)
    }
    #
    #suppressMessages(TableToRegress <- TableToRegress %>% group_by_at(setdiff(colnames(TableToRegress), c("Mutation","MutationNumber","ntAtRisk"))) %>%
    #                                     dplyr::summarize(MutationNumber = sum(MutationNumber), ntAtRisk = sum(ntAtRisk)))
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
##
##
writing_files = function(Gene_coeffs, directory, HGNC_symbol, Method){
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
##
##
