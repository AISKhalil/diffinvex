#' NOTE: this should be a drop-in replacement for the glm.nb, with sensible defaults in case of failure
#' Assumes that the features are categorical variables, with levels "0", "1", "2" etc...
#' 
#' Runs a glm.nb and returns coefficients and confidence intervals. 
#' Does not throw an exception if the regression doesn't work, but instead provides sensible defaults. 
#' 
#' 1) if glm.nb fails, returns NA for all coefficients.
#' 2) if glm.nb runs ok but confint() fails, returns lower/upper CI bounds that are equal to the coefficient itself.
#'
#' One row per coefficient. THESE ARE IN BASE E (i.e. not yet converted to base 2)
#' 
#' TODO could make a wrapper function to parallelize this
#' 
#' @param providePvals  Is FALSE by default for compatibility reasons; if true, gives one additional column with two-tail P-values (from "summary" function that uses a Z-test)

`%+%` <- function(a, b) paste0(a, b)
library(dplyr)
library(data.table)
library(MASS)
library(stringr)
library(fastglm)
library(arm)
library(pscl)
.fastglm_fit <- function(x,
                         y,
                         w = rep(1, NROW(y)),
                         family = gaussian(),
                         offset = rep(0, NROW(y)),
                         start = NULL,
                         mustart = NULL,
                         etastart = NULL,
                         control = list(maxit = 100, epsilon = 1e-07),
                         intercept) {
                                    fastglmPure(
                                      x = x,
                                      y = y,
                                      weights = w,
                                      family = family,
                                      offset = offset,
                                      start = start,
                                      etastart = etastart,
                                      mustart = mustart,
                                      method = 2L,
                                      tol = control$epsilon,
                                      maxit = control$maxit
                                      )
                        }

prior.scale_parameter = 1 #5/1 by default
.bayesglm_fit <- function(x,
                          y,
                          w = rep(1, NROW(y)),
                          family = gaussian(),
                          offset = rep(0, NROW(y)),
                          start = NULL,
                          mustart = NULL,
                          etastart = NULL,
                          control = list(maxit = 100, epsilon = 1e-07),
                          intercept) {
                                      bayesglm.fit(
                                        x=x,
                                        y=y,
                                        weights = w,
                                        start = start,
                                        etastart = etastart,
                                        mustart = mustart,
                                        offset = offset,
                                        family = family,
                                        control = control,
                                        prior.scale = prior.scale_parameter, #5 by default
                                        prior.mean.for.intercept = -14
                                      )
                          }
###
###
###
glm.safe = function(data,
                    family,
                    targetVar,
                    offsetVar,
                    featureVars,
                    interVars = NULL,
                    providePval = FALSE,                                           
                    contrasts,
                    provideConfint = FALSE,
                    baseIterationNumber = 100 ) {
  ##
  ##
  if  (family %!in% c("loglin.poisson", "loglin.bayes.poisson", "Fisher.OR", "Log.OR")) {
            # for linear regression with offset
            aData = data
            aData[[offsetVar]] = log(aData[[offsetVar]])
            
            aFormula = targetVar %>% 
              paste0(" ~ offset(") %>% 
              paste0(offsetVar) %>% paste0(") + ") %>% 
              paste0(paste(featureVars, collapse=" + " ))
            
            if (!is.null(interVars)) {
                maxN = length(interVars)
                interVars = lapply(2:maxN, function(x) combn(interVars, x)) 
                interVars = lapply(interVars, function(x) {
                    y=split(x, rep(1:ncol(x), each = nrow(x)))
                    y = lapply(y, function(x) paste(x, collapse="` : `")) %>% unlist()
                    }) %>% unlist()
                aFormula = aFormula %>% 
                            paste0(" + `", paste(interVars, collapse=" + `"), "`" )
              }
  } else {
            # for LLA and Fisher OR exact test
            bData_mutated = data
            bData_mutated$isMutated = "1"
            bData_not_mutated = data
            bData_not_mutated$isMutated = "0"
            bData_not_mutated[[targetVar]] = bData_not_mutated[[offsetVar]] - bData_not_mutated[[targetVar]]
            bData = rbind(bData_mutated, bData_not_mutated)
            rm(bData_mutated, bData_not_mutated)
            bData[[offsetVar]] = NULL
            
            bFormula = targetVar %>% 
              paste0(" ~ ") %>% 
              paste0("", paste(featureVars, collapse=" + " ))
            interVars = featureVars[str_detect(featureVars, "isTarget")]
            
            if (!is.null(interVars)) {
              interVars = lapply(interVars, function(intervar) {
                singleintervar = str_split(intervar, ":") %>% unlist() %>% unique()
                singleintervar = paste(singleintervar, "isMutated", sep=":")
                intervar = c(singleintervar, paste(intervar, "isMutated", sep=":"))
              }) %>% unlist() %>% unique()
              
              bFormula = bFormula %>% 
                paste0(" + ", paste(interVars, collapse=" + "))
              }  
  }
  ##
  ##
  Glm_model =	tryCatch({
    
    if (family == "negative.binomial") {
          invisible(glm.nb(as.formula(aFormula), 
                           data = aData, 
                           control = glm.control(maxit = baseIterationNumber), 
                           contrasts = contrasts, 
                           method = ".fastglm_fit"))
    } else if (family == "Fisher.OR") {
          bMatrix=as.matrix(dcast(bData, isTarget~isMutated, value.var="MutationNumber", fun.aggregate=mean, na.rm = TRUE),  rownames = "isTarget")
          pseudocounts = as.vector(rowSums(bMatrix))%*%t(as.vector(colSums(bMatrix)))*2/(sum(as.vector(rowSums(bMatrix)))^2)
          pseudocounts = 1
          bMatrix = bMatrix+pseudocounts
          fisher.test(bMatrix)
          #chisq.test(bMatrix, simulate.p.value = TRUE)
    } else if (family == "Log.OR") {
          bMatrix=as.matrix(dcast(bData, isTarget~isMutated, value.var="MutationNumber", fun.aggregate=mean, na.rm = TRUE),  rownames = "isTarget")
          pseudocounts = as.vector(rowSums(bMatrix))%*%t(as.vector(colSums(bMatrix)))*2/(sum(as.vector(rowSums(bMatrix)))^2)
          #  pseudocounts = 1
          bMatrix = bMatrix+pseudocounts
          bMatrix = as.data.table(bMatrix)
          logoddsRatio = bMatrix$"1"/ bMatrix$"0"
          logoddsRatio = logoddsRatio[2]/logoddsRatio[1]
          #
          pvalue = chisq.test(bMatrix, simulate.p.value = TRUE)$p.value
          list(pvalue=pvalue, OR = log(logoddsRatio))
    } else if (family == "poisson") {
          mm <- model.matrix(as.formula(aFormula), data = aData, contrasts.arg = contrasts)
          invisible(fastglm(x = mm,
                            y = aData[[targetVar]],
                            family = poisson(link = "log"),
                            data = aData,
                            maxit = baseIterationNumber, 
                            method = 2L))
    } else if (family == "bayes.negative.binomial") {
          invisible(glm.nb(as.formula(aFormula), 
                           data = aData, 
                           control = glm.control(maxit = baseIterationNumber), 
                           contrasts = contrasts, 
                           method = ".bayesglm_fit"))
    } else if (family == "bayes.poisson") {
          invisible(bayesglm(as.formula(aFormula), 
                         data = aData, 
                         family = poisson(link = "log"), 
                         control = glm.control(maxit = baseIterationNumber), 
                         prior.scale = prior.scale_parameter, # by default is 5
                         contrasts = contrasts))   
    } else if (family == "zeroinfl.negative.binomial") {
          invisible(zeroinfl(as.formula(paste(aFormula, "| 1")), 
                         data = aData, 
                         control = zeroinfl.control(maxit = baseIterationNumber), 
                         contrasts = contrasts,
                         dist = "negbin"))
    } else if (family == "loglin.poisson") {
          mm <- model.matrix(as.formula(bFormula), data = bData, contrasts.arg = contrasts)
          invisible(fastglm(x = mm,
                          y = bData[[targetVar]],
                          family = poisson(link = "log"),
                          data = bData,
                          maxit = baseIterationNumber,
                          method = 2L))
      
    } else if (family == "loglin.bayes.poisson") {          
          invisible(bayesglm(as.formula(bFormula), 
                             data = bData, 
                             family = poisson(link = "log"), 
                             prior.scale = prior.scale_parameter, # by default is 5
                             control = glm.control(maxit = baseIterationNumber), 
                             contrasts = contrasts))   
    } else {
          stop("Family not implemented.")
    }
  },
  error=function(err){ 
    NULL
  },  
  finally = {}
  )
  
  # Regression did not make it. Return table with NAs
  if (is.null(Glm_model)) {
    result <- data.table(coefName = "(Intercept)", 
                         coef = NA_real_, 
                         `2.5 %` = NA_real_, 
                         `97.5 %` = NA_real_, 
                         pVal = NA_real_, 
                         #Z.Value = NA_real_, 
                         upperTailPval = NA_real_,
                         lowerTailPval = NA_real_, 
                         Std.error = NA_real_, 
                         #drop = NA_real_,
                         converged = FALSE) 
    return(result)
    
  } else {
    
    if (family=="Fisher.OR") {
      result = data.table(coefName = "isTarget1", 
                          coef = log(Glm_model$estimate), 
                          `2.5 %`= log(Glm_model$conf.int[1]), 
                          `97.5 %` = log(Glm_model$conf.int[2]), 
                          pVal = Glm_model$p.value)
      result[, upperTailPval := ifelse(coef>=0, pVal/2, 1-pVal/2)]
      result[, lowerTailPval :=1-upperTailPval]
      if (providePval) return(result)
      else return(result[, -"pVal", with=F] )
    }
    else if (family=="Log.OR") {
      result = data.table(coefName = "isTarget1", 
                          coef = Glm_model$OR, 
                          `2.5 %`= NA_real_, 
                          `97.5 %` = NA_real_, 
                          pVal = Glm_model$pvalue)
      result[, upperTailPval := ifelse(coef>=0, pVal/2, 1-pVal/2)]
      result[, lowerTailPval :=1-upperTailPval]
      if (providePval) return(result)
      else return(result[, -"pVal", with=F] )
      
    } else {
      if (!stringr::str_detect(family,  paste(c("bayes", "zeroinfl"), collapse="|"))) {
        class(Glm_model) <- "fastglm"
      }
      
      
      # Retrieve coefficient and test for significance
      if (stringr::str_detect(family,  paste(c("zeroinfl"), collapse="|"))) {
        temp <- summary(Glm_model)$coefficients$count
      } else {
        temp <- summary(Glm_model)$coefficients
      }
      aCoef = temp[, 1L]
      aStd.error = temp[, 2L]
      aZ.Value = temp[, 3L]
      aPval = temp[, 4L]
      if (colnames(temp)[3L] == "t value") {
        # NB
        aupperTailPval = pt(temp[, 3L], df = Glm_model$df.residual, lower.tail = FALSE)
      } else {
        # Poisson
        aupperTailPval = pnorm(temp[, 3L], lower.tail = FALSE)
      }
      alowerTailPval = 1 - aupperTailPval
      
      
      result = data.table(coefName = rownames(temp), 
                          coef = aCoef, 
                          `2.5 %`= NA_real_, 
                          `97.5 %` = NA_real_, 
                          pVal = aPval, 
                          #Z.Value = aZ.Value, 
                          upperTailPval = aupperTailPval,
                          lowerTailPval = alowerTailPval, 
                          Std.error = aStd.error,
                          #drop = NA_real_,
                          converged = Glm_model$converged)
      rm(temp)
      if (providePval) return(result)
      else return(result[, -"pVal", with=F] )  # if needed for compatibility reasons, remove the pVal
      
    }
  }
}
