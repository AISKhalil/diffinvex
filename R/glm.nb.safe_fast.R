#' Fit Generalized Linear Models Using FastGLM
#'
#' This function fits a generalized linear model (GLM) using the `fastglmPure` function, which is designed for efficient GLM fitting. The function provides an interface to handle various model parameters and control settings.
#'
#' @param x A numeric matrix or data frame of predictor variables.
#' @param y A numeric vector or factor of response variables. Must be compatible with the specified family.
#' @param w A numeric vector of weights for each observation. Defaults to a vector of ones.
#' @param family A description of the error distribution and link function to be used in the model. Default is `gaussian()`.
#' @param offset A numeric vector of offset values to be included in the model. Defaults to a vector of zeros.
#' @param start Optional starting values for the coefficients.
#' @param mustart Optional starting values for the mean of the response variable.
#' @param etastart Optional starting values for the linear predictor.
#' @param control A list of control parameters. 
#'   - `maxit`: Maximum number of iterations for the algorithm (default is 100).
#'   - `epsilon`: Convergence tolerance for the algorithm (default is 1e-07).
#' @param intercept Logical flag indicating whether an intercept should be included in the model. (Note: `intercept` is not used in the current implementation but might be intended for future use or included in `fastglmPure`).
#'
#' @return An object returned by `fastglmPure`, typically including fitted coefficients, convergence status, and other model-related information.
#'
#' @details
#' This function leverages the `fastglmPure` function for fitting GLMs, which is optimized for speed and efficiency. It allows for customization of the model fitting process through parameters like weights, offsets, and starting values. The `control` argument allows tuning the convergence criteria and maximum iterations.
#'
#' @importFrom fastglm fastglmPure
#' @export
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
    fastglmPure(x = x,
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



#' Fit Bayesian Generalized Linear Models
#'
#' This function fits a Bayesian generalized linear model (GLM) using the `bayesglm.fit` function, with specific prior settings for the model. The function provides an interface to customize various parameters related to the GLM fitting process.
#'
#' @param x A numeric matrix or data frame of predictor variables.
#' @param y A numeric vector or factor of response variables. Must be compatible with the specified family.
#' @param w A numeric vector of weights for each observation. Defaults to a vector of ones.
#' @param family A description of the error distribution and link function to be used in the model. Default is `gaussian()`.
#' @param offset A numeric vector of offset values to be included in the model. Defaults to a vector of zeros.
#' @param start Optional starting values for the coefficients.
#' @param mustart Optional starting values for the mean of the response variable.
#' @param etastart Optional starting values for the linear predictor.
#' @param control A list of control parameters.
#'   - `maxit`: Maximum number of iterations for the algorithm (default is 100).
#'   - `epsilon`: Convergence tolerance for the algorithm (default is 1e-07).
#' @param intercept Logical flag indicating whether an intercept should be included in the model. (Note: `intercept` is not used in the current implementation but might be intended for future use or included in `bayesglm.fit`).
#'
#' @return An object returned by `bayesglm.fit`, typically including fitted coefficients, convergence status, and other model-related information.
#'
#' @details
#' This function leverages the `bayesglm.fit` function for fitting Bayesian GLMs. It includes default prior settings for the model, with a prior scale of 1 and a prior mean for the intercept of -14. The `control` argument allows tuning the convergence criteria and maximum iterations.
#'
#' @importFrom bayesglm bayesglm.fit
#' @export
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
  bayesglm.fit(x=x,
              y=y,
              weights = w,
              start = start,
              etastart = etastart,
              mustart = mustart,
              offset = offset,
              family = family,
              control = control,
              prior.scale = 1, # default value = 5
              prior.mean.for.intercept = -14
  )
}



#' Fit Generalized Linear Models with Various Families
#'
#' This function fits a generalized linear model (GLM) using various families and link functions. It supports multiple types of models including negative binomial, Poisson, Bayesian models, and specialized tests like Fisher's exact test and Log Odds Ratio. The function also allows for including interaction terms and provides options for p-values and confidence intervals.
#'
#' @param data A data frame containing the variables used in the model.
#' @param family A character string specifying the family of the GLM. Valid options include "negative.binomial", "poisson", "bayes.negative.binomial", "bayes.poisson", "zeroinfl.negative.binomial", "loglin.poisson", "loglin.bayes.poisson", "Fisher.OR", and "Log.OR".
#' @param targetVar A string specifying the name of the target (response) variable.
#' @param offsetVar A string specifying the name of the offset variable. Only used for certain families.
#' @param featureVars A vector of strings specifying the names of the feature variables to include in the model.
#' @param interVars A vector of strings specifying the names of interaction variables to include. Optional.
#' @param providePval A logical indicating whether to include p-values in the output. Defaults to `FALSE`.
#' @param contrasts A list of contrasts to apply to the categorical variables in the model.
#' @param provideConfint A logical indicating whether to provide confidence intervals in the output. Defaults to `FALSE`.
#' @param baseIterationNumber An integer specifying the maximum number of iterations for fitting the model. Defaults to 100.
#'
#' @return A data.table with the model coefficients and associated statistics. The exact contents depend on the family used:
#' \itemize{
#'   \item For "Fisher.OR", the result includes the odds ratio and p-value.
#'   \item For "Log.OR", the result includes the log odds ratio and p-value.
#'   \item For other families, the result includes coefficients, standard errors, p-values, and confidence intervals if requested.
#' }
#'
#' @details
#' The function fits various GLMs based on the specified family:
#' \itemize{
#'   \item **Negative Binomial and Bayesian Negative Binomial**: Uses `glm.nb` from the `MASS` package with options for the number of iterations and control parameters.
#'   \item **Poisson**: Uses `fastglm` for efficient GLM fitting.
#'   \item **Bayesian Poisson**: Uses `bayesglm` with default prior settings.
#'   \item **Zero-Inflated Negative Binomial**: Uses `zeroinfl` from the `pscl` package.
#'   \item **Log-Linear Poisson and Bayesian Log-Linear Poisson**: Uses `fastglm` and `bayesglm` respectively.
#'   \item **Fisher's Exact Test and Log Odds Ratio**: Custom implementations for contingency table analysis.
#' }
#'
#' @importFrom data.table data.table
#' @importFrom dplyr `%>%` 
#' @importFrom MASS glm.nb
#' @importFrom pscl zeroinfl
#' @importFrom fastglm fastglm
#' @importFrom bayesglm bayesglm
#' @export
glm.safe <- function(data,
                    family,
                    targetVar,
                    offsetVar,
                    featureVars,
                    interVars = NULL,
                    providePval = FALSE,
                    contrasts,
                    provideConfint = FALSE,
                    baseIterationNumber = 100) {
  #
  prior.scale_parameter <- 1 # default value = 5
  #
  if  (! family %in% c("loglin.poisson", "loglin.bayes.poisson", "Fisher.OR", "Log.OR")) {
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
