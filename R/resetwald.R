#' Wald RESET Test for glm or fixest Models
#'
#' Performs a RESET Wald test for model misspecification in glm or fixest models,
#' using Taylor or Fourier series augmentation with optional robust variance assumption.
#'
#' @param mod A model object of class \code{"glm"} or \code{"fixest"}.
#' @param data A data frame containing the relevant model data.
#' @param aug.terms An integer (1–4) specifying the number of augmenting terms (multiples of 2 or 4 for Fourier).
#' @param robust Logical. Whether to use the robust LM test.
#' @param fourier Logical. If \code{TRUE}, Fourier series terms are used for augmentation.
#' @param sin.link Logical. If using Fourier terms, determines the transformation of the linear predictor.
#'
#' @return A list with:
#'   \item{reset}{The Wald test statistic.}
#'   \item{pval}{The associated p-value.}
#'
#' @details
#' The test evaluates the null that the model is correctly specified. Simulations
#' suggest that the LM versions are superior.
#'
#' @export

resetWald<- function(mod, data, aug.terms = 4, robust = T, fourier = T, sin.link = F) {
  
  if (!inherits(mod, "glm") && !inherits(mod, "fixest")) {
    stop("Error: 'mod' must be a 'glm' or 'fixest' object.")
  }
  
  etahat <- predict(mod, type = "link")
  
  if (!is.null(mod$obs_selection$obsRemoved) && length(mod$obs_selection$obsRemoved) > 0) {
    all_rows <- 1:nrow(data)
    removed_rows <- abs(mod$obs_selection$obsRemoved)
    used_rows <- setdiff(all_rows, removed_rows)
  } else {
    used_rows <- 1:nrow(data)
  }
  
  data$etahat <- NA_real_
  data$aux1resids <- data$aux2resids <- data$aux3resids <- data$aux4resids <- data$ones <- NA_real_
  
  data$etahat[used_rows] <- etahat
  
  z <- list()
  
  if (inherits(mod, "glm")) {
    original_formula <- formula(mod)
  }
  
  if (inherits(mod, "fixest")) {
    parts <- strsplit(paste(deparse(formula(mod)), collapse = " "), "\\|")[[1]]
    parts <- trimws(parts)
    parts <- sub("\\+$", "", parts)
    parts <- gsub("\\s+", " ", parts)
    rhs_main <- trimws(parts[1])
    fes_rhs <- trimws(parts[2])
    
    call_text <- deparse(mod$call)
    call_str <- paste(call_text, collapse = " ")
    cluster_rhs <- sub('.*cluster *= *(c\\([^\\)]*\\)).*', '\\1', call_str)
    
    cluster_clean <- gsub('\\"', '"', cluster_rhs)      
    cluster_clean <- gsub(', *', ', ', cluster_clean)    
    cluster_vars <- eval(parse(text = cluster_clean))
  }
  
  if(fourier == T) {
    
    if(sin.link == T) {
      data$v <- 2 * pi * sin(etahat) ^ 2 - pi
    }
    
    if(sin.link == F) {
      data$v <- pi * (2 * (etahat) - (max(etahat) + min(etahat))) / (max(etahat) - min(etahat))
    }
    
    if (inherits(mod, "glm")) {
      if(aug.terms == 1 | aug.terms == 3) {stop("Augmentation terms must be a multiple of 2 for the Fourier transform")}
      if(aug.terms == 2) {new_formula <- update(original_formula, . ~ . + sin(v) + cos(v))}
      if(aug.terms == 4) {new_formula <- update(original_formula, . ~ . + sin(v) + cos(v) + I(sin(2*v)) + I(cos(2*v)))}
    }
    
    if (inherits(mod, "fixest")) {
      if(aug.terms == 1 | aug.terms == 3) {stop("Augmentation terms must be a multiple of 2 for the Fourier transform")}
      if(aug.terms == 2) {aug_terms <- "sin(v) + cos(v)"}
      if(aug.terms == 4) {aug_terms <- "sin(v) + cos(v) + I(sin(2*v)) + I(cos(2*v))"}
      rhs_aug <- paste(rhs_main, aug_terms, sep = " + ")
      new_formula <- as.formula(paste(rhs_aug, "|", fes_rhs))
    }
  }
  
  if(fourier == F) {
    aug_terms <- paste0("I(etahat^", 2:(aug.terms + 1), ")")
    if (inherits(mod, "glm")) {new_formula <- update(original_formula, paste(". ~ . +", paste(aug_terms, collapse = " + ")))}
    if (inherits(mod, "fixest")) {
      rhs_aug <- paste(rhs_main, paste(aug_terms, collapse = " + "), sep = " + ")
      new_formula <- as.formula(paste(rhs_aug, "|", fes_rhs))
    }
  }
  
  if (inherits(mod, "glm")) {
    aux <- glm(new_formula, data = data, family = stats::quasipoisson(link = "log"))
  }
  
  if (inherits(mod, "fixest")) {
    aux <- feglm(new_formula, data = data, cluster = cluster_vars, family = poisson())
  }
  
  if (isTRUE(robust)) {
    vcov_fun <- switch(
      class(mod)[1],
      glm = function(x) sandwich::vcovHC(x, type = "HC0"),
      fixest = vcov,
      stop("Unsupported model class for robust testing")
    )
    temp_result <- try(lmtest::waldtest(mod, aux, test = "Chisq", vcov = vcov_fun(aux)), silent = TRUE)
  } else {
    temp_result <- try(lmtest::waldtest(mod, aux, test = "Chisq"), silent = TRUE)
  }
  
  if (inherits(temp_result, "try-error")) {
    z$reset <- NA_real_
    z$pval <- NA_real_
  } else {
    z$reset <- temp_result$Chisq[2]
    z$pval <- temp_result$"Pr(>Chisq)"[2]
  }
  return(z)
}