#' Combined Wald RESET Test with Taylor and Fourier Augmentation
#'
#' Computes a four-term Wald RESET test combining both Fourier and Taylor series
#' augmentation terms. Applies to glm or fixest models.
#'
#' @param mod A model object of class \code{"glm"} or \code{"fixest"}.
#' @param data A data frame containing the model variables.
#' @param robust Logical. Whether to use the robust LM test.
#' @param sin.link Logical. Controls transformation used in the Fourier series terms.
#'
#' @return A list containing:
#'   \item{reset}{The Wald test statistic.}
#'   \item{pval}{The associated p-value.}
#'   
#' @details
#' The test evaluates the null that the model is correctly specified. Simulations
#' suggest that the LM versions are superior.
#' 
#' @export

resetWaldTF <- function(mod, data, robust = T, sin.link = F) {
  
  if (!inherits(mod, "glm") && !inherits(mod, "fixest")) {
    stop("Error: 'mod' must be a 'glm' or 'fixest' object.")
  }
  
  etahat <- predict(mod, type = "link")
  data$etahat <- NA_real_
  data$etahat[seq_along(etahat)] <- etahat
  
  z <- list()
  aug.terms <- 4
  
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
  
  if (inherits(mod, "glm")) {
    original_formula <- formula(mod)
  }
  
  if (inherits(mod, "fixest")) {
    parts <- strsplit(paste(deparse(formula(mod)), collapse = " "), "\\|")[[1]]
    parts <- trimws(sub("\\+$", "", gsub("\\s+", " ", parts)))
    rhs_main <- parts[1]
    fes_rhs <- parts[2]
    
    call_text <- paste(deparse(mod$call), collapse = " ")
    cluster_rhs <- sub(".*cluster *= *(c\\([^\\)]*\\)).*", "\\1", call_text)
    cluster_vars <- eval(parse(text = gsub('\\"', '"', cluster_rhs)))
  }
  
  if(sin.link == T) {
    data$v <- 2 * pi * ((sin(data$etahat)) ^ 2) - pi
  }
  
  if(sin.link == F) {
    data$v <- pi * (2 * (data$etahat) - (max(etahat) + min(etahat))) / (max(etahat) - min(etahat))
  }
  
  aug_terms <- "sin(v) + cos(v) + I(etahat^2) + I(etahat^3)"
  
  if (inherits(mod, "glm")) {
    new_formula <- update(original_formula, paste(". ~ . +", aug_terms))
    aux <- glm(new_formula, data = data, family = stats::quasipoisson(link = "log"))
  }
  
  if (inherits(mod, "fixest")) {
    rhs_aug <- paste(rhs_main, aug_terms, sep = " + ")
    new_formula <- as.formula(paste(rhs_aug, "|", fes_rhs))
    aux <- feglm(new_formula, data = data, cluster = cluster_vars, family = poisson())
  }
  
  temp_result <- try({
    if (robust) {
      vcov_fn <- if (inherits(mod, "glm")) {
        function(x) sandwich::vcovHC(x, type = "HC0")
      } else if (inherits(mod, "fixest")) {
        vcov
      }
      lmtest::waldtest(mod, aux, test = "Chisq", vcov = vcov_fn(aux))
    } else if (robust == F) {
      lmtest::waldtest(mod, aux, test = "Chisq")
    }
  }, silent = TRUE)
  
  if (inherits(temp_result, "try-error")) {
    z$reset <- NA_real_
    z$pval <- NA_real_
  } else {
    z$reset <- temp_result$Chisq[2]
    z$pval <- temp_result$"Pr(>Chisq)"[2]
  }
  
  return(z)
}