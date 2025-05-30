#' Lagrange Multiplier RESET Test for glm or fixest Models
#'
#' Performs a RESET Lagrange Multiplier test for model misspecification in glm or fixest models,
#' using Taylor or Fourier series augmentation with optional robust variance assumption.
#'
#' @param mod A model object of class \code{"glm"} or \code{"fixest"}.
#' @param data A data frame containing the relevant model data.
#' @param aug.terms An integer (1–4) specifying the number of augmenting terms (multiples of 2 or 4 for Fourier).
#' @param robust Logical. Whether to use the robust LM test.
#' @param fourier Logical. If \code{TRUE}, Fourier series terms are used for augmentation.
#' @param sin.link Logical. If using Fourier terms, determines the transformation of the linear predictor.
#'
#' @return A list with components:
#'   \item{reset}{The LM test statistic.}
#'   \item{pval}{The associated p-value.}
#'
#' @details
#' The test evaluates the null that the model is correctly specified. The robust version should almost
#' always be used. Simulations suggest that the Fourier approximation with linear link and four
#' augmenting terms work best, and is the current default.
#'
#' @export

resetlm <- function(mod, data, aug.terms = 4, robust = T, fourier = T, sin.link = F) {
  
  if (!inherits(mod, "glm") && !inherits(mod, "fixest")) {
    stop("Error: 'mod' must be a 'glm' or 'fixest' object.")
  }
  
  yhat <- mod$fitted.values
  etahat <- predict(mod, type = "link")
  y <- mod$y
  uhat <- y - yhat
  w <- sqrt(yhat)
  util <- uhat / w
  
  if (!is.null(mod$obs_selection$obsRemoved) && length(mod$obs_selection$obsRemoved) > 0) {
    all_rows <- 1:nrow(data)
    removed_rows <- abs(mod$obs_selection$obsRemoved)
    used_rows <- setdiff(all_rows, removed_rows)
  } else {
    used_rows <- 1:nrow(data)
  }
  
  data$etahat <- data$w <- data$util <- NA_real_
  data$aux1resids <- data$aux2resids <- data$aux3resids <- data$aux4resids <- data$ones <- NA_real_
  
  data$etahat[used_rows] <- etahat
  data$w[used_rows] <- w
  data$util[used_rows] <- util
  
    n <- length(used_rows)
  z <- list()
  
  if (inherits(mod, "glm")) {
    original_formula <- mod$formula
    response_var <- all.vars(original_formula)[1]
    predictor_vars <- all.vars(original_formula)[-1]
  }
  
  if (inherits(mod, "fixest")) {
    rhs <- strsplit(strsplit(paste(deparse(formula(mod)), collapse = " "), "\\|")[[1]][1], "~")[[1]][2]
    rhs <- trimws(rhs)
    predictor_vars <- trimws(strsplit(rhs, "\\+")[[1]])
  }
  
  if(fourier == T) {
    
    if(sin.link == T) {
      data$v <- 2 * pi * ((sin(data$etahat)) ^ 2) - pi
    }
    
    if(sin.link == F) {
      data$v <- pi * (2 * (data$etahat) - (max(etahat) + min(etahat))) / (max(etahat) - min(etahat))
    }
    
    if(robust == T) {
      
      if(aug.terms == 2) {
        new_response_var1 <- "I(w * sin(v))"
        new_formula_string1 <- paste(new_response_var1, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula1 <- as.formula(new_formula_string1)
        data$aux1resids[used_rows] <- lm(new_formula1, data = data)$residuals
        
        new_response_var2 <- "I(w * cos(v))"
        new_formula_string2 <- paste(new_response_var2, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula2 <- as.formula(new_formula_string2)
        data$aux2resids[used_rows] <- lm(new_formula2, data = data)$residuals
        
        data$ones[used_rows] <- rep(1, n)
        aux <- lm(ones ~ I(util * aux1resids) + I(util * aux2resids) -1, data = data)
        
        z$reset <- n - sum(aux$residuals ^ 2)
        z$pval <- 1 - pchisq(z$reset, aug.terms)
      }
      
      if(aug.terms == 4) {
        new_response_var1 <- "I(w * sin(v))"
        new_formula_string1 <- paste(new_response_var1, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula1 <- as.formula(new_formula_string1)
        data$aux1resids[used_rows] <- lm(new_formula1, data = data)$residuals
        
        new_response_var2 <- "I(w * cos(v))"
        new_formula_string2 <- paste(new_response_var2, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula2 <- as.formula(new_formula_string2)
        data$aux2resids[used_rows] <- lm(new_formula2, data = data)$residuals
        
        new_response_var3 <- "I(w * sin(2*v))"
        new_formula_string3 <- paste(new_response_var3, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula3 <- as.formula(new_formula_string3)
        data$aux3resids[used_rows] <- lm(new_formula3, data = data)$residuals
        
        new_response_var4 <- "I(w * cos(2*v))"
        new_formula_string4 <- paste(new_response_var4, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula4 <- as.formula(new_formula_string4)
        data$aux4resids[used_rows] <- lm(new_formula4, data = data)$residuals
        
        data$ones[used_rows] <- rep(1, n)
        aux <- lm(ones ~ I(util * aux1resids) + I(util * aux2resids) + I(util * aux3resids) + I(util * aux4resids) -1, data = data)
        
        z$reset <- n - sum(aux$residuals ^ 2)
        z$pval <- 1 - pchisq(z$reset, aug.terms)
      }
    }
    
    if(robust == F) {
      new_response_var <- "util"
      
      if(aug.terms == 2) {
        new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * sin(v)) + I(w * cos(v)) -1")
        new_formula <- as.formula(new_formula_string)
        
        aux <- lm(new_formula, data = data)
        
        z$reset <- n * summary(aux)$r.squared
        z$pval <- 1 - pchisq(z$reset, aug.terms)
      }
      
      if(aug.terms == 4) {
        new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * sin(v)) + I(w * cos(v)) + I(w * sin(2*v)) + I(w * cos(2*v)) -1")
        new_formula <- as.formula(new_formula_string)
        
        aux <- lm(new_formula, data = data)
        
        z$reset <- n * summary(aux)$r.squared
        z$pval <- 1 - pchisq(z$reset, aug.terms)
      }
    }
  }
  
  if(fourier == F) {
    
    if(robust == T) {
      
      if(aug.terms == 1) {
        new_response_var1 <- "I(w * etahat ^ 2)"
        new_formula_string1 <- paste(new_response_var1, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula1 <- as.formula(new_formula_string1)
        data$aux1resids[used_rows] <- lm(new_formula1, data = data)$residuals
        
        data$ones[used_rows] <- rep(1, n)
        aux <- lm(ones ~ I(util * aux1resids) -1, data = data)
        
        z$reset <- n - sum(aux$residuals ^ 2)
        z$pval <- 1 - pchisq(z$reset, aug.terms)
      }
      
      if(aug.terms == 2) {
        new_response_var1 <- "I(w * etahat ^ 2)"
        new_formula_string1 <- paste(new_response_var1, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula1 <- as.formula(new_formula_string1)
        data$aux1resids[used_rows] <- lm(new_formula1, data = data)$residuals
        
        new_response_var2 <- "I(w * etahat ^ 3)"
        new_formula_string2 <- paste(new_response_var2, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula2 <- as.formula(new_formula_string2)
        data$aux2resids[used_rows] <- lm(new_formula2, data = data)$residuals
        
        data$ones[used_rows] <- rep(1, n)
        aux <- lm(ones ~ I(util * aux1resids) + I(util * aux2resids) -1, data = data)
        
        z$reset <- n - sum(aux$residuals ^ 2)
        z$pval <- 1 - pchisq(z$reset, aug.terms)
      }
      
      if(aug.terms == 3) {
        new_response_var1 <- "I(w * etahat ^ 2)"
        new_formula_string1 <- paste(new_response_var1, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula1 <- as.formula(new_formula_string1)
        data$aux1resids[used_rows] <- lm(new_formula1, data = data)$residuals
        
        new_response_var2 <- "I(w * etahat ^ 3)"
        new_formula_string2 <- paste(new_response_var2, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula2 <- as.formula(new_formula_string2)
        data$aux2resids[used_rows] <- lm(new_formula2, data = data)$residuals
        
        new_response_var3 <- "I(w * etahat ^ 4)"
        new_formula_string3 <- paste(new_response_var3, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula3 <- as.formula(new_formula_string3)
        data$aux3resids[used_rows] <- lm(new_formula3, data = data)$residuals
        
        data$ones[used_rows] <- rep(1, n)
        aux <- lm(ones ~ I(util * aux1resids) + I(util * aux2resids) + I(util * aux3resids) -1, data = data)
        
        z$reset <- n - sum(aux$residuals ^ 2)
        z$pval <- 1 - pchisq(z$reset, aug.terms)
      }
      
      if(aug.terms == 4) {
        new_response_var1 <- "I(w * etahat ^ 2)"
        new_formula_string1 <- paste(new_response_var1, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula1 <- as.formula(new_formula_string1)
        data$aux1resids[used_rows] <- lm(new_formula1, data = data)$residuals
        
        new_response_var2 <- "I(w * etahat ^ 3)"
        new_formula_string2 <- paste(new_response_var2, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula2 <- as.formula(new_formula_string2)
        data$aux2resids[used_rows] <- lm(new_formula2, data = data)$residuals
        
        new_response_var3 <- "I(w * etahat ^ 4)"
        new_formula_string3 <- paste(new_response_var3, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula3 <- as.formula(new_formula_string3)
        data$aux3resids[used_rows] <- lm(new_formula3, data = data)$residuals
        
        new_response_var4 <- "I(w * etahat ^ 5)"
        new_formula_string4 <- paste(new_response_var4, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "- 1")
        new_formula4 <- as.formula(new_formula_string4)
        data$aux4resids[used_rows] <- lm(new_formula4, data = data)$residuals
        
        data$ones[used_rows] <- rep(1, n)
        aux <- lm(ones ~ I(util * aux1resids) + I(util * aux2resids) + I(util * aux3resids) + I(util * aux4resids) -1, data = data)
        
        z$reset <- n - sum(aux$residuals ^ 2)
        z$pval <- 1 - pchisq(z$reset, aug.terms)
      }
    }
    
    if(robust == F) {
      new_response_var <- "util"
      if(aug.terms == 1) {
        new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * etahat ^ 2) - 1")
      } else if(aug.terms == 2) {
        new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * etahat ^ 2) + I(w * etahat ^ 3) - 1")
      } else if(aug.terms == 3) {
        new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * etahat ^ 2) + I(w * etahat ^ 3) + I(w * etahat ^ 4) - 1")
      } else if(aug.terms == 4) {
        new_formula_string <- paste(new_response_var, "~ w +", paste("I(w *", predictor_vars, ")", collapse=" + "), "+ I(w * etahat ^ 2) + I(w * etahat ^ 3) + I(w * etahat ^ 4) + I(w * etahat ^ 5) - 1")
      }
      new_formula <- as.formula(new_formula_string)
      
      aux <- lm(new_formula, data = data)
      
      z$reset <- n * summary(aux)$r.squared
      z$pval <- 1 - pchisq(z$reset, aug.terms)
    }
  }
  return(z)
}
