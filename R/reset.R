#' Run Multiple Variants of the RESET
#'
#' Applies a suite of RESET-style specification tests to a glm or fixest model,
#' covering both Wald and Lagrange Multiplier (LM) tests, with options for robust
#' variance assumptions, Taylor or Fourier series augmentations, and combinations
#' thereof.
#'
#' @param mod A model object of class \code{"glm"} or \code{"fixest"}.
#' @param data A data frame containing the variables used in the model.
#'
#' @return A named list of test results. Each element is itself a list with entries:
#'   \item{reset}{The test statistic.}
#'   \item{pval}{The corresponding p-value.}
#'
#' @details
#' This function calls \code{resetWald}, \code{resetWaldTF}, \code{resetlm}, and
#' \code{resetlmTF} to compute variants of the RESET test under different assumptions
#' and augmentation schemes.
#' 
#' @importFrom stats predict formula update as.formula glm lm pchisq poisson vcov
#' @importFrom sandwich vcovHC
#' @importFrom lmtest waldtest
#' @importFrom fixest feglm
#' 
#' @examples
#' # Perform all RESETs on the gravity model from Silva and Tenreyro (2006)
#' data <- read.csv("https://rtgodwin.com/data/gravity.csv")
#' formula <- trade ~ lypex + lypim + lyex + lyim + ldist + border + comlang + 
#'   colony + landl_ex + landl_im + lremot_ex + lremot_im + comfrt_wto + open_wto
#' gravity <- glm(formula, data, family = stats::quasipoisson(link = "log"))
#' reset(gravity, data)
#' 
#' # Perform the recommended LM test on a gravity model with fixed effects and clustering
#' femod <- fixest::feglm(Euros ~ log(dist_km) | Origin + Destination + Year, 
#'   data = fixest::trade, cluster = c("Origin", "Destination"), family = poisson())
#' resetlm(femod, fixest::trade)
#' 
#' @references
#' Silva, J. M. C. Santos and Tenreyro, S. (2006).
#' “The Log of Gravity.”
#' *The Review of Economics and Statistics*, 88(4), 641–658.
#' \doi{10.1162/rest.88.4.641}
#' 
#' @export
#'
#' @seealso 
#'   \code{\link{resetWald}}, 
#'   \code{\link{resetlm}}, 
#'   \code{\link{resetWaldTF}}, 
#'   \code{\link{resetlmTF}}

reset <- function(mod, data) {
  z <- list()
  
  for (test_fun in c("resetWald", "resetlm")) {
    for (robust in c(FALSE, TRUE)) {
      for (fourier in c(FALSE, TRUE)) {
        for (aug.terms in 1:4) {
          
          # Skip invalid Fourier aug.terms (must be even)
          if (fourier && aug.terms %% 2 != 0) next
          
          for (sin.link in if (fourier) c(TRUE, FALSE) else NA) {
            
            tag <- paste0(
              if (test_fun == "resetWald") "Wald" else "lm", ".",
              if (robust) "robust" else "GLM", ".",
              if (fourier) paste0("Fourier.", if (sin.link) "sinlink" else "linlink") else "Taylor",
              ".aug", aug.terms
            )
            
            z[[tag]] <- tryCatch(
              do.call(test_fun, list(
                mod = mod,
                data = data,
                aug.terms = aug.terms,
                robust = robust,
                fourier = fourier,
                sin.link = if (!is.na(sin.link)) sin.link else NULL
              )),
              warning = function(w) "GLM failed to converge",
              error = function(e) "GLM failed to converge"
            )
          }
        }
      }
    }
  }
  
  # Combined Fourier-Taylor tests (only aug.terms = 4 used internally)
  for (test_fun in c("resetWaldTF", "resetlmTF")) {
    for (robust in c(FALSE, TRUE)) {
      for (sin.link in c(TRUE, FALSE)) {
        tag <- paste0(
          if (test_fun == "resetWaldTF") "Wald" else "lm", ".",
          if (robust) "robust" else "GLM",
          ".combined.", if (sin.link) "sinlink" else "linlink"
        )
        
        z[[tag]] <- tryCatch(
          do.call(test_fun, list(
            mod = mod,
            data = data,
            robust = robust,
            sin.link = sin.link
          )),
          warning = function(w) "GLM failed to converge",
          error = function(e) "GLM failed to converge"
        )
      }
    }
  }
  
  return(z)
}
