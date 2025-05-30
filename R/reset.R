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
  
  #GLM versions
  
  #Wald
  z$Wald.GLM.Taylor.aug1 <- resetWald(mod, data, aug.terms = 1, robust = F, fourier = F)
  z$Wald.GLM.Taylor.aug2 <- resetWald(mod, data, aug.terms = 2, robust = F, fourier = F)
  z$Wald.GLM.Taylor.aug3 <- tryCatch(
    resetWald(mod, data, aug.terms = 3, robust = F, fourier = F),
    warning = function(w) {z$Wald.GLM.Taylor.aug3 <<- "GLM failed to converge"})
  z$Wald.GLM.Taylor.aug4 <- resetWald(mod, data, aug.terms = 4, robust = F, fourier = F)
  
  z$Wald.GLM.Fourier.sinlink.aug2 <- resetWald(mod, data, aug.terms = 2, robust = F, fourier = T, sin.link = T)
  z$Wald.GLM.Fourier.sinlink.aug4 <- resetWald(mod, data, aug.terms = 4, robust = F, fourier = T, sin.link = T)
  
  z$Wald.GLM.Fourier.linlink.aug2 <- resetWald(mod, data, aug.terms = 2, robust = F, fourier = T, sin.link = F)
  z$Wald.GLM.Fourier.linlink.aug4 <- resetWald(mod, data, aug.terms = 4, robust = F, fourier = T, sin.link = F)
  
  z$Wald.GLM.combined.sinlink <- resetWaldTF(mod, data, robust = F, sin.link = T)
  z$Wald.GLM.combined.linlink <- resetWaldTF(mod, data, robust = F, sin.link = F)
  
  #LM
  z$lm.GLM.Taylor.aug1 <- resetlm(mod, data, aug.terms = 1, robust = F, fourier = F)
  z$lm.GLM.Taylor.aug2 <- resetlm(mod, data, aug.terms = 2, robust = F, fourier = F)
  z$lm.GLM.Taylor.aug3 <- resetlm(mod, data, aug.terms = 3, robust = F, fourier = F)
  z$lm.GLM.Taylor.aug4 <- resetlm(mod, data, aug.terms = 4, robust = F, fourier = F)
  
  z$lm.GLM.Fourier.sinlink.aug2 <- resetlm(mod, data, aug.terms = 2, robust = F, fourier = T, sin.link = T)
  z$lm.GLM.Fourier.sinlink.aug4 <- resetlm(mod, data, aug.terms = 4, robust = F, fourier = T, sin.link = T)
  
  z$lm.GLM.Fourier.linlink.aug2 <- resetlm(mod, data, aug.terms = 2, robust = F, fourier = T, sin.link = F)
  z$lm.GLM.Fourier.linlink.aug4 <- resetlm(mod, data, aug.terms = 4, robust = F, fourier = T, sin.link = F)
  
  z$lm.GLM.combined.sinlink <- resetlmTF(mod, data, robust = F, sin.link = T)
  z$lm.GLM.combined.linlink <- resetlmTF(mod, data, robust = F, sin.link = F)
  
  #robust versions
  
  #Wald
  z$Wald.robust.Taylor.aug1 <- resetWald(mod, data, aug.terms = 1, robust = T, fourier = F)
  z$Wald.robust.Taylor.aug2 <- resetWald(mod, data, aug.terms = 2, robust = T, fourier = F)
  z$Wald.robust.Taylor.aug3 <- resetWald(mod, data, aug.terms = 3, robust = T, fourier = F)
  z$Wald.robust.Taylor.aug4 <- resetWald(mod, data, aug.terms = 4, robust = T, fourier = F)
  
  z$Wald.robust.Fourier.sinlink.aug2 <- resetWald(mod, data, aug.terms = 2, robust = T, fourier = T, sin.link = T)
  z$Wald.robust.Fourier.sinlink.aug4 <- resetWald(mod, data, aug.terms = 4, robust = T, fourier = T, sin.link = T)
  
  z$Wald.robust.Fourier.linlink.aug2 <- resetWald(mod, data, aug.terms = 2, robust = T, fourier = T, sin.link = F)
  z$Wald.robust.Fourier.linlink.aug4 <- resetWald(mod, data, aug.terms = 4, robust = T, fourier = T, sin.link = F)
  
  z$Wald.robust.combined.sinlink <- resetWaldTF(mod, data, robust = T, sin.link = T)
  z$Wald.robust.combined.linlink <- resetWaldTF(mod, data, robust = T, sin.link = F)
  
  #LM
  z$lm.robust.Taylor.aug1 <- resetlm(mod, data, aug.terms = 1, robust = T, fourier = F)
  z$lm.robust.Taylor.aug2 <- resetlm(mod, data, aug.terms = 2, robust = T, fourier = F)
  z$lm.robust.Taylor.aug3 <- resetlm(mod, data, aug.terms = 3, robust = T, fourier = F)
  z$lm.robust.Taylor.aug4 <- resetlm(mod, data, aug.terms = 4, robust = T, fourier = F)
  
  z$lm.robust.Fourier.sinlink.aug2 <- resetlm(mod, data, aug.terms = 2, robust = T, fourier = T, sin.link = T)
  z$lm.robust.Fourier.sinlink.aug4 <- resetlm(mod, data, aug.terms = 4, robust = T, fourier = T, sin.link = T)
  
  z$lm.robust.Fourier.linlink.aug2 <- resetlm(mod, data, aug.terms = 2, robust = T, fourier = T, sin.link = F)
  z$lm.robust.Fourier.linlink.aug4 <- resetlm(mod, data, aug.terms = 4, robust = T, fourier = T, sin.link = F)
  
  z$lm.robust.combined.sinlink <- resetlmTF(mod, data, robust = T, sin.link = T)
  z$lm.robust.combined.linlink <- resetlmTF(mod, data, robust = T, sin.link = F)
  
  return(z)
}
