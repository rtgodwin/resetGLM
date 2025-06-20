% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resetlm.R
\name{resetlm}
\alias{resetlm}
\title{Lagrange Multiplier RESET Test for glm or fixest Models}
\usage{
resetlm(mod, data, aug.terms = 4, robust = T, fourier = T, sin.link = F)
}
\arguments{
\item{mod}{A model object of class \code{"glm"} or \code{"fixest"}.}

\item{data}{A data frame containing the relevant model data.}

\item{aug.terms}{An integer (1–4) specifying the number of augmenting terms (multiples of 2 or 4 for Fourier).}

\item{robust}{Logical. Whether to use the robust LM test.}

\item{fourier}{Logical. If \code{TRUE}, Fourier series terms are used for augmentation.}

\item{sin.link}{Logical. If using Fourier terms, determines the transformation of the linear predictor.}
}
\value{
A list with components:
\item{reset}{The LM test statistic.}
\item{pval}{The associated p-value.}
}
\description{
Performs a RESET Lagrange Multiplier test for model misspecification in glm or fixest models,
using Taylor or Fourier series augmentation with optional robust variance assumption.
}
\details{
The test evaluates the null that the model is correctly specified. The robust version should almost
always be used. Simulations suggest that the Fourier approximation with linear link and four
augmenting terms work best, and is the current default.
}
