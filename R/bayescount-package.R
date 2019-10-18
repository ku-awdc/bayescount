#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL

#' @name bayescount-package
#' @aliases bayescountpackage bayescount
#' @docType package
#' @title Statistical Analyses and Sample Size Calculations for Over-dispersed Count Data
#'
#' @description
#'
#' A full description of the main functions in this package is given in the package vignette:
#' vignette("bayescount", package="bayescount")
#'
#' For quick introductions to the two main functions of the package, see:
#'
#' \emph{Data analysis}
#'
#' See ?efficacy_analysis
#'
#' \emph{Sample size calculations}
#'
#' See ?efficacy_frequencies
#'
#'
#' @references
#' M. J. Denwood, S. W. J. Reid, S. Love, M. K. Nielsen, L. Matthews, I. J. McKendrick, and G. T. Innocent. Comparison of three alternative methods for analysis of equine Faecal Egg Count Reduction Test data. Prev. Vet. Med. (2009), doi:10.1016/j.prevetmed.2009.11.009
#'
#' Denwood, M. J. (2010). A quantitative approach to improving the analysis of faecal worm egg count data. University of Glasgow. Retrieved from http://www.gla.ac.uk/media/media_149338_en.pdf
#'
#' @seealso \code{\link{launch_shiny}}, \code{\link{efficacy_analysis}}, \code{\link{efficacy_frequencies}}
#' @import Rcpp
#' @useDynLib bayescount
NULL
