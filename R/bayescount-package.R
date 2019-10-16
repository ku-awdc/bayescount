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
#' \emph{Data analysis}
#'
#' See the analyse_fecrt function
#' 
#' \emph{Sample size calculations}
#' 
#' TODO
#' 
#' \emph{Publication results}
#' 
#' Code is provided to replicate the simulation results in section 6 of "A hypothesis testing framework for the ratio of means of two negative binomial distributions" via the "htfsim" vignette:  vignette("htfsim", package="bayescount")
#' 
#' @references
#' M. J. Denwood, S. W. J. Reid, S. Love, M. K. Nielsen, L. Matthews, I. J. McKendrick, and G. T. Innocent. Comparison of three alternative methods for analysis of equine Faecal Egg Count Reduction Test data. Prev. Vet. Med. (2009), doi:10.1016/j.prevetmed.2009.11.009
#' 
#' Denwood, M. J. (2010). A quantitative approach to improving the analysis of faecal worm egg count data. University of Glasgow. Retrieved from http://www.gla.ac.uk/media/media_149338_en.pdf
#' 
#' @seealso \code{\link{launch_shiny}}, \code{\link{analyse_fecrt}}
#' @import Rcpp
#' @useDynLib bayescount
NULL
