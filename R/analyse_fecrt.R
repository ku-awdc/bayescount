#' analyse_fecrt
#'
#' @description
#' This function analyses a given dataset using different methods
#'
#' @param data_1
#' @param data_2
#' @param H0_I
#' @param H0_A
#' @param paired
#' @param tail
#' @param mean_ratio
#' @param conjugate_priors
#' @param delta
#' @param beta_iters
#' @param approx
#' @param dobson_cl
#' @param dobson_priors
#'
#' @return
#' @export
#'
#' @examples
#' data1 <- rnbinom(10, mu=20, size=1)
#' data2 <- rnbinom(10, mu=2, size=0.75)
#' analyse_fecrt(data1, data2, H0_I=0.95, H0_A=0.9, paired=FALSE)
analyse_fecrt <- function(data_1, data_2, H0_I, H0_A, paired, tail=0.025, mean_ratio=1, conjugate_priors=c(0,0), delta=1, beta_iters=10^4, approx=1, dobson_cl=c(0.005,0.995), dobson_priors=c(1,1)) {
	
	# TODO: input checks
	
   # .Call('_bayescount_analyse_fecrt', PACKAGE = 'bayescount', as.integer(data_1), as.integer(data_2), as.logical(paired), as.double(mean_ratio), as.double(H0_I), as.double(H0_A), as.double(conjugate_priors), as.integer(delta), as.integer(beta_iters), as.integer(approx), as.double(tail), as.double(dobson_cl), as.double(dobson_priors))
   
   RCPP_analyse_fecrt(as.integer(data_1), as.integer(data_2), as.logical(paired), as.double(mean_ratio), as.double(H0_I), as.double(H0_A), as.double(conjugate_priors), as.integer(delta), as.integer(beta_iters), as.integer(approx), as.double(tail), as.double(dobson_cl), as.double(dobson_priors))
   
}

