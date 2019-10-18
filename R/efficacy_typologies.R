#' Obtain typology classifications relating to a theoretical dataset that may be encountered in practice
#'
#' @description This function provides the typology classifications (for each statistical method chosen) relating to a theoretical dataset with the given properties
#'
#' @param sum a length-2 integer representing the total count in the control/pre-treatment data and the treatment/post-treatment data
#' @param N a length-2 integer representing the number of individuals in the control/pre-treatment data and the treatment/post-treatment data
#' @param k a length-2 strictly positive continuous number giving the over-dispersion of the control/pre-treatment and treatment/post-treatment data samples
#' @param cor a length-1 number between 0 and 1 giving the correlation between pre- and post-treatment data (ignored for an unpaired study)
#' @param paired logical flag for a paired or unpaired study design assumption on which to simulate (and analyse) the data
#' @param T_I the threshold for inferioirty (target efficacy of the intervention)
#' @param T_A the threshold for non-inferioirty (target efficacy of the intervention minus a non-inferioirty margin delta)
#' @param R a length-2 integer giving the number of replicates used for control/pre-treatment and treatment/post-treatment individuals
#' @param S a length-2 numeric variable giving the counting sensitivity used for control/pre-treatment and treatment/post-treatment individuals
#' @param method a character vector indicating the statistical methods to show results for, where 'all' means all available methods
#' @param tail the significance level to use for the classification, where 0.025 corresponds to 95\% CI and 0.05 corresponds to 90\% CI
#' @param bnb_priors priors to use for the BNB method
#' @param use_delta logical flag to use the delta method approximation for the BNB method (NA means to use it unless it fails)
#' @param beta_iters number of iterations to use for the Monte Carlo approximation of the beta distribution transformation for the BNB method (when use_delta==FALSE)
#' @param binomial_priors the priors to use for the Binomial method
#' @param binomial_cl_adj the adjustment in confidence level to use for the Binomial method relative to the others (the default corresponds to e.g. 99\% rather than 95\% CI) - where 0 > adj >= 1
#'
#' @return A data frame as would be returned by \code{\link{efficacy_analysis}} for a theoretical dataset with the given properties
#' @export
#' @seealso \code{\link{bayescount}}, \code{\link{efficacy_frequencies}}
#'
#' @examples
#' efficacy_typologies(sum=c(20,1), N=c(10,10), k=c(1,1), cor=0.25,
#' paired = TRUE, T_I = 0.99, T_A = 0.95)

efficacy_typologies <- function(sum, N, k=c(1,1), cor=0.1, paired=TRUE, T_I=0.99, T_A=0.95, R=c(1,1), S=c(1,1), method='all', tail=0.025, bnb_priors=c(0,0), use_delta=NA, beta_iters=10^4, binomial_priors=c(1,1), binomial_cl_adj=0.2){

	# TODO: input checks

	if(!identical(R, c(1,1))){
		stop("Values of R other than 1 are not currently supported")
	}
	if(!identical(S, c(1,1))){
		stop("Values of S other than 1 are not currently supported")
	}
	mean_ratio <- 1

	# Turn delta into 0 or 2:
	use_delta <- as.integer(use_delta*2L)
	if(is.na(use_delta)){
		use_delta <- 1L
	}

	binomial_cl <- c(tail*binomial_cl_adj, 1-(tail*binomial_cl_adj))
	stopifnot(length(binomial_cl)==2 && all(!is.na(binomial_cl)) && all(binomial_cl > 0) && all(binomial_cl < 1) && binomial_cl[1] < binomial_cl[2])

	results <- RCPP_typology_analysis(as.integer(sum[1]), as.integer(sum[2]), as.integer(N[1]), as.integer(N[2]), as.double(k[1]), as.double(k[2]), as.double(cor), as.logical(paired), as.double(mean_ratio), as.double(T_I), as.double(T_A), as.double(bnb_priors), as.integer(use_delta), as.integer(beta_iters), 1L, as.double(tail), as.double(binomial_cl), as.double(binomial_priors))

	typgrp <- gsub("[[:alpha:]]","",results$Typology)
	results$Classification <- sapply(typgrp, switch, "1"="Reduced", "2"="Inconclusive", "3"="Borderline", "4"="Adequate", "Error with method")

	return(results)
}
