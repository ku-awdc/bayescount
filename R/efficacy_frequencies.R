#' Obtain a Monte Carlo estimation of the expected frequency of typologies associated with the given parameter values
#'
#' @description this calls underlying C++ code to repeatedly simulate data from a (paired) set of negative binomial distributions with given parameters, analyse the resultant data with each of four/five statistical methods, and return the expected frequencies of each typology encountered using each method
#'
#' @param r a vector representing a set of true population reductions to simulate.  If NA then values of r over an appropriate range of interest (based on T_I and T_A) are automatically generated.
#' @param paired logical flag for a paired or unpaired study design assumption on which to simulate (and analyse) the data
#' @param T_I the threshold for inferioirty (target efficacy of the intervention)
#' @param T_A the threshold for non-inferioirty (target efficacy of the intervention minus a non-inferioirty margin delta)
#' @param N a length-2 integer giving the number of control and treatment individuals for an unpaired study, or a length-1 integer giving the number of individuals for a paired study
#' @param R a length-2 integer giving the number of replicates used for control/pre-treatment and treatment/post-treatment individuals
#' @param S a length-2 numeric variable giving the counting sensitivity used for control/pre-treatment and treatment/post-treatment individuals
#' @param mean the control/pre-treatment mean count
#' @param k a length-2 strictly positive continuous number giving the over-dispersion of the control/pre-treatment and treatment/post-treatment data
#' @param cor a length-1 number between 0 and 1 giving the correlation between pre- and post-treatment data (ignored for an unpaired study)
#' @param iterations the number of iterations to use for the Monte Carlo approximation
#' @param alpha the significance level to use for the classification, where 0.025 corresponds to 95\% CI and 0.05 corresponds to 90\% CI
#' @param bnb_priors priors to use for the BNB method
#' @param use_delta logical flag to use the delta method approximation for the BNB method (NA means to use it unless it fails)
#' @param beta_iters number of iterations to use for the Monte Carlo approximation of the beta distribution transformation for the BNB method (when use_delta==FALSE)
#' @param use_ml a logical flag controlling if maximum likelihood is used (rather than less accurate but faster methods) to estimate the over-dispersion parameter k
#' @param binomial_priors the priors to use for the Binomial method
#' @param binomial_cl_adj the adjustment in confidence level to use for the Binomial method relative to the others (the default corresponds to e.g. 99\% rather than 95\% CI) - where 0 > adj >= 1
#'
#' @return A data frame with columns reflecting the method, population reduction r, typology, classification, frequency, and proportion within method & r.  The parameter values used for the simulation are given as an attribute list.
#' @export
#' @import dplyr
#' @importFrom rlang .data
#' @seealso \code{\link{bayescount}}
#'
#' @examples
#' # Monte Carlo estimates of typology frequencies for all five methods
#' # at the critical values corresponding to T_I and T_A:
#' ( critical_values <- efficacy_frequencies(r=c(0.95, 0.99), paired = TRUE, T_I = 0.99, T_A = 0.95) )
#' # The same for an unpaired analysis (excluding the Binimial method):
#' ( critical_values <- efficacy_frequencies(r=c(0.95, 0.99), paired = FALSE, T_I = 0.99, T_A = 0.95) )
efficacy_frequencies <- function(r = NA, paired = TRUE, T_I = 0.99, T_A = 0.95, N = c(20,20), R = c(1,1), S = c(1,1), mean = 20, k = c(1, 0.7), cor = 0.5, iterations=10^3, alpha=0.025, bnb_priors=c(0,0), use_delta=NA, beta_iters=10^4, use_ml=TRUE, binomial_priors=c(1,1), binomial_cl_adj=0.2){

	#TODO: parameter checks
	#TODO: remove clash for r and R arguments

	checksingleprob(T_I)
	checksingleprob(T_A)
	if(! T_I >= T_A){
		stop("The value supplied for T_I must be greater than or equal to T_A", call.=FALSE)
	}

	if(!identical(R, c(1,1))){
		stop("Values of R other than 1 are not currently supported")
	}
	if(!identical(S, c(1,1))){
		stop("Values of S other than 1 are not currently supported")
	}

	if(identical(r, NA)){
		delta <- T_I - T_A
		# Get lower and upper as function of ti and ta:
		upper <- 1 - max(0, T_A-delta)
		lower <- 1 - min(1, T_I+delta)

		# Make by a sensible number to give around 25 values depending on upper-lower:
		by <- seq(0,upper-lower,length.out=25)[2]
		choices <- c(0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05)
		by <- choices[which.min(abs(by-choices))]

		rvs <- seq(lower,upper,by=by)
		if(0 %in% rvs)
			rvs <- unique(c(0, choices[1:6], rvs))
		rvs <- unique(c(rvs, (1-T_A)+c(-by/2,0,by/2), (1-T_I)+c(-by/2,0,by/2)))

	}else{
		sapply(r, checksingleprob)
		rvs <- 1-r
	}

	# Turn delta into 0 or 2:
	use_delta <- as.integer(use_delta*2L)
	if(is.na(use_delta)){
		use_delta <- 1L
	}

	binomial_cl <- c(alpha*binomial_cl_adj, 1-(alpha*binomial_cl_adj))
	stopifnot(length(binomial_cl)==2 && all(!is.na(binomial_cl)) && all(binomial_cl > 0) && all(binomial_cl < 1) && binomial_cl[1] < binomial_cl[2])

	# Controls if the large-sample approximation is used (1 means if necessary):
	approx <- 1

	if(paired){
		# Calculate ks based on correlation:
		kc <- adjust_k(k[1], k[2], cor)$correlated_k

		if(length(N)==1) N <- rep(N,2)

		results <- RCPP_efficacy_frequencies_paired(as.integer(iterations), as.double(rvs), as.integer(N[1]), as.double(mean), as.double(k[1]), as.double(k[2]), as.double(kc), matrix(as.double(c(T_A,T_I)), ncol=2, byrow=TRUE), as.double(bnb_priors), as.integer(use_delta), as.integer(beta_iters), as.integer(approx), as.double(alpha), as.logical(use_ml), as.double(binomial_cl), as.double(binomial_priors))

	}else{
		results <- RCPP_efficacy_frequencies_unpaired(as.integer(iterations), as.double(rvs), as.integer(N[1]), as.integer(N[2]), as.double(mean), as.double(k[1]), as.double(k[2]), matrix(as.double(c(T_A,T_I)), ncol=2, byrow=TRUE), as.double(bnb_priors), as.integer(use_delta), as.integer(beta_iters), as.integer(approx), as.double(alpha), as.logical(use_ml))
	}

	sumres <- results %>%
		mutate(Efficacy = 1-.data$Reduction, Classification="") %>%
		group_by(.data$Method, .data$Efficacy, .data$Typology, .data$Classification) %>%
		summarise(Frequency = n(), Proportion=.data$Frequency/iterations) %>%
		ungroup()

	# Add any non-observed typologies explicitly:
	sumres <- full_join(sumres, expand.grid(Typology=c('1ab','1c','2a','2b','2c','3','4a','4bc'), Method=unique(sumres$Method), Efficacy=unique(sumres$Efficacy), stringsAsFactors=FALSE), by=c("Typology","Method","Efficacy"))
	sumres$Frequency[is.na(sumres$Frequency)] <- 0
	sumres$Proportion[is.na(sumres$Proportion)] <- 0

	typgrp <- gsub("[[:alpha:]]","",sumres$Typology)
	sumres$Classification <- factor(sapply(typgrp, switch, "1"="Reduced", "2"="Inconclusive", "3"="Borderline", "4"="Adequate", "Method_Failure"), levels=c("Reduced","Inconclusive","Borderline","Adequate","Method_Failure"))

	# Remove unwanted attributes:
	attr(sumres, "out.attrs") <- NULL

	# Add simulation parameters:
	attr(sumres, "parameters") <- list(paired = paired, T_I = T_I, T_A = T_A, N = N, R = R, S = S, mean = mean, k = k, cor = cor, iterations=iterations, alpha=alpha, bnb_priors=bnb_priors, use_delta=use_delta, beta_iters=beta_iters, use_ml=use_ml, binomial_priors=binomial_priors, binomial_cl_adj=binomial_cl_adj)

	return(sumres)
}
