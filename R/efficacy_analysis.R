#' efficacy_analysis
#'
#' @description
#' This function analyses a given count reduction dataset (e.g. egg reduction rate or faecal egg count reduction) using five different methods of analysis
#'
#' @param data_1 either the pre-treatment or control data - should either be an integer vector, or an integer matrix with replicates from the same individual in columns
#' @param data_2 either the post-treatment or control data - should either be an integer vector, or an integer matrix with replicates from the same individual in columns
#' @param paired logical flag for a paired or unpaired study design - if paired then length or nrow of data_1 must equal that of data_2
#' @param T_I the threshold for inferioirty (target efficacy of the intervention)
#' @param T_A the threshold for non-inferioirty (target efficacy of the intervention minus a non-inferioirty margin delta)
#' @param alpha the significance level to use for the classification, where 0.025 corresponds to 95\% CI and 0.05 corresponds to 90\% CI
#' @param S a length-2 numeric variable giving the counting sensitivity used for control/pre-treatment and treatment/post-treatment individuals
#' @param bnb_priors priors to use for the BNB method
#' @param use_delta logical flag to use the delta method approximation for the BNB method (NA means to use it unless it fails)
#' @param beta_iters number of iterations to use for the Monte Carlo approximation of the beta distribution transformation for the BNB method (when use_delta==FALSE)
#' @param use_ml a logical flag controlling if maximum likelihood is used (rather than less accurate but faster methods) to estimate the over-dispersion parameter k
#' @param binomial_priors the priors to use for the Binomial method
#' @param binomial_cl_adj the adjustment in confidence level to use for the Binomial method relative to the others (the default corresponds to e.g. 99\% rather than 95\% CI) - where 0 > adj >= 1
#'
#' @return A data frame with columns Method, LCI, UCI, pI, pA, Typology, and Classification
#' @export
#' @seealso \code{\link{bayescount}}, \code{\link{efficacy_frequencies}}
#'
#' @examples
#' (data1 <- rnbinom(10, mu=20, size=1))
#' (data2 <- rnbinom(10, mu=2, size=0.75))[1:5]
#' efficacy_analysis(data1, data2, paired=FALSE, T_I=0.95, T_A=0.9)
efficacy_analysis <- function(data_1, data_2, paired, T_I=0.99, T_A=0.95, S=c(1,1), alpha=0.05, bnb_priors=c(0,0), use_delta=NA, beta_iters=10^4, use_ml=TRUE, binomial_priors=c(1,1), binomial_cl_adj=0.2, known_ks=c(1.0, 1.0), k_ratio=1.0) {

	# TODO: input checks

	if(is.matrix(data_1)){
		if(!is.matrix(data_2)){
			stop("The type (vector vs matrix) of data_1 must match that of data_2")
		}
		if(ncol(data_1) != ncol(data_2)){
			stop("It is currently a requirement that the number of columns of data_1 and data_2 must match")
		}

		# TODO: work out which way mean_rato should go and implement in CI functions in C++:
		mean_ratio <- 1
		data_1 <- apply(data_1, 1, sum)
		data_2 <- apply(data_2, 1, sum)

	}else{
		mean_ratio <- 1
	}
	if(!identical(S, c(1,1))){
		stop("Values of S other than 1 are not currently supported")
	}

	if(paired && length(data_1)!=length(data_2)){
		stop("The number of observations in pre- and post-treatment data must match for paired data")
	}

	# Turn delta into 0 or 2:
	use_delta <- as.integer(use_delta*2L)
	if(is.na(use_delta)){
		use_delta <- 1L
	}

	binomial_cl <- c(alpha*binomial_cl_adj, 1-(alpha*binomial_cl_adj))
	stopifnot(length(binomial_cl)==2 && all(!is.na(binomial_cl)) && all(binomial_cl > 0) && all(binomial_cl < 1) && binomial_cl[1] < binomial_cl[2])

	stopifnot(length(known_ks)==2, all(!is.na(known_ks)), all(known_ks>0.0))
	stopifnot(length(k_ratio)==1, all(!is.na(k_ratio)), all(k_ratio>0.0))

	results <- RCPP_efficacy_analysis(as.integer(data_1), as.integer(data_2), as.logical(paired), as.double(mean_ratio), as.double(T_I), as.double(T_A), as.double(bnb_priors), as.integer(use_delta), as.integer(beta_iters), as.logical(use_ml), 1L, as.double(alpha), as.double(binomial_cl), as.double(binomial_priors), as.double(known_ks), as.double(k_ratio))

	typgrp <- gsub("[[:alpha:]]","",results$Typology)
	results$Classification <- sapply(typgrp, switch, "1"="Reduced", "2"="Inconclusive", "3"="Borderline", "4"="Adequate", "Method_Failure")

	results$Method <- factor(results$Method, levels=c('WAAVP','Gamma','Binomial','Asymptotic','BNB_KnownKs','BNB_FixK2','BNB'))
	results$Typology <- factor(results$Typology, levels=c("1ab", "1c", "2a", "2b","2c", "3", "4a","4bc","100%red", "<0%red", "error"))
	results$Classification <- factor(results$Classification, levels=c("Reduced", "Inconclusive", "Borderline", "Adequate","Method_Failure"), labels=c("Resistant", "Inconclusive", "(Low) Resistant", "Susceptible","Method_Failure"))

	return(results)
}

