#' Obtain a Monte Carlo estimation of the expected frequency of typologies associated with the given parameter values
#'
#' @description this calls underlying C++ code to repeatedly simulate data from a (paired) set of negative binomial distributions with given parameters, analyse the resultant data with each of four/five statistical methods, and return the expected frequencies of each typology encountered using each method
#'
#' @param N the sample sizes to examine for type 1 and 2 error rates, either as an integer vector where it is presumed that N_1 = N_2, or a 2-column matrix of paired N_1 and N_2 values to examine at each combination given by row
#' @param method the statistical method(s) for which to return type 1 and type 2 error rates
#' @param paired logical flag for a paired or unpaired study design assumption on which to simulate (and analyse) the data
#' @param T_I the threshold for inferioirty (target efficacy of the intervention, and also the value of efficacy tested for the non-inferiority test)
#' @param T_A the threshold for non-inferioirty (target efficacy of the intervention minus a non-inferioirty margin delta, and also the value of efficacy tested for the inferiority test)
#' @param ... other parameters to be passed to the underlying \code{\link{efficacy_frequencies}} and \code{\link[parallel]{mclapply}} functions
#'
#' @return A data frame with columns reflecting the method, N_1, N_2, hypothesis test, error type, and error rate.
#' @export
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom tidyr gather
#' @importFrom parallel mclapply
#' @seealso \code{\link{bayescount}}
#'
#' @examples
#' rates <- efficacy_power(N = seq(5,50,by=5), paired = TRUE, T_I = 0.99, T_A = 0.95)
#' if(requireNamespace("ggplot2")){
#'   ggplot(rates, aes(x=N_1, y=Rate, ymin=LCI, ymax=UCI, col=Test)) +
#'     geom_errorbar() + geom_line(lty='dashed') + geom_point() +
#'     geom_hline(aes(yintercept=h), data.frame(Type=c("Power","Type1Error"), h=c(0.8,0.025)), lty='dashed') +
#'     facet_wrap(~Type, scales='free_y') + xlab("N")
#' }

efficacy_power <- function(N = seq(5,50,by=5), method = "BNB", paired = TRUE, T_I = 0.99, T_A = 0.95, ...){

	#TODO: parameter checks
	#TODO: find and use mclapply arguments from ...

	checksingleprob(T_I)
	checksingleprob(T_A)
	if(! T_I >= T_A){
		stop("The value supplied for T_I must be greater than or equal to T_A", call.=FALSE)
	}

	if(!is.matrix(N)){
		stopifnot(is.numeric(N))
		N <- cbind(N_1=N,N_2=N)
	}
	stopifnot(is.matrix(N))
	
	possmethods <- c("BNB", "Gamma", "WAAVP", "Asymptotic", "Binomial")
	if('all' %in% method){
		method <- possmethods
	}
	stopifnot(all(method %in% possmethods))

		
	# Get output from underlying function:
	outputs <- mclapply(1:nrow(N), function(i) efficacy_frequencies(r = c(T_I,T_A), paired = paired, T_I = T_I, T_A = T_A, N = N[i,], ...))
	
	# Harvest constant parameter value attribute:
	parameters <- attr(outputs[[1]],'parameters')
	parameters$N <- NULL
	iters <- parameters$iterations
	
	# Convert to individual tests and error rates:
	allfreqs <- do.call('rbind', lapply(1:nrow(N), function(i) cbind(data.frame(N_1=as.integer(N[i,1]), N_2=as.integer(N[i,2]), outputs[[i]])))) %>%
		filter(.data$Method %in% method) %>%
		group_by(.data$N_1, .data$N_2, .data$Efficacy, .data$Method, .data$Classification) %>%
		summarise(Frequency = sum(.data$Frequency), Proportion = sum(.data$Proportion)) %>%
		ungroup() %>%
		mutate(Inferiority = .data$Classification %in% c("Reduced","Borderline"), NonInferiority = .data$Classification %in% c("Adequate","Borderline")) %>%
		gather("Test","Result", -.data$N_1, -.data$N_2, -.data$Efficacy, -.data$Method, -.data$Classification, -.data$Frequency, -.data$Proportion) %>%
		group_by(.data$N_1, .data$N_2, .data$Efficacy, .data$Method, .data$Test, .data$Result) %>%
		summarise(Frequency = sum(.data$Frequency), Proportion = sum(.data$Proportion)) %>%
		ungroup() %>%
		filter(Result == TRUE) %>%
		mutate(Type = ifelse((Efficacy == T_I & Test == "NonInferiority") | (Efficacy == T_A & Test == "Inferiority"), "Power", "Type1Error")) %>%
		mutate(LCI = qbeta(0.025, Frequency+1, (iters-Frequency)+1), UCI = qbeta(0.975, Frequency+1, (iters-Frequency)+1)) %>%
		select(Method, N_1, N_2, Test, Type, Rate=Proportion, LCI, UCI)
	
	# Add simulation parameters:
	attr(allfreqs, "parameters") <- parameters
	
	return(allfreqs)
}
