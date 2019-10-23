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
#' ss <- efficacy_samplesize(minN = 5, maxN = 50, byN = 5, paired = TRUE,
#'     T_I = 0.99, T_A = 0.95, plot = requireNamespace("ggplot2"))
#' ss$Recommendation
#' if(requireNamespace("ggplot2") && interactive()){
#'   ss$Plot
#' }

efficacy_samplesize <- function(method = "BNB", power = 0.8, minN = 5, maxN = 50, byN = 5, plot = FALSE, ...){

	# TODO: auto-select N
	# First verify an appropriate range of N:
	# initial <- efficacy_power()

	Ns <- seq(minN, maxN, by=byN)
	stopifnot(length(method)==1)
	stopifnot(length(Ns) >= 10)
	allfreqs <- efficacy_power(N = Ns, method = method, ...)

	# Try to find a suggested sample size for each test:
	fd <- data.frame(N = seq(1,maxN*10,by=1))
	pp <- 8

	fitdata <- allfreqs %>%
		filter(.data$Test == "Inferiority", .data$Type=="Power") %>%
		mutate(N = .data$N_1)
	fit <- lm(Rate ~ poly(N, pp), data=fitdata)
	fd$Inferiority <- predict(fit, newdata=fd)

	fitdata <- allfreqs %>%
		filter(.data$Test == "NonInferiority", .data$Type=="Power") %>%
		mutate(N = .data$N_1)
	fit <- lm(Rate ~ poly(N, pp), data=fitdata)
	fd$NonInferiority <- predict(fit, newdata=fd)

	infN <- max(5, min(fd %>% filter(.data$Inferiority > power) %>% pull(.data$N)))
	noninfN <- max(5, min(fd %>% filter(.data$NonInferiority > power) %>% pull(.data$N)))

	rec <- list(Inferiority = infN, NonInferiority = noninfN, Both = max(infN, noninfN))

	rates <- allfreqs %>% mutate(N=N_1) %>% select(Method, N, Test, Type, Rate, LCI, UCI)

	pt <- NULL
	if(plot){
		if(!requireNamespace("ggplot2")){
			stop("The ggplot2 package is required to produce plots", call.=FALSE)
		}

		linedat <- fd %>% filter(N >= minN, N <= maxN) %>% gather("Test","Rate",-N) %>% mutate(Rate = pmin(Rate, 1.0)) %>% mutate(Rate = pmax(Rate, 0.0))
		pt <- ggplot(mapping = aes(x=N, y=Rate, col=Test)) +
			# geom_line(data=linedat, lty='dashed') + 
			geom_line(data=rates %>% filter(Type=="Power"), lty='dashed') + 
			geom_errorbar(aes(x=N, ymin=LCI, ymax=UCI, col=Test), rates %>% filter(Type=="Power")) +
			geom_point(data=rates %>% filter(Type=="Power")) +
			geom_hline(yintercept=power, lty='dashed') +
			geom_vline(aes(xintercept=v, col=Test), data.frame(Type=c("Power","Power"), Test=c("Inferiority","NonInferiority"), v=c(infN,noninfN)), lty='solid') +
			xlab("N") + ylab("Power") + ggtitle(paste0("Overall recommended N: ", rec$Both))
	}

	return(list(Recommendation = rec, Power = rates, Plot=pt))

}
