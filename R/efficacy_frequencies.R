#' Title
#'
#' @param rvs
#' @param N
#' @param mu
#' @param k1
#' @param k2
#' @param ta
#' @param ti
#' @param iters
#' @param approx
#' @param priors
#' @param useml
#' @param pm
#' @param plot
#'
#' @return
#' @export
#'
#' @examples
efficacy_frequencies <- function(r = NA, T_I = 0.99, T_A = 0.95, paired = TRUE, N = c(20,20), mean = 20, k = c(1, 0.7), cor = 0.5, tail=0.025, bnb_priors=c(0,0), use_delta=NA, beta_iters=10^4, useml=TRUE, binomial_priors=c(1,1), binomial_cl_adj=0.2){
	
	#TODO: parameter checks
	
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
		rvs <- unique(c(rvs, (1-tsp$ta)+c(-by/2,0,by/2), (1-tsp$ti)+c(-by/2,0,by/2)))

	}

	args <- as.list((tsp %>% select(N, mu, k1, k2, kc, ta, ti)))
	args <- c(args, list(rvs=rvs, iters=citers, approx=1, priors=c(0,0), useml=TRUE, plot=FALSE))

	res <- do.call('powersim_paired', args)
	
	
}


powersim_unpaired <- function(rvs=seq(0,0.2,by=0.001), N=20, mu=20, k1=1, k2=0.7, ta=0.9, ti=0.95, iters=10^3, approx=1, priors=c(0,0), useml=FALSE, pm=c('BNB','Levecke','MLE','WAAVP'), plot=TRUE){

	stopifnot(ta >= min(1-rvs) && ta <= max(1-rvs))
	stopifnot(ti >= min(1-rvs) && ti <= max(1-rvs))

	stopifnot(length(priors)==2)

	t <- Sys.time()
	res <- RCPP_fecrt_sim_unpaired(as.integer(iters), as.double(rvs), as.integer(N), as.integer(N), as.double(mu), as.double(k1), as.double(k2), matrix(as.double(c(ta,ti)), ncol=2, byrow=TRUE), as.double(priors), as.integer(1), as.integer(10000), as.integer(approx), as.double(0.025), as.logical(useml))
	print(difftime(Sys.time(), t))

	if(!plot){
		return(res)
	}

	library('tidyverse')
	theme_set(theme_light())

	sumres <- res %>%
		mutate(Class = gsub("[[:alpha:]]","",Classification), InfTest = Class %in% c('1','3'), NonInfTest = Class %in% c('3','4')) %>%
		select(Method, Reduction, InfTest, NonInfTest) %>%
		gather(TestType, TestResult, -Method, -Reduction) %>%
		group_by(Method, Reduction, TestType) %>%
		summarise(Probability=sum(TestResult)/n())

	pt <- ggplot(sumres[sumres$Method %in% pm, ], aes(x=Reduction, y=Probability, col=Method)) +
		geom_line() +
		geom_hline(yintercept=0.025) +
		geom_vline(aes(xintercept=thresh), data=data.frame(TestType=c('InfTest', 'NonInfTest'), thresh=1-c(ti,ta))) +
		facet_wrap( ~ TestType)

	print(pt)

	invisible(res)
}

powersim_paired <- function(rvs=seq(0,0.2,by=0.001), N=20, mu=20, k1=1, k2=0.7, kc=1.2, ta=0.9, ti=0.95, iters=10^3, approx=1, priors=c(0,0), useml=FALSE, dobson_cl=c(0.005,0.995), dobson_priors=c(1,1), pm=c('BNB','Levecke','MLE','WAAVP','Dobson'), plot=TRUE){

	stopifnot(ta >= min(1-rvs) && ta <= max(1-rvs))
	stopifnot(ti >= min(1-rvs) && ti <= max(1-rvs))

	stopifnot(length(priors)==2)
	stopifnot(length(dobson_cl)==2)
	stopifnot(length(dobson_priors)==2)

	t <- Sys.time()
	res <- RCPP_fecrt_sim_paired(as.integer(iters), as.double(rvs), as.integer(N), as.double(mu), as.double(k1), as.double(k2), as.double(kc), matrix(as.double(c(ta,ti)), ncol=2, byrow=TRUE), as.double(priors), as.integer(1), as.integer(10000), as.integer(approx), as.double(0.025), as.logical(useml), as.double(dobson_cl), as.double(dobson_priors))
	print(difftime(Sys.time(), t))

	if(!plot){
		return(res)
	}

	sumres <- res %>%
		mutate(Class = gsub("[[:alpha:]]","",Classification), InfTest = Class %in% c('1','3'), NonInfTest = Class %in% c('3','4')) %>%
		select(Method, Reduction, InfTest, NonInfTest) %>%
		gather(TestType, TestResult, -Method, -Reduction) %>%
		group_by(Method, Reduction, TestType) %>%
		summarise(Probability=sum(TestResult)/n())

	library('tidyverse')
	theme_set(theme_light())

	pt <- ggplot(sumres[sumres$Method %in% pm, ], aes(x=Reduction, y=Probability, col=Method)) +
		geom_line() +
		geom_hline(yintercept=0.025) +
		geom_vline(aes(xintercept=thresh), data=data.frame(TestType=c('InfTest', 'NonInfTest'), thresh=1-c(ti,ta))) +
		facet_wrap( ~ TestType)

	print(pt)

	invisible(res)
}
