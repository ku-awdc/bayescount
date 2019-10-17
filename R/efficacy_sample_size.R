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
efficacy_frequencies <- function(){
	
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
