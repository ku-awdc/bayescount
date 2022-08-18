checksingleprob <- function(x, na.ok=FALSE){
	if(length(x)!=1){
		stop(paste(substitute(x), 'must be of length 1'), call.=FALSE)
	}
	if(!is.null(dim(x))){
		stop(paste(substitute(x), 'must not have dimensions'), call.=FALSE)
	}
	if(is.na(x)){
		if(!na.ok){
			stop(paste(substitute(x), 'must be non-missing'), call.=FALSE)
		}
	}else{
		if(!is.numeric(x)){
			stop(paste(substitute(x), 'must be numeric'), call.=FALSE)
		}
		if(! x>=0){
			stop(paste(substitute(x), 'must be >=0'), call.=FALSE)
		}
		if(! x<=1){
			stop(paste(substitute(x), 'must be <=1'), call.=FALSE)
		}
	}
}

checksinglelogical <- function(x, na.ok=FALSE){
	if(length(x)!=1){
		stop(paste(substitute(x), 'must be of length 1'), call.=FALSE)
	}
	if(!is.null(dim(x))){
		stop(paste(substitute(x), 'must not have dimensions'), call.=FALSE)
	}
	if(is.na(x)){
		if(!na.ok){
			stop(paste(substitute(x), 'must be non-missing'), call.=FALSE)
		}
	}else{
		if(!is.logical(x)){
			stop(paste(substitute(x), 'must be logical'), call.=FALSE)
		}
		if(! x>=0){
			stop(paste(substitute(x), 'must be >=0'), call.=FALSE)
		}
		if(! x<=1){
			stop(paste(substitute(x), 'must be <=1'), call.=FALSE)
		}
	}
}

checksingleposint <- function(x, na.ok=FALSE){
	if(length(x)!=1){
		stop(paste(substitute(x), 'must be of length 1'), call.=FALSE)
	}
	if(!is.null(dim(x))){
		stop(paste(substitute(x), 'must not have dimensions'), call.=FALSE)
	}
	if(is.na(x)){
		if(!na.ok){
			stop(paste(substitute(x), 'must be non-missing'), call.=FALSE)
		}
	}else{
		if(!is.numeric(x)){
			stop(paste(substitute(x), 'must be numeric'), call.=FALSE)
		}
		if(! x>0){
			stop(paste(substitute(x), 'must be >0'), call.=FALSE)
		}
		if(x%%1 !=0){
			stop(paste(substitute(x), 'must be an integer'), call.=FALSE)
		}
	}
}

checksingleposdouble <- function(x, na.ok=FALSE){
	if(length(x)!=1){
		stop(paste(substitute(x), 'must be of length 1'), call.=FALSE)
	}
	if(!is.null(dim(x))){
		stop(paste(substitute(x), 'must not have dimensions'), call.=FALSE)
	}
	if(is.na(x)){
		if(!na.ok){
			stop(paste(substitute(x), 'must be non-missing'), call.=FALSE)
		}
	}else{
		if(!is.numeric(x)){
			stop(paste(substitute(x), 'must be numeric'), call.=FALSE)
		}
		if(! x>0){
			stop(paste(substitute(x), 'must be >0'), call.=FALSE)
		}
	}
}

checklen2posdouble <- function(x, na.ok=FALSE){
	if(length(x)!=2){
		stop(paste(substitute(x), 'must be of length 2'), call.=FALSE)
	}
	if(!is.null(dim(x))){
		stop(paste(substitute(x), 'must not have dimensions'), call.=FALSE)
	}
	if(any(is.na(x))){
		if(!na.ok){
			stop(paste(substitute(x), 'must all be non-missing'), call.=FALSE)
		}
	}else{
		if(!is.numeric(x)){
			stop(paste(substitute(x), 'must be numeric'), call.=FALSE)
		}
		if(! all(x>0)){
			stop(paste(substitute(x), 'must all be >0'), call.=FALSE)
		}
	}
}

checkintlimit <- function(x){
	if(x > .Machine$integer.max || x < -.Machine$integer.max){
		stop(paste(substitute(x), 'is beyond the integer limits of your machine so cannot be passed to the underlying C++ code'), call.=FALSE)
	}
	return(as.integer(x))
}



# TODO: clean up
powersim_unpaired <- function(rvs=seq(0,0.2,by=0.001), N=20, mu=20, k1=1, k2=0.7, ta=0.9, ti=0.95, iters=10^3, approx=1, priors=c(0,0), useml=FALSE, pm=c('BNB','Levecke','MLE','WAAVP'), plot=TRUE){

	stopifnot(ta >= min(1-rvs) && ta <= max(1-rvs))
	stopifnot(ti >= min(1-rvs) && ti <= max(1-rvs))

	stopifnot(length(priors)==2)

	t <- Sys.time()
	res <- fecrt_sim_unpaired(as.integer(iters), as.double(rvs), as.integer(N), as.integer(N), as.double(mu), as.double(k1), as.double(k2), matrix(as.double(c(ta,ti)), ncol=2, byrow=TRUE), as.double(priors), as.integer(1), as.integer(10000), as.integer(approx), as.double(0.025), as.logical(useml))
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

#' @export
powersim_paired_LP <- function(rvs=seq(0,0.2,by=0.001), N=20, mu=20, cv1=1, cv2=0.7, ncor=0.2, ta=0.9, ti=0.95, iters=10^3, approx=1, priors=c(0,0), useml=FALSE, uselp=TRUE, dobson_cl=c(0.005,0.995), dobson_priors=c(1,1), pm=c('BNB','Levecke','MLE','WAAVP','Dobson'), plot=TRUE){

	stopifnot(ta >= min(1-rvs) && ta <= max(1-rvs))
	stopifnot(ti >= min(1-rvs) && ti <= max(1-rvs))

	stopifnot(length(priors)==2)
	stopifnot(length(dobson_cl)==2)
	stopifnot(length(dobson_priors)==2)

	# Note: cvs are unadjusted, ncor is correlation in bivariate normal

	t <- Sys.time()
	res <- fecrt_sim_paired_LP(as.integer(iters), as.double(rvs), as.integer(N), as.double(mu), as.double(cv1), as.double(cv2), as.double(ncor), matrix(as.double(c(ta,ti)), ncol=2, byrow=TRUE), as.double(priors), as.integer(1), as.integer(10000), as.integer(approx), as.double(0.025), as.logical(useml), as.double(dobson_cl), as.double(dobson_priors), as.logical(uselp))
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


findtheta <- function(data, min=0.001, max=20, tol=0.001){
	# Equivalent to:
	# print(optimize(function(theta) sum(dnbinom(data,theta,mu=mean(data),log=TRUE)), interval=c(min,max), maximum=TRUE, tol=tol)$maximum)
	find_theta(as.integer(data), as.double(mean(data)), as.double(min), as.double(max), as.double(tol))
}



summarise_fecr <- function(data_1, data_2, label=""){

	stopifnot(sum(is.na(data_1))==0)
	stopifnot(sum(is.na(data_2))==0)
	stopifnot(length(data_1)==length(data_2))
	N <- length(data_1)

	stopifnot(all(label==label[1]))

	mean_1 <- mean(data_1)
	var_1 <- var(data_1)
	mean_2 <- mean(data_2)
	var_2 <- var(data_2)
	cov_12 <- cov(data_1, data_2)

	# Efficacy:
	eff <- 1 - mean(data_2)/mean(data_1)

	# Remove Poisson variation and calculate adjusted correlation:
	var_eff_1 <- var_1 - mean_1
	var_eff_2 <- var_2 - mean_2
	correlation <- cov_12 / (sqrt(var_eff_1) *sqrt(var_eff_2))

	# Get estimated unpaired=overall and paired ks:
	ko <- estimate_k(as.double(mean_1), as.double(var_1), as.double(mean_2), as.double(var_2), as.double(0.0), as.logical(FALSE))
	kp <- estimate_k(as.double(mean_1), as.double(var_1), as.double(mean_2), as.double(var_2), as.double(cov_12), as.logical(TRUE))

	# Estimate kc (correlated k) from ko and correlation:
	kc <- sqrt(ko[1]*ko[2]) / correlation#cor(data_1,data_2)
	# Estimate ku (uncorrelated k) from kc and ko:
	ku <- ((kc+1)*ko)/(kc-ko)

	df <- data.frame(label=label[1], N=N, mean_1=mean_1, mean_2=mean_2, efficacy=eff, var_1=var_1, var_eff_1=var_eff_1, var_2=var_2, var_eff_2=var_eff_2, cov=cov_12, cor=cor(data_1, data_2), cor_eff=correlation, ko_1=ko[1], ko_2=ko[2], kc=kc, ku_1=ku[1], ku_2=ku[2], kp_1=kp[1], kp_2=kp[2])

	if(identical(label,"")){
		df$label <- NULL
	}

	return(df)
}

analyse_fecr <- function(data_1, data_2, paired=FALSE, ta=0.9, ti=0.95, dobson_cl=c(0.005,0.995), dobson_priors=c(1,1)){


	stopifnot(sum(is.na(data_1))==0)
	stopifnot(sum(is.na(data_2))==0)
	if(paired){
		stopifnot(length(data_1)==length(data_2))
	}

	mean_1 <- mean(data_1)
	var_1 <- var(data_1)
	mean_2 <- mean(data_2)
	var_2 <- var(data_2)

	ks_p <- estimate_k_ml(as.double(data_1), as.double(mean_1), as.double(var_1), as.double(data_2), as.double(mean_2), as.double(var_2), as.double(cov(data_1, data_2)), as.logical(paired))
	ks_u <- estimate_k_ml(as.double(data_1), as.double(mean_1), as.double(var_1), as.double(data_2), as.double(mean_2), as.double(var_2), as.double(0.0), as.logical(paired))


	# Efficacy:
	eff <- 1 - mean(data_2)/mean(data_1)

	# Get estimated unpaired and paired ks:

	# Rcpp::NumericMatrix methodcomp(long long sum1, int N1, double K1, double mu1, double var1, long long sum2, int N2, double K2, double mu2, double var2, double cov12, double mean_ratio, double H0_1, double H0_2, Rcpp::NumericVector conjugate_priors, int delta, int beta_iters, int approx, double tail)

	rv <- methodcomp(as.integer(sum(data_1)), as.integer(length(data_1)), as.double(ks_p[1]), as.double(ks_u[1]), as.double(mean(data_1)), as.double(var(data_1)), as.integer(sum(data_2)), as.integer(length(data_2)), as.double(ks_p[2]), as.double(ks_u[2]), as.double(mean(data_2)), as.double(var(data_2)), as.double(cov(data_1,data_2)), as.double(1.0), as.double(ta), as.double(ti), as.double(c(1.0,1.0)), as.integer(1), as.integer(10000), as.integer(1), as.double(0.025), as.double(dobson_cl), as.double(dobson_priors))

	if(!paired){
		return(rv[seq(2,8,by=2),])
	}else{
		return(rv[seq(1,9,by=2),])
	}

}




if(FALSE){


	## Testing pbnb:

	for(i in 0:10){
		noninc <- ifelse(i==0, 0, pbnb(i-1, 1, 20, 20, TRUE, TRUE))
		inc <- pbnb(i, 1, 20, 20, TRUE, TRUE)
		print(inc - noninc)
	}
	for(i in 0:10){
		noninc <- pbnb(i, 1, 20, 20, TRUE, FALSE)
		inc <- pbnb(i, 1, 20, 20, TRUE, TRUE)
		print(inc - noninc)
	}

	for(i in 0:10){
		inc <- 1 - ifelse(i==0, 0, pbnb(i-1, 1, 20, 20, TRUE, TRUE))
		noninc <- 1 - pbnb(i, 1, 20, 20, TRUE, TRUE)
		print(inc - noninc)
	}
	for(i in 0:10){
		noninc <- pbnb(i, 1, 20, 20, FALSE, FALSE)
		inc <- pbnb(i, 1, 20, 20, FALSE, TRUE)
		print(inc - noninc)
	}

}


# This is replaced by the c function now:
estimate_ks_old <- function(pre_data, post_data, true_k=NA, k_change=NA){
	if(!is.na(true_k)){
		preK <- true_k
		if(is.na(k_change)){
			postK <- preK
		}else{
			postK <- preK * k_change
		}
	}else{
		preK <- mean(pre_data)^2 / (var(pre_data) - mean(pre_data))

		# If preK can't be calculated then set to 10 with a warning:
		if(var(pre_data) < mean(pre_data)){
			warning('The mean of the untreated data is higher than its variance: using an estimate of k=10', call.=FALSE)
			preK <- 10
		}
		if(preK > 10){
			warning('The variance of the untreated data is very close to its variance: using an estimate of k=10', call.=FALSE)
			preK <- 10
		}
		# This is consistent with the C++ code in fecrt_power

		postK <- mean(post_data)^2 / (var(post_data) - mean(post_data))

		# If k_change is set always use it and ignore above calculation:
		if(!is.na(k_change)){
			postK <- preK * k_change

		# Otherwise if postK can't be calculated:
		}else if(sum(post_data)==0 || var(post_data) < mean(post_data) || postK > 10){
			postK <- preK
		}
	}
	stopifnot(preK > 0 & postK > 0)
	return(list(K_untx=preK, K_tx=postK))
}

# Not using as impossible to correctly estimate post k from pre k if replicates changes anyway - but keep for prosperity:
estimate_ks_old_old <- function(pre_data, post_data, true_k=NA, k_change=NA, rep_pre=1, rep_post=1, pool_pre=1, pool_post=1){
	if(!is.na(true_k)){
		preK <- true_k
		if(is.na(k_change)){
			postK <- (true_k/(rep_pre * pool_pre)) * rep_post * pool_post
		}else{
			postK <- (true_k/(rep_pre * pool_pre)) * rep_post * pool_post * k_change
		}
	}else{
		preK <- mean(pre_data)^2 / (var(pre_data) - mean(pre_data))

		# If preK can't be calculated then set to 10 with a warning:
		if(var(pre_data) < mean(pre_data)){
			warning('The mean of the pre-treatment/control data is higher than its variance: using an estimate of k=10', call.=FALSE)
			preK <- 10
		}
		if(preK > 10){
			warning('The variance of the pre-treatment/control data is very close to its variance: using an estimate of k=10', call.=FALSE)
			preK <- 10
		}
		# This is consistent with the C++ code in fecrt_power

		postK <- mean(post_data)^2 / (var(post_data) - mean(post_data))

		# If k_change is set always use it and ignore above calculation:
		if(!is.na(k_change)){
			postK <- (preK/(rep_pre * pool_pre)) * rep_post * pool_post * k_change

		# Otherwise if postK can't be calculated:
		}else if(sum(post_data)==0 || var(post_data) < mean(post_data) || postK > 10){
			postK <- (preK/(rep_pre * pool_pre)) * rep_post * pool_post
		}
	}
	stopifnot(preK > 0 & postK > 0)
	return(list(preK=preK, postK=postK))
}

fecrt_analyses <- function(untreated_data, treated_data, edt_ratio=1, replicates_ratio=1, H0_1=0.9, H0_2=0.95, tail=0.025, k_change=NA, true_k=NA, delta_method=1, beta_iters=10^6, conjugate_priors=c(1,1), dobson_ci_lower=0.005, dobson_ci_upper=0.995, dobson_priors=c(1,1), positive_covariance=FALSE){

	##### Argument checks

	warning('k_change is different here than in the C++ code - it should be changed to a default of 1, only to be used when post k is based on pre k')

	# The untreated_data and treated_data must be integer (egg counts), and given as vectors reflecting the sum of eggs counted
	stopifnot(is.numeric(untreated_data))
	stopifnot(is.null(dim(untreated_data)))
	stopifnot(is.numeric(treated_data))
	stopifnot(is.null(dim(treated_data)))

	# Data might be paired only if the lengths are the same and non are removed as missing:
	paired <- length(untreated_data) == length(treated_data)
	if(any(is.na(untreated_data))){
		warning('Removing missing values in the untreated data', call.=FALSE)
		untreated_data <- as.numeric(na.omit(untreated_data))
		paired <- FALSE
	}
	if(any(is.na(treated_data))){
		warning('Removing missing values in the treated data', call.=FALSE)
		treated_data <- as.numeric(na.omit(treated_data))
		paired <- FALSE
	}
	N1 <- length(untreated_data)
	N2 <- length(treated_data)

	if(N1 < 2 || N2 < 2){
		stop('Either the (non-missing) treated or untreated data is less than length 2', call.=FALSE)
	}
	if(any(untreated_data < 0)){
		stop('Negative values are not permitted in the (untreated) data', call.=FALSE)
	}
	if(any(treated_data < 0)){
		stop('Negative values are not permitted in the (treated) data', call.=FALSE)
	}
	if(!sum(untreated_data) > 0){
		stop('The untreated data must sum to >=1', call.=FALSE)
	}

	# Check integers:
	if(any((untreated_data %% 1) > 0)){
		stop('Non-integer values are not permitted in the (untreated) data', call.=FALSE)
	}
	if(any((treated_data %% 1) > 0)){
		stop('Non-integer values are not permitted in the (treated) data', call.=FALSE)
	}

	# Check other arguments:

	# edt_ratio is treated threshold / untreated threshold - so expected to be <=1
	checksingleposdouble(edt_ratio)
	if(edt_ratio > 1){
		cat('Note:  an increasing edt_ratio implies that a less sensitive diagnostic technique was used for the treated data!\n')
	}
	# replicates_ratio is treated replicates / untreated replicates - so expected to be >=1
	checksingleposdouble(replicates_ratio)
	if(replicates_ratio < 1){
		cat('Note:  a decreasing replicates_ratio implies that fewer replicate samples were taken for the treated data!\n')
	}
	checksingleprob(H0_1)
	checksingleprob(H0_2)
	checksingleprob(tail)
	checklen2posdouble(conjugate_priors)

	checksingleposdouble(k_change, na.ok=TRUE)
	checksingleposdouble(true_k, na.ok=TRUE)

	# Note: delta_method can take 3 values:  0=never, 1=unless_fails, 2=always
	checksingleposint(delta_method)
	stopifnot(delta_method %in% 0:2)

	# Note: bnb_approximation can take 3 values:  0=never, 1=if_necessary, 2=always
	bnb_approximation <- 1   # Currently this is actually ignored and the bnb function is run 3 times for comparison/checking
	checksingleposint(bnb_approximation)
	stopifnot(bnb_approximation %in% 0:2)

	checksingleposint(beta_iters)
	stopifnot(beta_iters >= 1000)
	stopifnot(beta_iters <= 10^8)
	# Max 10^8 - 10^6 is pretty instant and gives the same results

	checksingleprob(dobson_ci_lower)
	checksingleprob(dobson_ci_upper)
	checklen2posdouble(dobson_priors)
	checksinglelogical(positive_covariance)


	##### Further argument processing:

	# NB edt is a threshold, i.e opposite of replicates!
	mean_ratio <- replicates_ratio/edt_ratio
	checksingleposdouble(mean_ratio)
	# Need to adjust treated data for WAAVP etc methods according to the change in mean from replicates and edt:
	treated_corrected <- treated_data / mean_ratio

	sum_untx <- sum(untreated_data)
	sum_tx <- sum(treated_data)
	ks <- estimate_ks(untreated_data, treated_data, true_k=true_k, k_change=k_change)

	mu_untx <- mean(untreated_data)
	mu_tx <- mean(treated_corrected)
	var_untx <- var(untreated_data)
	if(sum(treated_corrected)>0){
		var_tx <- var(treated_corrected)
		covar <- 0
		if(paired){
			covar <- cov(untreated_data, treated_corrected)
			if(positive_covariance && covar < 0){
				covar <- 0
			}
		}
	}else{
		var_tx <- 0
		covar <- 0
	}


	##### Calculate outputs:

	warning('Conjugate Beta could also be made possibly paired by adjusting K')

	methods <- c('BNB unpaired (no approximation)', 'BNB unpaired (with approximation)', 'BNB unpaired', 'WAAVP unpaired', 'Levecke unpaired', 'MLE unpaired', 'BNB paired (no approximation)', 'BNB paired (with approximation)', 'BNB paired', 'WAAVP paired', 'Levecke paired', 'MLE paired', 'Conjugate Beta', 'Dobson')
	rownums <- 1:length(methods)
	names(rownums) <- methods
	types <- c('Paired','Unpaired')[c(2,2,2,2,2,2,1,1,1,1,1,1,1,2)]
	alloutputs <- data.frame(Method=methods, Type=types, LowerCI=as.numeric(NA), UpperCI=as.numeric(NA), p1=as.numeric(NA), p2=as.numeric(NA), Classification=factor(NA, levels=c('Reduced','Inconclusive','Marginal','Adequate','UnableToAnalyse')))

	## Can be done for all data:

	# For unpaired just set cov12 to 0:
	cov12 <- 0

	bnb_approximation <- 0
	bnb_strict <- .C_OLD(C_bnb_pval_wrap, sum1=checkintlimit(sum_untx), N1=checkintlimit(N1), K1=as.double(ks$K_untx), mu1=as.double(mu_untx), var1=as.double(var_untx), sum2=checkintlimit(sum_tx), N2=checkintlimit(N2), K2=as.double(ks$K_tx), mu2=as.double(mu_tx), var2=as.double(var_tx), cov12 = as.double(cov12), mean_ratio=as.double(mean_ratio), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), conjugate_priors=as.numeric(conjugate_priors), delta=checkintlimit(delta_method), beta_iters=checkintlimit(beta_iters), approx=checkintlimit(bnb_approximation), p_1=as.numeric(-1), p_2=as.numeric(-1))
	# void bnb_pval_wrap(int *sum1, int *N1, double *K1, double *mu1, double *var1, int *sum2, int *N2, double *K2, double *mu2, double *var2, double *cov12, double *mean_ratio, double *H0_1, double *H0_2, double *conjugate_priors, int *delta, int *beta_iters, int *approx, double *p_1, double *p_2);

	bnb_approximation <- 1
	bnb_auto <- .C_OLD(C_bnb_pval_wrap, sum1=checkintlimit(sum_untx), N1=checkintlimit(N1), K1=as.double(ks$K_untx), mu1=as.double(mu_untx), var1=as.double(var_untx), sum2=checkintlimit(sum_tx), N2=checkintlimit(N2), K2=as.double(ks$K_tx), mu2=as.double(mu_tx), var2=as.double(var_tx), cov12 = as.double(cov12), mean_ratio=as.double(mean_ratio), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), conjugate_priors=as.numeric(conjugate_priors), delta=checkintlimit(delta_method), beta_iters=checkintlimit(beta_iters), approx=checkintlimit(bnb_approximation), p_1=as.numeric(-1), p_2=as.numeric(-1))
	# void bnb_pval_wrap(int *sum1, int *N1, double *K1, double *mu1, double *var1, int *sum2, int *N2, double *K2, double *mu2, double *var2, double *cov12, double *mean_ratio, double *H0_1, double *H0_2, double *conjugate_priors, int *delta, int *beta_iters, int *approx, double *p_1, double *p_2);

	bnb_approximation <- 2
	bnb_approx <- .C_OLD(C_bnb_pval_wrap, sum1=checkintlimit(sum_untx), N1=checkintlimit(N1), K1=as.double(ks$K_untx), mu1=as.double(mu_untx), var1=as.double(var_untx), sum2=checkintlimit(sum_tx), N2=checkintlimit(N2), K2=as.double(ks$K_tx), mu2=as.double(mu_tx), var2=as.double(var_tx), cov12 = as.double(cov12), mean_ratio=as.double(mean_ratio), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), conjugate_priors=as.numeric(conjugate_priors), delta=checkintlimit(delta_method), beta_iters=checkintlimit(beta_iters), approx=checkintlimit(bnb_approximation), p_1=as.numeric(-1), p_2=as.numeric(-1))
	# void bnb_pval_wrap(int *sum1, int *N1, double *K1, double *mu1, double *var1, int *sum2, int *N2, double *K2, double *mu2, double *var2, double *cov12, double *mean_ratio, double *H0_1, double *H0_2, double *conjugate_priors, int *delta, int *beta_iters, int *approx, double *p_1, double *p_2);

	conjbeta <- .C_OLD(C_conjbeta_ci_wrap, N1=checkintlimit(N1), sum1=checkintlimit(sum_untx), K1=as.double(ks$K_untx), N2=checkintlimit(N2), sum2=checkintlimit(sum_tx), K2=as.double(ks$K_tx), mean_ratio=as.double(mean_ratio), iters=checkintlimit(beta_iters), conjugate_priors=as.numeric(conjugate_priors), tail=as.double(tail), ci_l=as.double(-1), ci_u=as.double(-1))
	# void conjbeta_ci_wrap(int *N1, int *sum1, double *K1, int *N2, int *sum2, double *K2, double *mean_ratio, int *iters, double *conjugate_priors, double *tail, double *ci_l, double *ci_u);

	dobson <- .C_OLD(C_dobson_ci_wrap, sum1=checkintlimit(sum_untx), sum2=checkintlimit(sum_tx), lci_c=as.double(dobson_ci_lower), uci_c=as.double(dobson_ci_upper), dobson_priors=as.double(dobson_priors), ci_l=as.double(-1), ci_u=as.double(-1))
	# void dobson_ci_wrap(int *sum1, int *sum2, double *lci_c, double *uci_c, double *dobson_priors, double *ci_l, double *ci_u);

	alloutputs[rownums['BNB unpaired (no approximation)'],c('p1','p2')] <- c(bnb_strict$p_1, bnb_strict$p_2)
	alloutputs[rownums['BNB unpaired (with approximation)'],c('p1','p2')] <- c(bnb_approx$p_1, bnb_approx$p_2)
	alloutputs[rownums['BNB unpaired'],c('p1','p2')] <- c(bnb_auto$p_1, bnb_auto$p_2)
	alloutputs[rownums['Conjugate Beta'],c('LowerCI','UpperCI')] <- c(conjbeta$ci_l, conjbeta$ci_u)
	alloutputs[rownums['Dobson'],c('LowerCI','UpperCI')] <- c(dobson$ci_l, dobson$ci_u)


	## Can only be done with <100% reduction:

	if(sum_tx > 0){
		waavp_u <- .C_OLD(C_waavp_u_ci_wrap, mu1=as.double(mu_untx), mu2=as.double(mu_tx), var1=as.double(var_untx), var2=as.double(var_tx), N1=checkintlimit(N1), N2=checkintlimit(N2), tail=as.double(tail), ci_l=as.double(-1), ci_u=as.double(-1))
		# void waavp_u_ci_wrap(double *mu1, double *mu2, double *var1, double *var2, int *N1, int *N2, double *tail, double *ci_l, double *ci_u);

		levecke_u <- .C_OLD(C_levecke_u_ci_wrap, mu1=as.double(mu_untx), mu2=as.double(mu_tx), var1=as.double(var_untx), var2=as.double(var_tx), N1=checkintlimit(N1), N2=checkintlimit(N2), tail=as.double(tail), ci_l=as.double(-1), ci_u=as.double(-1))
		# void levecke_u_ci_wrap(double *mu1, double *mu2, double *var1, double *var2, int *N1, int *N2, double *tail, double *ci_l, double *ci_u);

		mle_u <- .C_OLD(C_mle_u_ci_wrap, mu1=as.double(mu_untx), mu2=as.double(mu_tx), var1=as.double(var_untx), var2=as.double(var_tx), N1=checkintlimit(N1), N2=checkintlimit(N2), tail=as.double(tail), ci_l=as.double(-1), ci_u=as.double(-1))
		# void mle_u_ci_wrap(double *mu1, double *mu2, double *var1, double *var2, int *N1, int *N2, double *tail, double *ci_l, double *ci_u);

		alloutputs[rownums['WAAVP unpaired'],c('LowerCI','UpperCI')] <- c(waavp_u$ci_l, waavp_u$ci_u)
		alloutputs[rownums['Levecke unpaired'],c('LowerCI','UpperCI')] <- c(levecke_u$ci_l, levecke_u$ci_u)
		alloutputs[rownums['MLE unpaired'],c('LowerCI','UpperCI')] <- c(mle_u$ci_l, mle_u$ci_u)
	}else{
		alloutputs$Classification[rownums['WAAVP unpaired']] <- 'UnableToAnalyse'
		alloutputs$Classification[rownums['Levecke unpaired']] <- 'UnableToAnalyse'
		alloutputs$Classification[rownums['MLE unpaired']] <- 'UnableToAnalyse'
	}


	## Can only be done for paired data:

	if(paired){
		# Fr paired cov12 is covar:
		cov12 <- covar

		bnb_approximation <- 0
		bnb_strict <- .C_OLD(C_bnb_pval_wrap, sum1=checkintlimit(sum_untx), N1=checkintlimit(N1), K1=as.double(ks$K_untx), mu1=as.double(mu_untx), var1=as.double(var_untx), sum2=checkintlimit(sum_tx), N2=checkintlimit(N2), K2=as.double(ks$K_tx), mu2=as.double(mu_tx), var2=as.double(var_tx), cov12 = as.double(cov12), mean_ratio=as.double(mean_ratio), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), conjugate_priors=as.numeric(conjugate_priors), delta=checkintlimit(delta_method), beta_iters=checkintlimit(beta_iters), approx=checkintlimit(bnb_approximation), p_1=as.numeric(-1), p_2=as.numeric(-1))
		# void bnb_pval_wrap(int *sum1, int *N1, double *K1, double *mu1, double *var1, int *sum2, int *N2, double *K2, double *mu2, double *var2, double *cov12, double *mean_ratio, double *H0_1, double *H0_2, double *conjugate_priors, int *delta, int *beta_iters, int *approx, double *p_1, double *p_2);

		bnb_approximation <- 1
		bnb_auto <- .C_OLD(C_bnb_pval_wrap, sum1=checkintlimit(sum_untx), N1=checkintlimit(N1), K1=as.double(ks$K_untx), mu1=as.double(mu_untx), var1=as.double(var_untx), sum2=checkintlimit(sum_tx), N2=checkintlimit(N2), K2=as.double(ks$K_tx), mu2=as.double(mu_tx), var2=as.double(var_tx), cov12 = as.double(cov12), mean_ratio=as.double(mean_ratio), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), conjugate_priors=as.numeric(conjugate_priors), delta=checkintlimit(delta_method), beta_iters=checkintlimit(beta_iters), approx=checkintlimit(bnb_approximation), p_1=as.numeric(-1), p_2=as.numeric(-1))
		# void bnb_pval_wrap(int *sum1, int *N1, double *K1, double *mu1, double *var1, int *sum2, int *N2, double *K2, double *mu2, double *var2, double *cov12, double *mean_ratio, double *H0_1, double *H0_2, double *conjugate_priors, int *delta, int *beta_iters, int *approx, double *p_1, double *p_2);

		bnb_approximation <- 2
		bnb_approx <- .C_OLD(C_bnb_pval_wrap, sum1=checkintlimit(sum_untx), N1=checkintlimit(N1), K1=as.double(ks$K_untx), mu1=as.double(mu_untx), var1=as.double(var_untx), sum2=checkintlimit(sum_tx), N2=checkintlimit(N2), K2=as.double(ks$K_tx), mu2=as.double(mu_tx), var2=as.double(var_tx), cov12 = as.double(cov12), mean_ratio=as.double(mean_ratio), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), conjugate_priors=as.numeric(conjugate_priors), delta=checkintlimit(delta_method), beta_iters=checkintlimit(beta_iters), approx=checkintlimit(bnb_approximation), p_1=as.numeric(-1), p_2=as.numeric(-1))
		# void bnb_pval_wrap(int *sum1, int *N1, double *K1, double *mu1, double *var1, int *sum2, int *N2, double *K2, double *mu2, double *var2, double *cov12, double *mean_ratio, double *H0_1, double *H0_2, double *conjugate_priors, int *delta, int *beta_iters, int *approx, double *p_1, double *p_2);

		alloutputs[rownums['BNB paired (no approximation)'],c('p1','p2')] <- c(bnb_strict$p_1, bnb_strict$p_2)
		alloutputs[rownums['BNB paired (with approximation)'],c('p1','p2')] <- c(bnb_approx$p_1, bnb_approx$p_2)
		alloutputs[rownums['BNB paired'],c('p1','p2')] <- c(bnb_auto$p_1, bnb_auto$p_2)
	}


	## Can only be done for paired data and if <100% reduction:

	if(paired && sum_tx > 0){

		N <- N1
		stopifnot(N==N2)

		waavp_p <- .C_OLD(C_waavp_p_ci_wrap, mu1=as.double(mu_untx), mu2=as.double(mu_tx), var1=as.double(var_untx), var2=as.double(var_tx), cov12=as.double(covar), N=checkintlimit(N), tail=as.double(tail), ci_l=as.double(-1), ci_u=as.double(-1))
		# void waavp_p_ci_wrap(double *mu1, double *mu2, double *var1, double *var2, double *cov12, int *N, double *tail, double *ci_l, double *ci_u);

		levecke_p <- .C_OLD(C_levecke_p_ci_wrap, mu1=as.double(mu_untx), mu2=as.double(mu_tx), var1=as.double(var_untx), var2=as.double(var_tx), cov12=as.double(covar), N=checkintlimit(N), tail=as.double(tail), ci_l=as.double(-1), ci_u=as.double(-1))
		# void levecke_p_ci_wrap(double *mu1, double *mu2, double *var1, double *var2, double *cov12, int *N, double *tail, double *ci_l, double *ci_u);

		mle_p <- .C_OLD(C_mle_p_ci_wrap, mu1=as.double(mu_untx), mu2=as.double(mu_tx), var1=as.double(var_untx), var2=as.double(var_tx), cov12=as.double(covar), N=checkintlimit(N), tail=as.double(tail), ci_l=as.double(-1), ci_u=as.double(-1))
		# void mle_p_ci_wrap(double *mu1, double *mu2, double *var1, double *var2, double *cov12, int *N, double *tail, double *ci_l, double *ci_u);

		alloutputs[rownums['WAAVP paired'],c('LowerCI','UpperCI')] <- c(waavp_p$ci_l, waavp_p$ci_u)
		alloutputs[rownums['Levecke paired'],c('LowerCI','UpperCI')] <- c(levecke_p$ci_l, levecke_p$ci_u)
		alloutputs[rownums['MLE paired'],c('LowerCI','UpperCI')] <- c(mle_p$ci_l, mle_p$ci_u)

	}else{
		alloutputs$Classification[rownums['WAAVP paired']] <- 'UnableToAnalyse'
		alloutputs$Classification[rownums['Levecke paired']] <- 'UnableToAnalyse'
		alloutputs$Classification[rownums['MLE paired']] <- 'UnableToAnalyse'
	}


	##### Process and return outputs:

	ps <- which(!is.na(alloutputs$p1) & !is.na(alloutputs$p2))
	ans <- alloutputs[ps,]
	ans$Classification[ans$p1 > tail & ans$p2 <= tail] <- 'Reduced'
	ans$Classification[ans$p1 > tail & ans$p2 > tail] <- 'Inconclusive'
	ans$Classification[ans$p1 <= tail & ans$p2 <= tail] <- 'Marginal'
	ans$Classification[ans$p1 <= tail & ans$p2 > tail] <- 'Adequate'
	stopifnot(all(!is.na(ans$Classification)))
	alloutputs$Classification[ps] <- ans$Classification

	cis <- which(!is.na(alloutputs$LowerCI) & !is.na(alloutputs$UpperCI))
	ans <- alloutputs[cis,]
	ans$Classification[ans$LowerCI < H0_1 & ans$UpperCI < H0_2] <- 'Reduced'
	ans$Classification[ans$LowerCI < H0_1 & ans$UpperCI >= H0_2] <- 'Inconclusive'
	ans$Classification[ans$LowerCI >= H0_1 & ans$UpperCI < H0_2] <- 'Marginal'
	ans$Classification[ans$LowerCI >= H0_1 & ans$UpperCI >= H0_2] <- 'Adequate'
	stopifnot(all(!is.na(ans$Classification)))
	alloutputs$Classification[cis] <- ans$Classification

	return(alloutputs)

}


#### TODO
# Modify and check fecrt_power_comparison
# Website interface and fecrt_power etc could have paired type with 2 k values - ani+pre then pre+efficacy+post - use 2 new k parameters to avoid confusion and set ones not using to negative to be sure
# Modify and check fecrt_power
# Remove functions below here or leave as shims for the function above


waavp_ci <- function(pre_data, post_data, tail=0.025){

	stopifnot(length(pre_data)>=1)
	stopifnot(length(post_data)>=1)
	checksingleprob(tail)

	rv <- .C_OLD(C_waavp_ci_wrap, presum=as.integer(sum(pre_data)), predata=as.integer(pre_data), preN=as.integer(length(pre_data)), postsum=as.integer(sum(post_data)), postdata=as.integer(post_data), postN=as.integer(length(post_data)), tail=as.double(tail), ci_l=as.double(0), ci_u=as.double(0))
	# void waavp_ci_wrap(int *presum, int *predata, int *preN, int *postsum, int *postdata, int *postN, double *tail, double *ci_l, double *ci_u);

	return(c(LowerCI=rv$ci_l*100, UpeprCI=rv$ci_u*100))
}

dobson_ci <- function(pre_data, post_data, dobson_ci_lower=0.005, dobson_ci_upper=0.995, dobson_priors=c(1,1)){

	stopifnot(length(pre_data)>=1)
	stopifnot(length(post_data)>=1)
	checksingleprob(dobson_ci_lower)
	checksingleprob(dobson_ci_upper)
	stopifnot(length(dobson_priors)==2)
	stopifnot(all(dobson_priors > 0))

	rv <- .C_OLD(C_dobson_ci_wrap, presum=as.integer(sum(pre_data)), postsum=as.integer(sum(post_data)), lci_c=as.double(dobson_ci_lower), uci_c=as.double(dobson_ci_upper), dobson_priors=as.double(dobson_priors), ci_l=as.double(0), ci_u=as.double(0))
	# void dobson_ci_wrap(int *presum, int *postsum, double *lci_c, double *uci_c, double *dobson_priors, double *ci_l, double *ci_u);

	return(c(LowerCI=rv$ci_l*100, UpeprCI=rv$ci_u*100))
}

conjbeta_ci <- function(pre_data, post_data, iterations=1000, tail=0.025, conjugate_priors=c(1,1)){

	stopifnot(length(pre_data)>=1)
	stopifnot(length(post_data)>=1)
	stopifnot(length(iterations)==1)
	stopifnot(iterations > 0)
	checksingleprob(tail)

	ks <- estimate_ks(pre_data, post_data)

	rv <- .C_OLD(C_conjbeta_ci_wrap, preN=as.integer(length(pre_data)), presum=as.integer(sum(pre_data)), preK=as.double(ks$preK), postN=as.integer(length(post_data)), postsum=as.integer(sum(post_data)), postK=as.double(ks$postK), iters=as.integer(iterations), conjugate_priors=as.double(conjugate_priors), tail=as.double(tail), ci_l=as.double(0), ci_u=as.double(0))
	# void conjbeta_ci_wrap(int *preN, int *presum, double *preK, int *postN, int *postsum, double *postK, int *iters, double *conjugate_priors, double *tail, double *ci_l, double *ci_u);

	return(c(LowerCI=rv$ci_l*100, UpeprCI=rv$ci_u*100))
}

asymptotic_ci <- function(pre_data, post_data, positive_covariance=TRUE, tail=0.025){

	stopifnot(length(pre_data)>=1)
	stopifnot(length(post_data)>=1)
	stopifnot(is.logical(positive_covariance) && length(positive_covariance)==1)
	checksingleprob(tail)
	stopifnot(length(pre_data)==length(post_data))
	N <- length(pre_data)

	rv <- .C_OLD(C_asymptotic_ci_wrap, N=as.integer(N), presum=as.integer(sum(pre_data)), predata=as.integer(pre_data), postsum=as.integer(sum(post_data)), postdata=as.integer(post_data), poscov=as.integer(positive_covariance), tail=as.double(tail), u_ci_l=as.double(0), u_ci_u=as.double(0), p_ci_l=as.double(0), p_ci_u=as.double(0))
	# void asymptotic_ci_wrap(int *N, int *presum, int *predata, int *postsum, int *postdata, double *tail, double *u_ci_l, double *u_ci_u, double *p_ci_l, double *p_ci_u);

	ans <- matrix(unlist(rv[c('u_ci_l','p_ci_l','u_ci_u','p_ci_u')])*100, ncol=2, dimnames=list(c('Unpaired','Paired'), c('LowerCI','UpperCI')))
	ans[,2] <- pmin(ans[,2], 100)
	return(ans)
}

fecrt_bnb <- function(pre_data, post_data, H0_1=0.9, H0_2=0.95, edt_pre=1, edt_post=1, rep_pre=1, rep_post=1, pool_pre=1, pool_post=1, conjugate_priors=c(1,1), k_change=NA, true_k=NA, delta_method=TRUE, beta_iters=10^6){

	# Also return conjbeta CIs
	# And check that the sum of the pre and post data are both below the max int for the system otherwise error
	# Also for anything passing iters or other possibly big ints check it is below maxint

	stopifnot(length(pre_data)>=1)
	stopifnot(length(post_data)>=1)

	# Note: delta can take 3 values:  0=never, 1=unless_fails, 2=always
	delta_method <- as.numeric(delta_method)  # default is unless_fails
	stopifnot(delta_method %in% 0:2)

	stopifnot(beta_iters >= 1000)
	stopifnot(beta_iters <= 10^8)
	# Max 10^8 - 10^6 is pretty instant and gives the same results

	stopifnot(all(pre_data%%1 ==0, na.rm=TRUE))
	stopifnot(all(post_data%%1 ==0, na.rm=TRUE))
	stopifnot(length(prob_priors)==2 && all(prob_priors > 0))

	# NB edt is a threshold, i.e opposite of replicates!
	edt_change <- (rep_post*edt_pre) / (rep_pre*edt_post)
	stopifnot(length(edt_change)==1 && !is.na(edt_change) && edt_change > 0)

	# Need to adjust H0 thresholds according to the change in mean from replicates and edt:
	H0_1 <- 1 - ((1 - H0_1) * edt_change)
	H0_2 <- 1 - ((1 - H0_2) * edt_change)

	stopifnot(length(true_k)==1)
	stopifnot(length(k_change)==1)
	stopifnot(is.na(k_change) || k_change > 0)
	stopifnot(sum(pre_data) > 0)

	ks <- estimate_ks(pre_data, post_data, true_k=true_k, k_change=k_change, rep_pre=rep_pre, rep_post=rep_post, pool_pre=pool_pre, pool_post=pool_post)

	rv <- .C_OLD(C_fecrt_bnb_wrap, preN=as.integer(length(pre_data)), presum=as.integer(sum(pre_data)), preK=as.numeric(ks$preK), premean=as.numeric(mean(pre_data)), prevar=as.numeric(var(pre_data)), postN=as.integer(length(post_data)), postsum=as.integer(sum(post_data)), postK=as.numeric(ks$postK), postmean=as.numeric(mean(post_data)), postvar=as.numeric(var(post_data)), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), prob_priors=as.numeric(prob_priors), delta=as.integer(delta_method), beta_iters=as.integer(beta_iters), p_1=as.numeric(-1), p_2=as.numeric(-1), ap_1=as.numeric(-1), ap_2=as.numeric(-1))
	# void fecrt_bnb_wrap(int *preN, int *presum, double *preK, double *premean, double *prevar, int *postN, int *postsum, double *postK, double *postmean, double *postvar, double *H0_1, double *H0_2, double *prob_priors, int *delta, int *beta_iters, double *p_1, double *p_2);

	pvs <- c(H0_1p=rv$p_1, H0_2p=rv$p_2, aH0_1p=rv$ap_1, aH0_2p=rv$ap_2)
	pvs[pvs < 0] <- NA
	return(pvs)
}

bnb_ps <- function(preN, presum, preK, postN, postsum, postK, H0_1, H0_2, prob_priors, delta, beta_iters){

	rv <- .C_OLD(C_fecrt_bnb_wrap, preN=as.integer(preN), presum=as.integer(presum), preK=as.numeric(preK), postN=as.integer(postN), postsum=as.integer(postsum), postK=as.numeric(postK), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), prob_priors=as.numeric(prob_priors), delta=as.integer(delta), beta_iters=as.integer(beta_iters), p_1=as.numeric(0), p_2=as.numeric(0))
	# void fecrt_bnb_wrap(int *preN, int *presum, double *preK, int *postN, int *postsum, double *postK, double *H0_1, double *H0_2, double *prob_priors, int *delta, int *beta_iters, double *p_1, double *p_2);

	return(c(H0_1p=rv$p_1, H0_2p=rv$p_2))
}

# Keep wrappers here:
fecrt_power_comparison_wrap <- function(iterations=10^4, preN=10, postN=preN, rep_pre=1, rep_post=1, edt_pre=1, edt_post=1, paired=FALSE, premean=50, reduction=0.99, animalk=Inf, efficacyk=Inf, prek=1, postk=0.75, H0_1=0.9, H0_2=0.95, tail=0.025, prob_priors=c(1,1), delta_method=TRUE, beta_iters=10^4, positive_covariance=TRUE, dobson_ci_lower=0.005, dobson_ci_upper=0.995, dobson_priors=c(1,1)){

	# More sanity checks of inputs required
	# Forbid things to be converted to int >.Machine$integer.max

	checksingleposint(iterations)
	checksingleposint(preN)
	checksingleposint(postN)
	checksingleposint(postN)
	checksingleposint(postN)
	checksingleposint(postN)
	checksingleposint(postN)
	stopifnot(is.logical(paired) && length(paired)==1)
	checksingleposdouble(premean)
	checksingleprob(reduction)
	checksingleposdouble(animalk)
	checksingleposdouble(efficacyk)
	checksingleposdouble(prek)
	checksingleposdouble(postk)
	checksingleprob(H0_1)
	checksingleprob(H0_2)
	checksingleprob(tail)
	stopifnot(is.numeric(prob_priors) && length(prob_priors)==2 && all(prob_priors > 0))
	checksingleprob(dobson_ci_lower)
	checksingleprob(dobson_ci_upper)
	stopifnot(is.numeric(dobson_priors) && length(dobson_priors)==2 && all(dobson_priors > 0))

	checksingleposint(beta_iters)
	stopifnot(beta_iters >= 1000)
	stopifnot(beta_iters <= 10^8)
	# Max 10^8 - 10^6 is pretty instant and gives the same results
	# Also may need to be aware of .Machine$integer.max on some systems

	stopifnot(is.logical(positive_covariance) && length(positive_covariance)==1)

	# Note: delta can take 3 values:  0=never, 1=unless_fails, 2=always
	delta_method <- as.numeric(delta_method)  # default is unless_fails
	stopifnot(delta_method %in% 0:2)

	if(animalk > 100) animalk <- 100
	if(efficacyk > 100) efficacyk <- 100
	if(prek > 100) prek <- 100
	if(postk > 100) postk <- 100


	# If not paired model then rep_pre and rep_post must be 1
	# Pooling must be accounted for in R

	# Calculate pair type:
	# 0) Simple unpaired model ignores animalk and efficacyk
	# 1) Unpaired model but with replicates so needs to use animalk and efficacyk
	# 2) Paired model
	# 3) Paired model but not using replicates (suspended as little is gained by avoiding a gamma draw, and power change from 1->2 replicates is unpredictable because the distribution changes)
	if(!paired && rep_pre==1 && rep_post==1){
		pair_type <- 0
	}else if(paired){
		pair_type <- 1
	}else{
		pair_type <- 2
	}

	# Total number of classifications is 5 * number of methods:
	methods <- c('BNB', 'WAAVP', 'Dobson', 'Asymp_Paired', 'Asymp_Unpaired')
	nclass <- 5 * length(methods)

	rv <- .C_OLD(C_fecrt_power_comparison, iters=as.integer(iterations), preN=as.integer(preN), postN=as.integer(postN), maxN=as.integer(max(preN, postN)), rep_pre=as.integer(rep_pre),
	rep_post=as.integer(rep_post), edt_pre=as.double(edt_pre), edt_post=as.double(edt_post), premean=as.double(premean), reduction=as.double(reduction), pairtype=as.integer(pair_type),
	animalk=as.double(animalk), efficacyk=as.double(efficacyk), prek=as.double(prek), postk=as.double(postk), H0_1=as.double(H0_1), H0_2=as.double(H0_2), tail=as.double(tail),
	prob_priors=as.double(prob_priors), delta=as.integer(delta_method), beta_iters=as.integer(beta_iters), poscov=as.integer(positive_covariance), lci_c=as.double(dobson_ci_lower),
	uci_c=as.double(dobson_ci_upper), dobson_priors=as.double(dobson_priors), classifications=as.integer(rep(0, nclass)), obsred=as.double(rep(0, iterations)), ndeltafail=as.integer(0))
	# void fecrt_power_comparison(int *iters, int *preN, int *postN, int *maxN, int *rep_pre, int *rep_post, double *edt_pre, double *edt_post, double *premean, double *reduction, int *pairtype, double *animalk, double *efficacyk, double *prek, double *postk, double *H0_1, double *H0_2, double *tail, double *prob_priors, int *delta, int *beta_iters, int poscov, double *lci_c, double *uci_c, double *dobson_priors, int *classifications, double *obsred);

	# Results are 1: reduced, 2: inconclusive, 3: marginal, 4: adequate, 5: method failure
	# ans <- matrix(unlist(rv$classifications), ncol=5, dimnames=list(methods, c('Reduced','Inconclusive','Marginal','Adequate','Failure')), byrow=TRUE)

	ans <- data.frame(TrueEfficacy = reduction, Method = rep(methods, each=5), Classification = rep(c('Reduced','Inconclusive','Marginal','Adequate','Failure'), length(methods)), Tally=unlist(rv$classifications))

	return(ans)
}

fecrt_power_wrap <- function(iterations=10^5, use_truek=FALSE, k_change=1, preN=10, postN=preN, rep_pre=1, rep_post=1, edt_pre=1, edt_post=1, paired_data=FALSE, paired_analysis=paired_data, approximation=1, premean=50, reduction=0.99, animalk=Inf, efficacyk=Inf, prek=1, postk=0.75, H0_1=0.9, H0_2=0.95, tail=0.025, prob_priors=c(1,1), delta_method=TRUE, beta_iters=10^4){

	# More sanity checks of inputs required

	checksingleposint(iterations)
	checksinglelogical(use_truek)
	checksingleposdouble(k_change)
	checksingleposint(preN)
	checksingleposint(postN)
	checksingleposint(postN)
	checksingleposint(postN)
	checksingleposint(postN)
	checksingleposint(postN)
	stopifnot(is.logical(paired_data) && length(paired_data)==1)
	checksingleposdouble(premean)
	checksingleprob(reduction)
	checksingleposdouble(animalk)
	checksingleposdouble(efficacyk)
	checksingleposdouble(prek)
	checksingleposdouble(postk)
	checksingleprob(tail)
	stopifnot(is.numeric(prob_priors) && length(prob_priors)==2 && all(prob_priors > 0))

	checksingleposint(beta_iters)
	stopifnot(beta_iters >= 1000)
	stopifnot(beta_iters <= 10^8)
	# Max 10^8 - 10^6 is pretty instant and gives the same results

	if(length(H0_1)!=length(H0_2)){
		stop('Lengths of H0_1 and H0_2 differ', call.=FALSE)
	}
	H0_N <- length(H0_1)
	if(H0_N==1){
		checksingleprob(H0_1)
		checksingleprob(H0_2)
	}else{
		stopifnot(is.numeric(H0_1) && all(H0_1 >=0) && all(H0_1 <=1))
		stopifnot(is.numeric(H0_2) && all(H0_2 >=0) && all(H0_2 <=1))
	}


	# Note: delta can take 3 values:  0=never, 1=unless_fails, 2=always
	delta_method <- as.numeric(delta_method)  # default is unless_fails
	stopifnot(delta_method %in% 0:2)

	print('Work out a better max')
	if(animalk > 1000) animalk <- 1000
	if(efficacyk > 1000) efficacyk <- 1000
	if(prek > 1000) prek <- 1000
	if(postk > 1000) postk <- 1000


	# If not paired model then rep_pre and rep_post must be 1
	# Pooling must be accounted for in R

	# Calculate pair type:
	# 0) Simple unpaired model ignores animalk and efficacyk
	# 1) Unpaired model but with replicates so needs to use animalk and efficacyk
	# 2) Paired model
	# 3) Paired model but not using replicates (suspended as little is gained by avoiding a gamma draw, and power change from 1->2 replicates is unpredictable because the distribution changes)
	if(!paired_data && rep_pre==1 && rep_post==1){
		pair_type <- 0
	}else if(!paired_data){
		pair_type <- 1
	}else{
		pair_type <- 2
	}

	paired_analysis <- as.integer(paired_analysis)

	rv <- .C_OLD(C_fecrt_power, iters=as.integer(iterations), paired_analysis=as.integer(paired_analysis), use_truek=as.integer(use_truek), kchange=as.double(k_change),
	preN=as.integer(preN), postN=as.integer(postN), maxN=as.integer(max(preN, postN)), rep_pre=as.integer(rep_pre), rep_post=as.integer(rep_post),
	edt_pre=as.double(edt_pre), edt_post=as.double(edt_post), premean=as.double(premean), reduction=as.double(reduction), pairtype=as.integer(pair_type),
	animalk=as.double(animalk), efficacyk=as.double(efficacyk), prek=as.double(prek), postk=as.double(postk), H0_1=as.double(H0_1), H0_2=as.double(H0_2),
	H0_N=as.integer(H0_N), prob_priors=as.double(prob_priors), delta=as.integer(delta_method), beta_iters=as.integer(beta_iters), approximation=as.integer(approximation),
	p_1=as.double(rep(0, iterations*H0_N)), p_2=as.double(rep(0, iterations*H0_N)), obsred=as.double(rep(0, iterations)), ndeltafail=as.integer(0))
	# void fecrt_power(int *iters, int *paired_analysis, int *use_truek, double *kchange, int *preN, int *postN, int *maxN, int *rep_pre, int *rep_post, double *edt_pre, double *edt_post, double *premean, double *reduction, int *pairtype, double *animalk, double *efficacyk, double *prek, double *postk, double *H0_1, double *H0_2, int *H0_N, double *prob_priors, int *delta, int *beta_iters, int *approximation, double *p_1, double *p_2, double *obsred, int *ndeltafail);

	stopifnot(length(rv$obsred) == iterations)
	stopifnot((iterations*H0_N) == length(rv$p_1))
	ans <- data.frame(H0_1=rep(H0_1, each=iterations), H0_2=rep(H0_2, each=iterations), TrueEfficacy=rep(reduction, iterations*H0_N), ObsEfficacy=rep(rv$obsred, each=H0_N), p_1=rv$p_1, p_2=rv$p_2, Classification=factor('Failure', levels=c('Reduced','Inconclusive','Marginal','Adequate','Failure')))
	ans$Classification[ans$p_1 > tail & ans$p_2 <= tail] <- 'Reduced'
	ans$Classification[ans$p_1 > tail & ans$p_2 > tail] <- 'Inconclusive'
	ans$Classification[ans$p_1 <= tail & ans$p_2 <= tail] <- 'Marginal'
	ans$Classification[ans$p_1 <= tail & ans$p_2 > tail] <- 'Adequate'

	return(ans)
}



##### OLDER


# Consolidate power functions for different tasks i.e. pair types??  Probably put elsewhere.
fecrt_power <- function(reduction, pair_type=0, preN=20, postN=20, poolsize_pre=1, poolsize_post=1, premean=10, animalk=2, efficacyk=2, prek=3, postk=prek, rep_pre=1, rep_post=1, edt_pre=1, edt_post=1, target=0.95, tol=0.05, tail=0.025, iterations=10000, prob_priors=c(1,1), delta_method=1, beta_iters=10^4){

	stopifnot(length(reduction)==1 && reduction <= 1 && reduction > -Inf)
	stopifnot(length(pair_type)==1 && pair_type >= 0)
	stopifnot(length(preN)==1 && preN > 0 && round(preN)==preN)
	stopifnot(length(postN)==1 && postN > 0 && round(postN)==postN)
	stopifnot(length(poolsize_pre)==1 && poolsize_pre > 0 && round(poolsize_pre)==poolsize_pre)
	stopifnot(length(poolsize_post)==1 && poolsize_post > 0 && round(poolsize_post)==poolsize_post)
	stopifnot(length(premean)==1 && premean > 0)
	stopifnot(length(animalk)==1 && animalk > 0)
	stopifnot(length(prek)==1 && prek > 0)
	stopifnot(length(postk)==1 && postk > 0)
	stopifnot(length(efficacyk)==1 && efficacyk > 0)
	stopifnot(length(rep_pre)==1 && rep_pre > 0 && round(rep_pre)==rep_pre)
	stopifnot(length(rep_post)==1 && rep_post > 0 && round(rep_post)==rep_post)
	stopifnot(length(edt_pre)==1 && edt_pre > 0)
	stopifnot(length(edt_post)==1 && edt_post > 0)
	stopifnot(length(target)==1 && target <= 1 && target > -Inf)
	stopifnot(length(tol)==1 && tol < 1)
	stopifnot(length(tail)==1 && tail > 0 && tail < 1)
	stopifnot(length(iterations)==1 && iterations > 0 && round(iterations)==iterations)
	stopifnot(length(prob_priors)==2 && is.numeric(prob_priors) && all(prob_priors > 0))

	# Note: delta can take 3 values:  0=never, 1=unless_fails, 2=always
	stopifnot(delta_method %in% 0:2)
	if(delta_method == 0){
		cat("Note:  Not using the delta method approximation at all will dramatically increase computation time\n")
	}
	beta_iters <- as.integer(beta_iters)
	stopifnot(beta_iters >= 1000)

	H0_1 <- target-tol
	H0_2 <- target
	reduction

	# The paired int is either:
	# 0:  no pairing, no animalk, no efficacyk, can't use replicates or pooling
	# 1:  no pairing, use animalk and efficacyk so we can use replicates or pooling
	# 2:  with pairing, use animalk and efficacyk so we can use replicates or pooling
	# suspended:  with pairing, use animalk but no efficacyk, can't use replicates or pooling - disadvantage is distribution change from combining gammas

	if(pair_type == 0){

		stopifnot(rep_pre==1 && rep_post==1)

		# C++ function does not combine variances so need to do this if animalk AND efficacyk not supplied as Inf:
		if(animalk == Inf || efficacyk == Inf){
			stopifnot(animalk==Inf && efficacyk==Inf)

			# Avoid NA/Inf in foreign call:
			animalk <- 0
			efficacyk <- 0
		}else{
			prek <- (animalk * prek) / (animalk + prek + 1)
			intk <- (animalk * efficacyk) / (animalk + efficacyk + 1)
			postk <- (intk * postk) / (intk + postk + 1)
			# Probably need a combine_k function!
		}

	}else if(pair_type == 1){

		# C++ function combines animal and efficacy variability sources for treatment group

	}else if(pair_type == 2){

		# All variability sources are separate

	}else if(pair_type == 3){

		stopifnot(rep_pre==1 && rep_post==1)

		# C++ function combines efficacyk and postk
		stop('pair_type is disabled')

	}else{
		stop('Unrecognised pair_type')
	}

	# Pooling works by reducing k sources:
	if(poolsize_pre > 1 || poolsize_post > 1){

		# NOTE that in the case of pooling, preN and postN becomes number of pools!

		# Currently only available for simple unpaired:
		stopifnot(pair_type == 0)

		prek <- prek * poolsize_pre
		postk <- postk * poolsize_post

		# When taking into account rep_pre and rep_post will need to adjust k and then set rep to 1
		# (assumes that repeat samples from the same animal are always pooled together)
		# Alternatively, rep_pre and rep_post could represent repeat samples from the same pool...

	}

	# Currently only 5*2 but will be 5*4:
	classifications <- 20

	results <- .C_OLD(C_fecrt_power_comparison, iters=as.integer(iterations), preN=as.integer(preN), postN=as.integer(postN), maxN=as.integer(max(preN, postN)), rep_pre=as.integer(rep_pre),
	rep_post=as.integer(rep_post), edt_pre=as.double(edt_pre), edt_post=as.double(edt_post), premean=as.double(premean), reduction=as.double(reduction), pair_type=as.integer(pair_type),
	animalk=as.double(animalk), efficacyk=as.double(efficacyk), prek=as.double(prek), postk=as.double(postk), H0_1=as.double(H0_1), H0_2=as.double(H0_2), tail=as.double(tail),
	prob_priors=as.double(prob_priors), delta=as.integer(delta_method), beta_iters=as.integer(beta_iters), predata=as.integer(rep(0, max(preN, postN))), postdata=as.integer(rep(0, max(preN, postN))), classifications=as.integer(rep(0, classifications)),
	obsred=as.double(rep(0, iterations)))

	return(results)

}

fecrt_power_paired <- function(reduction, preN=20, postN=20, animalk=2, efficacyk=2, prek=3, postk=prek, premean=10, rep_pre=1, rep_post=1, edt_pre=1, edt_post=1, target=0.95, tol=0.05, tail=0.025, iterations=1000, prob_priors=c(1,1)){

	stopifnot(length(reduction)==1 && reduction <= 1 && reduction > -Inf)
	stopifnot(length(preN)==1 && preN > 0 && round(preN)==preN)
	stopifnot(length(postN)==1 && postN > 0 && round(postN)==postN)
	stopifnot(length(animalk)==1 && animalk > 0)
	stopifnot(length(prek)==1 && prek > 0)
	stopifnot(length(postk)==1 && postk > 0)
	stopifnot(length(efficacyk)==1 && efficacyk > 0)
	stopifnot(length(premean)==1 && premean > 0)
	stopifnot(length(rep_pre)==1 && rep_pre > 0 && round(rep_pre)==rep_pre)
	stopifnot(length(rep_post)==1 && rep_post > 0 && round(rep_post)==rep_post)
	stopifnot(length(edt_pre)==1 && edt_pre > 0)
	stopifnot(length(edt_post)==1 && edt_post > 0)
	stopifnot(length(target)==1 && target <= 1 && target > -Inf)
	stopifnot(length(tol)==1 && tol < 1)
	stopifnot(length(tail)==1 && tail > 0 && tail < 1)
	stopifnot(length(iterations)==1 && iterations > 0 && round(iterations)==iterations)
	stopifnot(length(prob_priors)==2 && is.numeric(prob_priors) && all(prob_priors > 0))

	H0_1 <- target-tol
	H0_2 <- target
	reduction

	# The paired int is either:
	# 0:  no pairing, no animalk, no efficacyk, can't use replicates or pooling
	# 1:  no pairing, use animalk and efficacyk so we can use replicates or pooling
	# 2:  with pairing, use animalk but no efficacyk, can't use replicates or pooling
	# 3:  with pairing, use animalk and efficacyk so we can use replicates or pooling

	if(rep_pre > 1 || rep_post > 1){
		# The C++ function will combine efficacy and postk
		paired <- 2
	}else{
		paired <- 3
	}


	# Not using but keep for info:
	if(FALSE && paired==0){
		prek <- (animalk * prek) / (animalk + prek + 1)
		intk <- (animalk * efficacyk) / (animalk + efficacyk + 1)
		postk <- (intk * postk) / (intk + postk + 1)
	}

	# Currently only 5*2 but will be 5*4:
	classifications <- 20

	results <- .C_OLD(C_fecrt_power_comparison, iters=as.integer(iterations), preN=as.integer(preN), postN=as.integer(postN), maxN=as.integer(max(preN, postN)), rep_pre=as.integer(rep_pre),
	rep_post=as.integer(rep_post), edt_pre=as.double(edt_pre), edt_post=as.double(edt_post), premean=as.double(premean), reduction=as.double(reduction), paired=as.integer(paired),
	animalk=as.double(animalk), efficacyk=as.double(efficacyk), prek=as.double(prek), postk=as.double(postk), H0_1=as.double(H0_1), H0_2=as.double(H0_2), tail=as.double(tail),
	prob_priors=as.double(prob_priors), predata=as.integer(rep(0, max(preN, postN))), postdata=as.integer(rep(0, max(preN, postN))), classifications=as.integer(rep(0, classifications)),
	obsred=as.double(rep(0, iterations)))

	return(results)

}

fecrt_power_unpaired <- function(reduction, controlN=20, treatmentN=20, controlk=1, treatmentk=0.5, controlmean=10, edt_control=1, edt_treatment=1, target=0.95, tol=0.05, tail=0.025, iterations=1000, prob_priors=c(1,1)){

	paired <- FALSE
	# No replicates for simple unpaired:

	stopifnot(length(reduction)==1 && reduction <= 1 && reduction > -Inf)
	stopifnot(length(paired)==1 && is.logical(paired))
	stopifnot(length(controlN)==1 && controlN > 0 && round(controlN)==controlN)
	stopifnot(length(treatmentN)==1 && treatmentN > 0 && round(treatmentN)==treatmentN)
	stopifnot(length(controlk)==1 && controlk > 0)
	stopifnot(length(treatmentk)==1 && treatmentk > 0)
	stopifnot(length(controlmean)==1 && controlmean > 0)
	stopifnot(length(rep_pre)==1 && rep_pre > 0 && round(rep_pre)==rep_pre)
	stopifnot(length(rep_post)==1 && rep_post > 0 && round(rep_post)==rep_post)
	stopifnot(length(edt_control)==1 && edt_control > 0)
	stopifnot(length(edt_treatment)==1 && edt_treatment > 0)
	stopifnot(length(target)==1 && target <= 1 && target > -Inf)
	stopifnot(length(tol)==1 && tol < 1)
	stopifnot(length(tail)==1 && tail > 0 && tail < 1)
	stopifnot(length(iterations)==1 && iterations > 0 && round(iterations)==iterations)
	stopifnot(length(prob_priors)==2 && is.numeric(prob_priors) && all(prob_priors > 0))

	H0_1 <- target-tol
	H0_2 <- target

	# Currently only 5*2 but will be 5*4:
	classifications <- 20

	results <- .C_OLD(C_fecrt_power_comparison, iters=as.integer(iterations), controlN=as.integer(controlN), treatmentN=as.integer(treatmentN), maxN=as.integer(max(controlN, treatmentN)), rep_pre=as.integer(rep_pre), rep_post=as.integer(rep_post), edt_control=as.double(edt_control), edt_treatment=as.double(edt_treatment), controlmean=as.double(controlmean), reduction=as.double(reduction), paired=as.integer(paired), animalk=as.double(Inf), efficacyk=as.double(Inf), prek=as.double(controlk), postk=as.double(treatmentk), H0_1=as.double(H0_1), H0_2=as.double(H0_2), tail=as.double(tail), prob_priors=as.double(prob_priors), predata=as.integer(rep(0, max(controlN, treatmentN))), postdata=as.integer(rep(0, max(controlN, treatmentN))), classifications=as.integer(rep(0, classifications)), obsred=as.double(rep(0, iterations)))

	return(results)

}

fecrt_power_pooled <- function(reduction, pools_control=5, pools_treatment=5, Nperpool_control=4, Nperpool_treatment=4, controlk=1, treatmentk=0.5, controlmean=10, edt_control=1, edt_treatment=1, target=0.95, tol=0.05, tail=0.025, iterations=1000, prob_priors=c(1,1)){

	# Pools is for unpaired data only currently
	# (but can be extended to unpaired just by adjusting all k values ???  Check by simulation)
	paired <- FALSE

	# No replicates for pools currently
	# (although theoretically justifiable as long as all replicates from the same animal in)
	rep_pre <- 1
	rep_post <- 1

	#

	stopifnot(length(reduction)==1 && reduction <= 1 && reduction > -Inf)
	stopifnot(length(paired)==1 && is.logical(paired))
	stopifnot(length(controlN)==1 && controlN > 0 && round(controlN)==controlN)
	stopifnot(length(treatmentN)==1 && treatmentN > 0 && round(treatmentN)==treatmentN)
	stopifnot(length(controlk)==1 && controlk > 0)
	stopifnot(length(treatmentk)==1 && treatmentk > 0)
	stopifnot(length(controlmean)==1 && controlmean > 0)
	stopifnot(length(rep_pre)==1 && rep_pre > 0 && round(rep_pre)==rep_pre)
	stopifnot(length(rep_post)==1 && rep_post > 0 && round(rep_post)==rep_post)
	stopifnot(length(edt_control)==1 && edt_control > 0)
	stopifnot(length(edt_treatment)==1 && edt_treatment > 0)
	stopifnot(length(target)==1 && target <= 1 && target > -Inf)
	stopifnot(length(tol)==1 && tol < 1)
	stopifnot(length(tail)==1 && tail > 0 && tail < 1)
	stopifnot(length(iterations)==1 && iterations > 0 && round(iterations)==iterations)
	stopifnot(length(prob_priors)==2 && is.numeric(prob_priors) && all(prob_priors > 0))

	H0_1 <- target-tol
	H0_2 <- target

	# Currently only 5*2 but will be 5*4:
	classifications <- 20

	results <- .C_OLD(C_fecrt_power_comparison, iters=as.integer(iterations), controlN=as.integer(controlN), treatmentN=as.integer(treatmentN), maxN=as.integer(max(controlN, treatmentN)), rep_pre=as.integer(rep_pre), rep_post=as.integer(rep_post), edt_control=as.double(edt_control), edt_treatment=as.double(edt_treatment), controlmean=as.double(controlmean), reduction=as.double(reduction), paired=as.integer(paired), animalk=as.double(Inf), efficacyk=as.double(Inf), prek=as.double(controlk), postk=as.double(treatmentk), H0_1=as.double(H0_1), H0_2=as.double(H0_2), tail=as.double(tail), prob_priors=as.double(prob_priors), predata=as.integer(rep(0, max(controlN, treatmentN))), postdata=as.integer(rep(0, max(controlN, treatmentN))), classifications=as.integer(rep(0, classifications)), obsred=as.double(rep(0, iterations)))

	return(results)

}
