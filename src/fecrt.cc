/*
This file is part of the bayescount package:  https://cran.r-project.org/package=bayescount
Copyright of Matthew James Denwood, licensed GPL>=2
The code is for analysis and power calculations for FECRT studies
General note: long long is used for all sum1 and sum2 in case int is too short on
some systems - for this reason CXX_STD = CXX11 is specified in the Makevars file
[search https://cran.r-project.org/doc/manuals/r-release/R-exts.html for "long long"]
Last updated by MJD June 2018
*/

#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rcpp.h>

// For std::nth_element used in ratio of betas MC approximation:
#include <vector>
#include <algorithm>
#include <functional>

#include "dists.h"

#include "fecrt.h"


/**************************************************************/
/*   Beta negative binomial functions (not exported to R)     */
/**************************************************************/

/*
double pbnb_lower(long long q, double bnb_k, double bnb_alpha, double bnb_beta, bool inclusive){
	
	// Convert from beta-NB to generalised hypergeometric:
	double ghg_a = -bnb_beta;
	double ghg_k = -bnb_k;
	double ghg_N = bnb_alpha - 1.0;
	double p;
	
	// If an inclusive probability then redefine q and return 0 if necessary:
	if(inclusive){
		q--;
		if(q < 0L){
			return(0.0);
		}
	}
	
	hyperType variety = typeHyper(ghg_a, ghg_k, ghg_N);
	if (! checkHyperArgument(q, ghg_a, ghg_k, ghg_N, variety))
		p = NA_REAL;
	else if (variety==classic)
		p = phypergeometric(q, (int) ghg_a, (int) ghg_k, (int) ghg_N);
	else
		p = pgenhypergeometric(q, ghg_a, ghg_k, ghg_N, variety);

	return(p);
}
double pbnb_upper(long long q, double bnb_k, double bnb_alpha, double bnb_beta, bool inclusive){
	
	// Convert from beta-NB to generalised hypergeometric:
	double ghg_a = -bnb_beta;
	double ghg_k = -bnb_k;
	double ghg_N = bnb_alpha - 1.0;
	double p;
	
	// We redefine upper tail as inclusive, so if q=0:
	if(q==0L){
		return(1.0);
	}
	
	hyperType variety = typeHyper(ghg_a, ghg_k, ghg_N);
	if (! checkHyperArgument(q, ghg_a, ghg_k, ghg_N, variety))
		p = NA_REAL;
	else if (variety==classic)
		p = 1.0 - phypergeometric(q-1L, (int) ghg_a, (int) ghg_k, (int) ghg_N);
	else
		p = 1.0 - pgenhypergeometric(q-1L, ghg_a, ghg_k, ghg_N, variety);
	// Note q-1 here to make p inclusive
	
	return(p);
}

*/

double pbnb_lower(long long q, double bnb_k, double bnb_alpha, double bnb_beta){
	
	// Convert from beta-NB to generalised hypergeometric:
	double ghg_a = -bnb_beta;
	double ghg_k = -bnb_k;
	double ghg_N = bnb_alpha - 1.0;
	double p;
	
	hyperType variety = typeHyper(ghg_a, ghg_k, ghg_N);
	if (! checkHyperArgument(q, ghg_a, ghg_k, ghg_N, variety))
		p = NA_REAL;
	else if (variety==classic)
		p = phypergeometric(q, (int) ghg_a, (int) ghg_k, (int) ghg_N);
	else
		p = pgenhypergeometric(q, ghg_a, ghg_k, ghg_N, variety);

	return(p);
}

double pbnb_upper(long long q, double bnb_k, double bnb_alpha, double bnb_beta){
	
	// Convert from beta-NB to generalised hypergeometric:
	double ghg_a = -bnb_beta;
	double ghg_k = -bnb_k;
	double ghg_N = bnb_alpha - 1.0;
	double p;
	
	// We redefine upper tail as inclusive, so if q=0:
	if(q==0L){
		return(1.0);
	}
	
	hyperType variety = typeHyper(ghg_a, ghg_k, ghg_N);
	if (! checkHyperArgument(q, ghg_a, ghg_k, ghg_N, variety))
		p = NA_REAL;
	else if (variety==classic)
		p = 1.0 - phypergeometric(q-1L, (int) ghg_a, (int) ghg_k, (int) ghg_N);
	else
		p = 1.0 - pgenhypergeometric(q-1L, ghg_a, ghg_k, ghg_N, variety);
	// Note q-1 here to make p inclusive
	
	return(p);
}

double pbnb_2(int q, double k, double alpha, double beta, bool lower, bool inclusive){
	
	double p = NA_REAL;
	
	if(lower){
		if(!inclusive){
			if(q==0L){
				return 0.0;
			}
			q--;
		}
		p = pbnb_lower(q, k, alpha, beta);
	}else{
		if(inclusive){
			if(q==0L){
				return 1.0;
			}
			q--;
		}
		p = 1.0 - pbnb_lower(q, k, alpha, beta);
	}
	
	return p;
}


/**************************************************************/
/*   Derivative and utility functions (not exported to R)     */
/**************************************************************/

// Simple function to convert mean and variance to alpha and beta:
void beta_params(double mu, double s2, double &alpha, double &beta){
	
	double mu1m = mu * (1.0 - mu);

	// Limit s2 to be fractionally greater than mu * (1-mu):
	if(s2 >= mu1m){
//		printf("Adjusting s2 from %f ", s2);
		s2 = 0.999999 * mu1m;
//		printf("to %f\n", s2);
	}

	double ab = mu1m / s2 -1;
	alpha = mu * ab;
	beta = (1-mu)*ab;
	
}

// Non-linear transformation function:
double g_fun(double p, double r, double s, double t){
	return(p*r*t / (p*r*t - p*s + s));
}

// First 4 derivatives of g:
double gp1_fun(double p, double r, double s, double t){
	return(r*s*t / std::pow((p-1.0)*s - p*r*t, 2));
}
double gp2_fun(double p, double r, double s, double t){
	return(-2.0*r*s*t*(s - r*t) / std::pow((p-1.0)*s - p*r*t, 3));
}
double gp3_fun(double p, double r, double s, double t){
	return(6.0*r*s*t* std::pow(s - r*t, 2) / std::pow(p*r*t - p*s + s, 4));
}
double gp4_fun(double p, double r, double s, double t){
	return(24.0*r*s*t* std::pow(r*t - s, 3) / std::pow((p-1.0)*s - p*r*t, 5));
}

// Functions for expectations of powers 2-5 of a Beta distribution with parameters a(lpha), b(eta):
double epf2(double a, double b){
	return(((a+1.0)*a)  / ((a+b+1.0)*(a+b)));
}
double epf3(double a, double b){
	return(((a+2.0)*(a+1.0)*a)  / ((a+b+2.0)*(a+b+1.0)*(a+b)));
}
double epf4(double a, double b){
	return(((a+3.0)*(a+2.0)*(a+1.0)*a)  / ((a+b+3.0)*(a+b+2.0)*(a+b+1.0)*(a+b)));
}
double epf5(double a, double b){
	return(((a+4.0)*(a+3.0)*(a+2.0)*(a+1.0)*a)  / ((a+b+4.0)*(a+b+3.0)*(a+b+2.0)*(a+b+1.0)*(a+b)));
}

// The binomial expansions for E((x - mu)^t) for t in 2:5:
double expow2(double a, double b){
	double m = a / (a+b);
	return(epf2(a,b) - std::pow(m,2));
}
double expow3(double a, double b){
	double m = a / (a+b);
	return(epf3(a,b) - 3.0*m*epf2(a,b) + 2*std::pow(m,3));
}
double expow4(double a, double b){
	double m = a / (a+b);
	return(epf4(a,b) - 4.0*m*epf3(a,b) + 6.0*std::pow(m,2)*epf2(a,b) - 3.0*std::pow(m,4));
}
double expow5(double a, double b){
	double m = a / (a+b);
	return(epf5(a,b) - 5.0*m*epf4(a,b) + 10.0*std::pow(m,2)*epf3(a,b) - 10.0*std::pow(m,3)*epf2(a,b) + 4.0*std::pow(m,5));
}


// Delta method approximation to the change in mean (first 2 terms of Taylor series):
double delta_mean(double p_mu, double p_var, double r, double s, double t){
	
	// p_mu and p_var are calculated by calling function (for efficiency) as:
	// double p_mu = alpha / (alpha + beta);
	// double p_var = (alpha * beta) / (std::pow(alpha + beta, 2) * (alpha + beta + 1));

	double rm = g_fun(p_mu, r, s, t) + 0.5 * gp2_fun(p_mu, r, s, t) * p_var;
	return(rm);
}

// Delta method approximation to the change in variance (higher order Taylor series):
double delta_var(double alpha, double beta, double p_mu, double p_var, double r, double s, double t){
	
	// p_mu and p_var are calculated by calling function (for efficiency) as:
	// double p_mu = alpha / (alpha + beta);
	// double p_var = (alpha * beta) / (std::pow(alpha + beta, 2) * (alpha + beta + 1));
	
	double g1 = gp1_fun(p_mu, r, s, t);
	double g2 = gp2_fun(p_mu, r, s, t);
	double g3 = gp3_fun(p_mu, r, s, t);
	double g4 = gp4_fun(p_mu, r, s, t);
	
	double rv = std::pow(g1,2) * p_var +
				2.0 * g1 * g2/2.0 * expow3(alpha, beta) +
				(std::pow(g2,2)/4.0 + 2.0*g1 * g3/6.0) * expow4(alpha, beta) +
				(2.0*g1 * g4/24.0 + g2 * g3/6.0) * expow5(alpha, beta);
	return(rv);
}


/*********************************************************************/
/*   Underlying C++ hypothesis test function (not exported to R)     */
/*********************************************************************/

int bnb_pval(long long sum1, int N1, double K1, double mu1, double var1, long long sum2, int N2, double K2, double mu2, double var2, double cov12, double mean_ratio, double H0_1, double H0_2, double *conjugate_priors, int delta, int beta_iters, int approx, double *p_1, double *p_2){
	
	// Note: delta can take 3 values:  0=never, 1=unless_fails, 2=always
	// Note: approx can take 3 values:  0=never, 1=if_necessary, 2=always
	
	// This same function is used for paired and unpaired - estimate_ks sorts out the correction
	
	// To force K values:
	// K1 = 24.0;
	// K2 = 15.0;
	
	// Calculate beta parameters:
	double alpha1 = (double) sum1 + conjugate_priors[1];
	double beta1 = (double) N1 * K1 + conjugate_priors[0];
	
	double tmu = alpha1 / (alpha1 + beta1);
	double tvar = (alpha1 * beta1) / (std::pow(alpha1 + beta1, 2L) * (alpha1 + beta1 + 1.0));
	
	// effK takes account of the change in replicates and/or edt:
	double effK = K2 * mean_ratio * (double) N2;
	
	/*  This is now done in the Rcpp wrapper as no other C function calls this with delta=1	
	if(delta > 0){
		GetRNGstate();
	}
	*/
	
	bool dfailed = false;
	
	// Calculate the observed reduction
	double obsred = 1.0 - (mu2 / mu1);
	// Note: replicates and edt is dealt with by effK
	
	// First hypothesis:

	// Only calculate the pvalue when it makes sense to do so:
	if(obsred < H0_1){

		*p_1 = 1.0;

	}else{
		
		double meanchange = (1.0 - H0_1);
		if(meanchange < 0.0){
			meanchange = 0.0;
		}
		
		double newEprob = 0.0;
		double newVprob = 0.0;
		
		bool intensive=true;
		// Note: delta can take 3 values:  0=never, 1=unless_fails, 2=always
		if(delta > 0){
			
			newEprob = delta_mean(tmu, tvar, K1, K2, meanchange);
			newVprob = delta_var(alpha1, beta1, tmu, tvar, K1, K2, meanchange);
			
			// If we get a negative or zero variance then fall back to not using the delta method:
			intensive = delta > 0 && (newEprob <= 0 || newVprob <= 0);
			dfailed = newEprob <= 0 || newVprob <= 0;
		}

		if(intensive){
			double rmean=0;
			double rvar=0;
			double deltaval=0;
			for(int i=1; i<=beta_iters; i++){
				deltaval = g_fun(R::rbeta(alpha1, beta1), K1, K2, meanchange) - rmean;
				rmean += deltaval / (double) i;
				rvar += deltaval*deltaval;
			}
			newEprob = rmean;
			newVprob = rvar / (double) (beta_iters-1);
			
		}
		
		double alpha2 = 0.0;
		double beta2 = 0.0;
		beta_params(newEprob, newVprob, alpha2, beta2);
		
		// Note alpha and beta swapped as BNB is failures before successes and we want vice versa!
		*p_1 = pbnb_lower(sum2, effK, beta2, alpha2);
		// pbnb_lower is inclusive probability
		
		// Only use approximation if failing to calculate p-value
		// Note: approx can take 3 values:  0=never, 1=if_necessary, 2=always
		if(approx > 0 && (approx > 1 || !R_finite(*p_1))){
			
			if(beta2 <= 2.0){
				warning("Bad approximation due to small alpha2: ", beta2);
			}
			// Variance and mean of beta negative binomial - alpha and beta are swapped:
			double bb = alpha2;
			double aa = beta2;
			double totalvar = ((effK*(aa+effK-1.0)*bb*(aa+bb-1.0)) / ((aa-2.0)*(aa-1.0)*(aa-1.0)));
			double totalmean = (effK*bb) / (aa-1.0);
				
			// Get equivalent gamma parameters:
			double app_scale = totalvar/totalmean;
			double app_shape = totalmean/app_scale;
			
		    // Signature:  double pgamma (double x, double alph, double scale, int lower_tail, int log_p)
			*p_1 = R::pgamma(double(sum2)+0.5, app_shape, app_scale, true, false);
			// We add 0.5 here as we assume the double counts are rounded to int
			
		}
		
	}
	
	// Second hypothesis

	// Only calculate the pvalue when it makes sense to do so:
	if(obsred >= H0_2 || sum2 == 0){

		*p_2 = 1.0;

	}else{
		
		double meanchange = (1.0 - H0_2);
		if(meanchange < 0.0){
			meanchange = 0.0;
		}
		
		double newEprob = 0.0;
		double newVprob = 0.0;
		
		bool intensive=true;
		// Note: delta can take 3 values:  0=never, 1=unless_fails, 2=always
		if(delta > 0){
			
			newEprob = delta_mean(tmu, tvar, K1, K2, meanchange);
			newVprob = delta_var(alpha1, beta1, tmu, tvar, K1, K2, meanchange);
			
			// If we get a negative or zero variance then fall back to not using the delta method:
			intensive = delta > 0 && (newEprob <= 0 || newVprob <= 0);
			dfailed = dfailed || newEprob <= 0 || newVprob <= 0;
		}

		if(intensive){
			double rmean=0;
			double rvar=0;
			double deltaval=0;
			for(int i=1; i<=beta_iters; i++){
				deltaval = g_fun(R::rbeta(alpha1, beta1), K1, K2, meanchange) - rmean;
				rmean += deltaval / (double) i;
				rvar += deltaval*deltaval;
			}
			newEprob = rmean;
			newVprob = rvar / (double) (beta_iters-1);
		}
	
		double alpha2 = 0.0;
		double beta2 = 0.0;
		beta_params(newEprob, newVprob, alpha2, beta2);
		
		// Note alpha and beta swapped as BNB is failures before successes and we want vice versa!
		*p_2 = pbnb_upper(sum2, effK, beta2, alpha2);
		// pbnb_upper is inclusive probability - no deduction of 1L needed here

		// Only use approximation if failing to calculate p-value
		// Note: approx can take 3 values:  0=never, 1=if_necessary, 2=always
		if(approx > 0 && (approx > 1 || !R_finite(*p_2))){
			
			if(beta2 <= 2.0){
				warning("Bad approximation due to small beta2: ", beta2);
			}
			
			// Variance and mean of beta negative binomial - alpha and beta are swapped:
			double bb = alpha2;
			double aa = beta2;
			double totalvar = ((effK*(aa+effK-1.0)*bb*(aa+bb-1.0)) / ((aa-2.0)*(aa-1.0)*(aa-1.0)));
			double totalmean = (effK*bb) / (aa-1.0);
				
			// Get equivalent gamma parameters:
			double app_scale = totalvar/totalmean;
			double app_shape = totalmean/app_scale;
			
		    // Signature:  double pgamma (double x, double alph, double scale, int lower_tail, int log_p)
			*p_2 = R::pgamma(double(sum2)-1.5, app_shape, app_scale, false, false);
			// We take 1.5 off here because (a) sum2 must be LESS than expected to be inferior, and (b) assume the double counts are rounded to int

		}
		
	}
	
	
	
	// Check to see if we need to use a backup approximation:
	// Note: approx can take 3 values:  0=never, 1=if_necessary, 2=always
	
	// TODO: fix approximation:
	// if(approx > 0 && (approx > 1 || ! R_finite(*p_1) || ! R_finite(*p_2))){
	if(approx > 0 && (! R_finite(*p_1) || ! R_finite(*p_2))){
		
		// Rcpp::warning() gives a segfault!?!??!?!:
		warning("Incorrect approximation called");
		
		// TODO: either remove this or fall back on an MLE approximation
		
		// Bail if zero variance:
		if(var1 <= 0 && var2 <= 0){
			// Error handling:  https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Error-signaling
			if(approx == 3){
				warning("Unable to evaluate approximate p-values with zero variance in one or other dataset");	
			}else{
				warning("Both the BNB and approximations failed (the latter due to zero variance in one or other dataset)");	
			}
			*p_1 = NA_REAL;
			*p_2 = NA_REAL;
			return (int) dfailed;
		}
		
		double red_mu = mu2 / mu1;

		// Use N2 to calculate the variance:
		double Nd = (double) N2;
		double red_se = std::sqrt((std::pow(mu2, 2) / std::pow(mu1, 4) * var1 + std::pow(mu1, -2) * mu2 + 2.0 * mu2 / std::pow(mu1, 3) * cov12) / Nd);
		
		if(approx > 1 || !R_finite(*p_1)){
			double meanchange = (1.0 - H0_1);
			if(meanchange < 0.0){
				meanchange = 0.0;
			}
			// Signature: pnorm(double x, double mu, double sigma, int lower_tail, int log_p)
			*p_1 = R::pnorm(meanchange,red_mu,red_se,false,false);
		}
		
		if(approx > 1 || !R_finite(*p_2)){
			double meanchange = (1.0 - H0_2);
			if(meanchange < 0.0){
				meanchange = 0.0;
			}
			// Signature: pnorm(double x, double mu, double sigma, int lower_tail, int log_p)
			*p_2 = R::pnorm(meanchange,red_mu,red_se,true,false);
		}
	}
	
	/*  This is now done in the wrapper as no other C function calls this with delta=1	
	if(delta > 0){
		PutRNGstate();
	}
	*/
	
	// printf("%f, %f\n", *p_1, *p_2);
	
	// Return an indication if one or other delta method failed:
	return (int) dfailed;
	
}

/*******************************************************************************/
/*   Underlying C++ function for conjugate Beta method (not exported to R)     */
/*******************************************************************************/

void conjbeta_ci(int N1, long long sum1, double K1, int N2, long long sum2, double K2, double mean_ratio, int iters, double *conjugate_priors, double tail, double *ci_l, double *ci_u){
	
	std::vector<double> ratio(iters, 0);

//	Theoretically faster but makes no difference in practice:
//	std::vector<double> ratio;
//	ratio.reserve(iters);
	
	double alpha1 = (double) N1 * K1 + conjugate_priors[0];
	double beta1 = sum1 + conjugate_priors[1];
	double alpha2 = (double) N2 * K2 + conjugate_priors[0];
	double beta2 = sum2 + conjugate_priors[1];
	
	for(int i=0; i<iters; i++){
		double p1 = R::rbeta(alpha1, beta1);
		double m1 = K1 / p1 - K1;
		double p2 = R::rbeta(alpha2, beta2);
		double m2 = K2 / p2 - K2;

		ratio[i] = 1.0 - (m2 / (m1 * mean_ratio));
//		ratio.push_back(1.0 - m2/m1);
	}
	int in = (int) ((double) iters * tail);
	std::nth_element(ratio.begin(), ratio.begin() + in, ratio.end());
	*ci_l = ratio[in];

	in = (int) ((double) iters * tail);
	std::nth_element(ratio.begin(), ratio.begin() + in, ratio.end(), std::greater<double>());
	*ci_u = ratio[in];

}

/**********************************************************************/
/*   Underlying C++ function for WAAVP paired (not exported to R)     */
/**********************************************************************/

void waavp_p_ci(double mu1, double mu2, double var1, double var2, double cov12, int N, double tail, double *ci_l, double *ci_u){
	
	// Method B of Lyndal-Murphy, M., Swain, a J., & Pepper, P. M. (2014). Methods to determine resistance to anthelmintics when continuing larval development occurs. Veterinary Parasitology, 199(3–4), 191–200. https://doi.org/10.1016/j.vetpar.2013.11.002

	int df = N -1;
	// Signature of qt is:  double  qt(double, double, int, int);
	double tval = R::qt(1.0 - tail, (double) df, 1, 0);
	
	double Nd = (double) N;
	double varred = var1 / (Nd * mu1 * mu1) + var2 / (Nd * mu2 * mu2) - 2.0 * cov12 / (Nd * mu1 * mu2);
	
	*ci_u = 1 - (mu2 / mu1 * std::exp(-tval * std::sqrt(varred)));
	*ci_l = 1 - (mu2 / mu1 * std::exp(tval * std::sqrt(varred)));

}

/**********************************************************************/
/*   Underlying C++ function for WAAVP unpaired (not exported to R)   */
/**********************************************************************/

void waavp_u_ci(double mu1, double mu2, double var1, double var2, int N1, int N2, double tail, double *ci_l, double *ci_u){
	
	// Method A of Lyndal-Murphy, M., Swain, a J., & Pepper, P. M. (2014). Methods to determine resistance to anthelmintics when continuing larval development occurs. Veterinary Parasitology, 199(3–4), 191–200. https://doi.org/10.1016/j.vetpar.2013.11.002
	
	int df = N1 + N2 - 2;
	// Signature of qt is:  double  qt(double, double, int, int);
	double tval = R::qt(1.0 - tail, (double) df, 1, 0);
	
	double N1d = (double) N1;
	double N2d = (double) N2;
	double varred = var1 / (N1d * mu1 * mu1) + var2 / (N2d * mu2 * mu2);
	
	*ci_u = 1 - (mu2 / mu1 * std::exp(-tval * std::sqrt(varred)));
	*ci_l = 1 - (mu2 / mu1 * std::exp(tval * std::sqrt(varred)));

}

/***********************************************************************/
/*   Underlying C++ function for Dobson method (not exported to R)     */
/***********************************************************************/

void dobson_ci(long long sum1, long long sum2, double lci_c, double uci_c, double *dobson_priors, double *ci_l, double *ci_u){
	
	double postsumd = (double) sum2;
	double sumdiffd = (double) (sum1 - sum2);
	*ci_l = 1.0 - R::qbeta(lci_c, postsumd+dobson_priors[0], sumdiffd+dobson_priors[1], 0, 0);
	*ci_u = 1.0 - R::qbeta(uci_c, postsumd+dobson_priors[0], sumdiffd+dobson_priors[1], 0, 0);
	
}

/*******************************************************************************/
/*   Underlying C++ function for MLE unpaired (not exported to R)              */
/*******************************************************************************/

void mle_u_ci(double mu1, double mu2, double var1, double var2, int N1, int N2, double tail, double *ci_l, double *ci_u){
	
	// If N1 and N2 are unequal use the smallest:
	double Nd = N1 < N2 ? (double) N1 : (double) N2;

	double delta = std::sqrt((std::pow(mu2, 2) / std::pow(mu1, 4) * var1 + std::pow(mu1, -2) * mu2) / Nd);
	double r = mu2 / mu1;
	
	// Signature of qt is:  double  qt(double, double, int, int);
	// double critval = qt(0.975, Nd-1.0, 1, 0);
	// But it is MLE so normal is appropriate:
	// Signature of qnorm is:  qnorm(double, double, double, int, int)
	double critval = R::qnorm((double) 1.0 - tail, 0.0, 1.0, 1, 0);
	// double critval = 1.96;  

	*ci_l = 1.0 - (r + critval * delta);
	*ci_u = 1.0 - (r - critval * delta);
	
	// Detect if the var1 or var2 are 0:
	if(var1 <= 0.0 || var2 <= 0.0){
		*ci_l = NA_REAL;
		*ci_u = NA_REAL;
	}
	
}

/*******************************************************************************/
/*   Underlying C++ function for MLE paired (not exported to R)                */
/*******************************************************************************/

void mle_p_ci(double mu1, double mu2, double var1, double var2, double cov12, int N, double tail, double *ci_l, double *ci_u){
	
//	double delta = std::sqrt((std::pow(mu2, 2) / std::pow(mu1, 4) * var1 + std::pow(mu1, -2) * mu2 + 2.0 * mu2 / std::pow(mu1, 3) * cov12) / (double) N);
//	double delta = std::sqrt((std::pow(mu2, 2) / std::pow(mu1, 4) * var1 + std::pow(mu1, -2) * mu2 - 2.0 * mu2 / std::pow(mu1, 3) * cov12) / (double) N);
//	double delta = std::sqrt((std::pow(mu2, 2) / std::pow(mu1, 4) * var1 + std::pow(mu1, -2) * mu2) / (double) N);
	
	// Adjust the variance to the non-correlated variance:
	double correlation = cov12 / (std::sqrt(var1) * std::sqrt(var2));
	var1 *= (1.0 - correlation);
	var2 *= (1.0 - correlation);	
	
//	double delta = std::sqrt(((std::pow(mu2, 2) / std::pow(mu1, 4) * var1 + std::pow(mu1, -2) * mu2) - (2.0 * mu2 / std::pow(mu1, 3) * cov12)) / (double) N);
	
	double delta = std::sqrt((std::pow(mu2, 2) / std::pow(mu1, 4) * var1 + std::pow(mu1, -2) * mu2) / double(N));
	double r = mu2 / mu1;
	
	// Signature of qt is:  double  qt(double, double, int, int);
	// double critval = qt(0.975, Nd-1.0, 1, 0);
	// But it is MLE so normal is appropriate:
	// Signature of qnorm is:  qnorm(double, double, double, int, int)
	double critval = R::qnorm((double) 1.0 - tail, 0.0, 1.0, 1, 0);
	// double critval = 1.96;  
	
	*ci_l = 1.0 - (r + critval * delta);
	*ci_u = 1.0 - (r - critval * delta);
	
	// Detect if the var1 or var2 are 0:
	if(var1 <= 0.0 || var2 <= 0.0){
		*ci_l = NA_REAL;
		*ci_u = NA_REAL;
	}
	
}

/*******************************************************************************/
/*   Underlying C++ function for Levecke unpaired (not exported to R)          */
/*******************************************************************************/

void levecke_u_ci(double mu1, double mu2, double var1, double var2, int N1, int N2, double tail, double *ci_l, double *ci_u){
	
	double varred = std::pow(mu2/mu1, 2) * (var1/std::pow(mu1,2) + var2/std::pow(mu2,2));
	
	// If N1 and N2 are unequal use the smallest:
	double Nd = N1 < N2 ? (double) N1 : (double) N2;
	double r = mu2/mu1;
	double shape = (std::pow(r, 2) * Nd) / varred;
	double scale = varred / (r * Nd);
	
	// Signature of qgamma is:  qgamma(double p, double alpha, double scale, int lower_tail, int log_p)
	*ci_l = 1.0 - R::qgamma(1.0-tail, shape, scale, 1, 0);
	*ci_u = 1.0 - R::qgamma(tail, shape, scale, 1, 0);

}

/*******************************************************************************/
/*   Underlying C++ function for Levecke paired (not exported to R)            */
/*******************************************************************************/

//  From Levecke et al:  Assessment of Anthelmintic Efficacy of Mebendazole inSchool Children in Six Countries Where Soil-TransmittedHelminths Are Endemic

void levecke_p_ci(double mu1, double mu2, double var1, double var2, double cov12, int N, double tail, double *ci_l, double *ci_u){
	
//	double cor12 = cov12 / std::sqrt(var1 * var2);
//	double varred = std::pow(mu2/mu1, 2) * (var1/std::pow(mu1,2) + var2/std::pow(mu2,2) - 2.0 * cor12 * std::sqrt(var1 * var2) / (mu1 * mu2));	
	// Equivalent:
	double varred = std::pow(mu2/mu1, 2) * (var1/std::pow(mu1,2) + var2/std::pow(mu2,2) - 2.0 * cov12 / (mu1 * mu2));	

	double Nd = (double) N;
	double r = mu2/mu1;
	double shape = (std::pow(r, 2) * Nd) / varred;
	double scale = varred / (r * Nd);

	// Signature of qgamma is:  qgamma(double p, double alpha, double scale, int lower_tail, int log_p)
	*ci_l = 1.0 - R::qgamma(1.0-tail, shape, scale, 1, 0);
	*ci_u = 1.0 - R::qgamma(tail, shape, scale, 1, 0);

}

