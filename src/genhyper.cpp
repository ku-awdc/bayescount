/*
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rcpp.h>

#include "dists.h"

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
*/