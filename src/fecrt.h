/*
This file is part of the bayescount package:  https://cran.r-project.org/package=bayescount
Copyright of Matthew James Denwood, licensed GPL>=2
The code is for analysis and power calculations for FECRT studies
General note: long long is used for all sum1 and sum2 in case int is too short on
some systems - for this reason CXX_STD = CXX11 is specified in the Makevars file
[search https://cran.r-project.org/doc/manuals/r-release/R-exts.html for "long long"]
Last updated by MJD June 2018
*/

/*
double pbnb_lower(long long q, double bnb_k, double bnb_alpha, double bnb_beta, bool inclusive);

double pbnb_upper(long long q, double bnb_k, double bnb_alpha, double bnb_beta, bool inclusive);
*/

double pbnb_lower(long long q, double bnb_k, double bnb_alpha, double bnb_beta);

double pbnb_upper(long long q, double bnb_k, double bnb_alpha, double bnb_beta);

int bnb_pval(long long sum1, int N1, double K1, double mu1, double var1, long long sum2, int N2, double K2, double mu2, double var2, double cov12, double mean_ratio, double H0_1, double H0_2, double *conjugate_priors, int delta, int beta_iters, int approx, double *p_1, double *p_2);

void conjbeta_ci(int N1, long long sum1, double K1, int N2, long long sum2, double K2, double mean_ratio, int iters, double *conjugate_priors, double tail, double *ci_l, double *ci_u);

void waavp_p_ci(double mu1, double mu2, double var1, double var2, double cov12, int N, double tail, double *ci_l, double *ci_u);

void waavp_u_ci(double mu1, double mu2, double var1, double var2, int N1, int N2, double tail, double *ci_l, double *ci_u);

void dobson_ci(long long sum1, long long sum2, double lci_c, double uci_c, double *dobson_priors, double *ci_l, double *ci_u);

void mle_u_ci(double mu1, double mu2, double var1, double var2, int N1, int N2, double tail, double *ci_l, double *ci_u);

void mle_p_ci(double mu1, double mu2, double var1, double var2, double cov12, int N, double tail, double *ci_l, double *ci_u);

void levecke_u_ci(double mu1, double mu2, double var1, double var2, int N1, int N2, double tail, double *ci_l, double *ci_u);

void levecke_p_ci(double mu1, double mu2, double var1, double var2, double cov12, int N, double tail, double *ci_l, double *ci_u);
