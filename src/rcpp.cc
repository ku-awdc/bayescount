//#include "data_sim.h"

//#include "Rcpp.h"
// This is required for bivariate gamma and precludes Rcpp.h:
#include <RcppDist.h>
#include "distribution.h"

#include <algorithm>  // std::random_shuffle
#include <cmath>  // std::sqrt

#include "fecrt.h"


double pbnb(int q, double k, double alpha, double beta, bool lower, bool inclusive){

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


Rcpp::NumericVector estimate_k(const double mean_1, const double var_1, const double mean_2, const double var_2, const double cov_12, const bool paired){

	Rcpp::NumericVector retval(2);

	// Remove Poisson variation:
	double var_eff_1 = var_1 - mean_1;
	double var_eff_2 = var_2 - mean_2;

	// If using the paired model and both adjusted variances are positive:
	double correlation = 0.0;
	if(paired && var_eff_1 > 0.0 && var_eff_2 > 0.0){
		// Assume that none of the covariance is due to Poisson variation:
		// correlation = cov_12 / (std::sqrt(var_eff_1) * std::sqrt(var_eff_2));
	}
	// Empirically this seems to give the best results:
	if(paired && var_1 > 0.0 && var_2 > 0.0){
		// Assume that none of the covariance is due to Poisson variation:
		correlation = cov_12 / (std::sqrt(var_1) * std::sqrt(var_2));
	}
	// TODO: correlation can come out negative here - make min 0.0??
	// TODO: ks are over-estimated when covariance is huge

	if(correlation < 0.0){
		//Rcpp::warning("Negative correlation generated");
	}
	if(correlation >= 1.0){
		Rcpp::warning("correlation >= 1.0 corrected");
		correlation = 0.99;
	}


	// Can happen if the control variance is less than the control mean:
	if( var_eff_1 <= 0.0 ){
		retval[0] = (mean_1 * mean_1) / (var_1 * (1.0 - correlation));
	}else{
		retval[0] = (mean_1 * mean_1) / (var_eff_1 * (1.0 - correlation));
	}

	// Can happen if the treatment variance is either 0 or less than the control mean:
	if( var_eff_2 <= 0.0 ){
		retval[1] = retval[0];
	}else{
		retval[1] = (mean_2 * mean_2) / (var_eff_2 * (1.0 - correlation));
	}

	return retval;

}


// Forward declaration:
double find_theta(const Rcpp::IntegerVector data, const double mu, const double ax, const double bx, const double tol);

Rcpp::NumericVector estimate_k_ml(const Rcpp::IntegerVector data_1, const double mean_1, const double var_1, const Rcpp::IntegerVector data_2, const double mean_2, const double var_2, const double cov_12, const bool paired, const double k_ratio = 1.0, const double min=0.001, const double max=20.0, const double tol=0.01){

	Rcpp::NumericVector retval(5);

	// If min and max are negative then they are relative to the Poisson variance
	// TODO: do we need a min or can it always be 0???
	double min_1 = min;
	double min_2 = min;
	double max_1 = max;
	double max_2 = max;

	if(min < 0.0){
		Rcpp::stop("Not implemented!");
	}
	if(max < 0.0){
		// max (p) is expressed as the minimum proportion of the total variance due to the Gamma part:
		// solve p = (m^2 / k) / ((m^2 / k) + m) for k
		// k = m * (1/p -1)
		//double m = std::max(mean_1, mean_2);
		double p = -max;
		max_1 = mean_1 * (1.0/p -1.0);
		max_2 = mean_2 * (1.0/p -1.0);

		// Correct but not needed:
		// double pois_cv = 1.0 / std::sqrt(std::max(mean_1, mean_2));
	}

	// Get MLE estimates of two k values:
	double k1 = find_theta(data_1, mean_1, min_1, max_1, tol);
	double k2 = find_theta(data_2, mean_2, min_2, max_2, tol);

	retval[3] = k1;
	retval[4] = k2;

	// TODO: allow negative correlation??  Need to calculate and report it at least
//	if(paired && var_1 > 0.0 && var_2 > 0.0 && cov_12 > 0.0){
	if(paired && var_1 > 0.0 && var_2 > 0.0){
		double correlation = 0.0;

		// Assuming that none of the covariance is due to Poisson variation:
		// correlation = cov_12 / (std::sqrt(var_1 - mean_1) * std::sqrt(var_2 - mean_2));

		// Empirically this seems to give better results:
		correlation = cov_12 / (std::sqrt(var_1) * std::sqrt(var_2));
		// Major advantage is that it is limited to -1,1

		// New approach - use estimated k to derive gamma variance:
//		correlation = cov_12 / (std::sqrt(std::min(var_1, mean_1*mean_1 / k1)) * std::sqrt(std::min(var_2, mean_2*mean_2 / k2)));
//		correlation = cov_12 / ((mean_1 / std::sqrt(k1)) * (mean_2 / std::sqrt(k2)));
		// This doesn't work because there is a correlation >1 too often

		retval[2] = correlation;

		if(correlation >= 1.0){
			Rcpp::warning("correlation >= 1.0");
			correlation = 0.99;
		}
		if(correlation > 0.9){
			correlation = 0.9;
		}

		// Adjust k only if positive correlation:
		if(correlation > 0.0){
			k1 /= (1.0 - correlation);
			k2 /= (1.0 - correlation);
		}
	}

	// TODO: might be that var_2 > 0 but var_1 = 0??
	if(var_2 <= 0.0){
		k2 = k1 * k_ratio;
		retval[4] = NA_REAL;
	}

	if(!R_finite(k2) || k2 <= 0.0){
		Rcpp::Rcout << k2 << " - " << mean_2 << " - " << var_2 << " - " << cov_12 << "\n";
	}

	retval[0] = k1;
	retval[1] = k2;

	return retval;

}



Rcpp::String get_type_ci(const double eff, const double lci, const double uci, const double ta, const double ti){

	Rcpp::String tp;

	if(uci < ta){
		tp = "1ab";
	}else if(uci < ti && eff < ta){
		tp = "1ab";
	}else if(uci < ti && lci < ta){
		tp = "1c";
	}else if(lci < ta && uci >= ti && eff < ta){
		tp = "2a";
	}else if(lci < ta && uci >= ti && eff >= ta && eff < ti){
		tp = "2b";
	}else if(lci < ta && uci >= ti && eff >= ti){
		tp = "2c";
	}else if(lci >= ta && uci < ti){
		tp = "3";
	}else if(lci >= ta && uci >= ti && eff < ti){
		tp = "4a";
	}else if(lci >= ta && uci >= ti && eff >= ti){
		tp = "4bc";
	}else if(lci >= ti){
		tp = "4bc";
	}else{
//		Rcpp::Rcout << eff << " - " << lci << " - " << uci << " - " << ta << " - " << ti << "\n";
		tp = "error";
//		Rcpp::stop("Error assigning ci type");
	}

	return tp;
}

Rcpp::String get_type_pv(const double eff, const double p1, const double p2, const double tail, const double ta, const double ti){

	Rcpp::String tp;

	if(p2 <= tail && p1 > tail && eff < ta){
		tp = "1ab";
	}else if(p2 <= tail && p1 > tail && eff >= ta){
		tp = "1c";
	}else if(p2 > tail && p1 > tail && eff < ta){
		tp = "2a";
	}else if(p2 > tail && p1 > tail && eff >= ta && eff < ti){
		tp = "2b";
	}else if(p2 > tail && p1 > tail && eff >= ti){
		tp = "2c";
	}else if(p2 <= tail && p1 <= tail){
		tp = "3";
	}else if(p2 > tail && p1 <= tail && eff < ti){
		tp = "4a";
	}else if(p2 > tail && p1 <= tail && eff >= ti){
		tp = "4bc";
	}else{
		tp = "error";
//		Rcpp::stop("Error assigning pv type");
	}

	return tp;
}


Rcpp::DataFrame efficacy_analysis(Rcpp::IntegerVector data_1, Rcpp::IntegerVector data_2, bool paired, double mean_ratio, double H0_I, double H0_A, Rcpp::NumericVector conjugate_priors, int delta, int beta_iters, bool useml, int approx, double tail, Rcpp::NumericVector dobson_cl, Rcpp::NumericVector dobson_priors, Rcpp::NumericVector known_ks)
{
	using namespace Rcpp;

	const int nrows = paired ? 7L : 6L;
	Rcpp::DoubleVector td( nrows, NA_REAL );
	Rcpp::StringVector ts( nrows );
	Rcpp::StringVector output_Method = clone(ts);
	Rcpp::DoubleVector output_LCI = clone(td);
	Rcpp::DoubleVector output_UCI = clone(td);
	Rcpp::DoubleVector output_pI = clone(td);
	Rcpp::DoubleVector output_pA = clone(td);
	Rcpp::StringVector output_Classification = clone(ts);

	// Note:  Rcpp allows conversion from long long but not to double array directly
	double conjugate_priors_db[2] = { conjugate_priors[0], conjugate_priors[1] };
	double dobson_priors_db[2] = { dobson_priors[0], dobson_priors[1] };

	const int N_1 = data_1.size();
	const int N_2 = data_1.size();

	long long sum_1 = 0L;
	for(int i=0; i<N_1; ++i){
		sum_1 += data_1[i];
	}
	const double mean_1 = Rcpp::mean(data_1);
	const double var_1 = Rcpp::var(data_1);

	long long sum_2 = 0L;
	for(int i=0; i<N_2; ++i){
		sum_2 += data_2[i];
	}
	const double mean_2 = Rcpp::mean(data_2);
	const double var_2 = Rcpp::var(data_2);

	const double obsred = mean_2 / mean_1;

	double cov = 0.0;
	if(paired)
	{
		if(N_1 != N_2)
		{
			stop("Unequal length of pre- and post-treatment data");
		}
	    for(int i=0; i<N_1; i++){
			cov += (double(data_1[i]) - mean_1) * (double(data_2[i]) - mean_2);
		}
	    cov = cov / double(N_1 - 1L);
		if(cov > 0.9)
			{
				cov = 0.9;
			}
	}

	Rcpp::NumericVector ks;
	if(useml){
		ks = estimate_k_ml(data_1, mean_1, var_1, data_2, mean_2, var_2, cov, true);
	}else{
		ks = estimate_k(mean_1, var_1, mean_2, var_2, cov, true);
	}

	// Test:
	if(!R_finite(ks[0])){
		stop("Non-finite k1 generated");
	}
	if(ks[0] <= 0.0){
		stop("k1 below 0.0 generated");
	}
	if(!R_finite(ks[1])){
		stop("Non-finite k2 generated");
	}
	if(ks[1] <= 0.0){
		stop("k2 below 0.0 generated");
	}

	int row = 0L;

	output_Method[row] = "BNB";
	bnb_pval(sum_1, N_1, ks[0L], mean_1, var_1, sum_2, N_2, ks[1L], mean_2, var_2, cov, mean_ratio, H0_A, H0_I, conjugate_priors_db, delta, beta_iters, approx, &output_pA[row], &output_pI[row]);
	output_Classification[row] = get_type_pv(1.0-obsred, output_pA[row], output_pI[row], tail, H0_A, H0_I);
	row++;

	output_Method[row] = "BNB_FixK2";
	bnb_pval(sum_1, N_1, ks[0L], mean_1, var_1, sum_2, N_2, ks[0L], mean_2, var_2, cov, mean_ratio, H0_A, H0_I, conjugate_priors_db, delta, beta_iters, approx, &output_pA[row], &output_pI[row]);
	output_Classification[row] = get_type_pv(1.0-obsred, output_pA[row], output_pI[row], tail, H0_A, H0_I);
	row++;

	output_Method[row] = "BNB_KnownKs";
	// TODO: known_ks will be ignored when using approximation - refactor code!
	bnb_pval(sum_1, N_1, known_ks[0L], mean_1, var_1, sum_2, N_2, known_ks[1L], mean_2, var_2, cov, mean_ratio, H0_A, H0_I, conjugate_priors_db, delta, beta_iters, approx, &output_pA[row], &output_pI[row]);
	output_Classification[row] = get_type_pv(1.0-obsred, output_pA[row], output_pI[row], tail, H0_A, H0_I);
	row++;

	if(paired)
	{
		output_Method[row] = "Gamma";
		levecke_p_ci(mean_1, mean_2, var_1, var_2, cov, N_1, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

		output_Method[row] = "WAAVP";
		waavp_p_ci(mean_1, mean_2, var_1, var_2, cov, N_1, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

		output_Method[row] = "Asymptotic";
		mle_p_ci(mean_1, mean_2, var_1, var_2, cov, N_1, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

		output_Method[row] = "Binomial";
		dobson_ci(sum_1, sum_2, dobson_cl[0], dobson_cl[1], dobson_priors_db, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2>sum_1) ? "<0%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

	}
	else
	{
		output_Method[row] = "Gamma";
		levecke_u_ci(mean_1, mean_2, var_1, var_2, N_1, N_2, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

		output_Method[row] = "WAAVP";
		waavp_u_ci(mean_1, mean_2, var_1, var_2, N_1, N_2, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

		output_Method[row] = "Asymptotic";
		mle_u_ci(mean_1, mean_2, var_1, var_2, N_1, N_2, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

	}

	DataFrame output = DataFrame::create(Named("Method") = output_Method,
		Named("LCI") = output_LCI, Named("UCI") = output_UCI,
		Named("pI") = output_pI, Named("pA") = output_pA, Named("Typology") = output_Classification);

	return output;
}

Rcpp::DataFrame typology_analysis(const int sum_1, const int sum_2, const int N_1, const int N_2, const double k_1, const double k_2, const double cor, const bool paired, const double mean_ratio, const double H0_I, const double H0_A, const Rcpp::NumericVector conjugate_priors, const int delta, const int beta_iters, const int approx, const double tail, const Rcpp::NumericVector dobson_cl, const Rcpp::NumericVector dobson_priors)
{
	using namespace Rcpp;

	const int nrows = paired ? 5L : 4L;
	Rcpp::DoubleVector td( nrows, NA_REAL );
	Rcpp::StringVector ts( nrows );
	Rcpp::StringVector output_Method = clone(ts);
	Rcpp::DoubleVector output_LCI = clone(td);
	Rcpp::DoubleVector output_UCI = clone(td);
	Rcpp::DoubleVector output_pI = clone(td);
	Rcpp::DoubleVector output_pA = clone(td);
	Rcpp::StringVector output_Classification = clone(ts);

	// Note:  Rcpp allows conversion from long long but not to double array directly
	double conjugate_priors_db[2] = { conjugate_priors[0], conjugate_priors[1] };
	double dobson_priors_db[2] = { dobson_priors[0], dobson_priors[1] };

	const double mean_1 = static_cast<double>(sum_1) / static_cast<double>(N_1);
	const double var_1 = mean_1 + (mean_1*mean_1) / k_1;
	const double mean_2 = static_cast<double>(sum_2) / static_cast<double>(N_2);
	const double var_2 = mean_2 + (mean_2*mean_2) / k_2;

	const double obsred = mean_2 / mean_1;

	double cov = 0.0;
	if(paired)
	{
	    cov = cor * std::sqrt(var_1) * std::sqrt(var_2);
	}
	Rcpp::NumericVector ks = estimate_k(mean_1, var_1, mean_2, var_2, cov, true);

	int row = 0L;

	output_Method[row] = "BNB";
	bnb_pval(sum_1, N_1, ks[0L], mean_1, var_1, sum_2, N_2, ks[1L], mean_2, var_2, cov, mean_ratio, H0_A, H0_I, conjugate_priors_db, delta, beta_iters, approx, &output_pA[row], &output_pI[row]);
	output_Classification[row] = get_type_pv(1.0-obsred, output_pA[row], output_pI[row], tail, H0_A, H0_I);
	row++;

	if(paired)
	{
		output_Method[row] = "Gamma";
		levecke_p_ci(mean_1, mean_2, var_1, var_2, cov, N_1, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

		output_Method[row] = "WAAVP";
		waavp_p_ci(mean_1, mean_2, var_1, var_2, cov, N_1, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

		output_Method[row] = "Asymptotic";
		mle_p_ci(mean_1, mean_2, var_1, var_2, cov, N_1, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

		output_Method[row] = "Binomial";
		dobson_ci(sum_1, sum_2, dobson_cl[0], dobson_cl[1], dobson_priors_db, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2>sum_1) ? "<0%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

	}
	else
	{
		output_Method[row] = "Gamma";
		levecke_u_ci(mean_1, mean_2, var_1, var_2, N_1, N_2, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

		output_Method[row] = "WAAVP";
		waavp_u_ci(mean_1, mean_2, var_1, var_2, N_1, N_2, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

		output_Method[row] = "Asymptotic";
		mle_u_ci(mean_1, mean_2, var_1, var_2, N_1, N_2, tail, &output_LCI[row], &output_UCI[row]);
		output_Classification[row] = (sum_2==0) ? "100%red" : get_type_ci(1.0-obsred, output_LCI[row], output_UCI[row], H0_A, H0_I);
		row++;

	}

	DataFrame output = DataFrame::create(Named("Method") = output_Method,
		Named("LCI") = output_LCI, Named("UCI") = output_UCI,
		Named("pI") = output_pI, Named("pA") = output_pA, Named("Typology") = output_Classification);

	return output;
}


Rcpp::DataFrame efficacy_frequencies_unpaired(const int iters, const Rcpp::NumericVector red, const int N_ctl, const int N_tx, const double mu, const double k_tx, const double k_ctl, const Rcpp::NumericMatrix thresholds, const Rcpp::NumericVector conjugate_priors, const int delta, const int beta_iters, const int approx, const double tail, const bool useml)
{

	// Note: k_pre and k_post are the correlated k - I think - check!!!

	// Check that the thresholds have 2 columns and >0 rows:
	if(thresholds.ncol() != 2){
		Rcpp::stop("There must be exactly 2 columns");
	}
	if(thresholds.nrow() < 1){
		Rcpp::stop("There must be 1 or more rows");
	}
	for(int i=0; i<thresholds.nrow(); i++){
		if(thresholds(i,1) < thresholds(i,0)){
			Rcpp::stop("Second thresholds must be higher than the first");
		}
		if(thresholds(i,0) < 0 || thresholds(i,0) > 1){
			Rcpp::stop("Thresholds must be between 0 and 1");
		}
		if(thresholds(i,1) < 0 || thresholds(i,1) > 1){
			Rcpp::stop("Thresholds must be between 0 and 1");
		}
	}
	// Check that N has 2 columns and >0 rows:


	// Use the expand.grid R function:
	Rcpp::Function expGrid("expand.grid");

	Rcpp::StringVector methods = { "BNB", "Gamma", "WAAVP", "Asymptotic" };
	Rcpp::IntegerVector tseq = Rcpp::seq(1L, thresholds.nrow());
	Rcpp::DoubleVector td = { NA_REAL };
	Rcpp::StringVector ts = { NA_STRING };
//	Rcpp::IntegerVector ti = { NA_INTEGER };
	Rcpp::IntegerVector iteration = Rcpp::seq(1L, iters);

	Rcpp::DataFrame output = expGrid(methods, red, tseq, td, td, iteration, td, td, td, ts, Rcpp::Named("stringsAsFactors")=false);
	output.names() = Rcpp::StringVector::create("Method", "Reduction", "ThresholdIndex", "ThresholdLower", "ThresholdUpper", "Iteration", "ObsReduction", "Stat1", "Stat2", "Typology");

	// Get references to list items:
	Rcpp::StringVector output_Method = output[0L];
	Rcpp::DoubleVector output_Reduction = output[1L];
	Rcpp::IntegerVector output_ThresholdIndex = output[2L];
	Rcpp::DoubleVector output_ThresholdLower = output[3L];
	Rcpp::DoubleVector output_ThresholdUpper = output[4L];
	Rcpp::IntegerVector output_Iteration = output[5L];
	Rcpp::DoubleVector output_ObsReduction = output[6L];
	Rcpp::DoubleVector output_Stat1 = output[7L];
	Rcpp::DoubleVector output_Stat2 = output[8L];
	Rcpp::StringVector output_Classification = output[9L];

	// Note:  Rcpp allows conversion from long long but not to double array directly
	double conjugate_priors_db[2] = { conjugate_priors[0], conjugate_priors[1] };

	// To control the row index to write to:
	int row = 0L;

	// TODO: tidy up check and error:
	if((double(N_ctl)*double(mu)) > (0.1 * INT_MAX) || (double(N_ctl)*double(mu)) > (0.1 * INT_MAX)){
		Rcpp::stop("Possible int overflow");
	}

	// First loop over control datasets to simulate:
	for(int i=0; i<iters; i++){
		Rcpp::IntegerVector ctl_data;
		long long sum_ctl = 0L;
		do {
			sum_ctl = 0L;
			ctl_data = Rcpp::rnbinom_mu(N_ctl, k_ctl, mu);
			for(int i=0; i<N_tx; i++){
				sum_ctl += ctl_data[i];
			}
		} while (sum_ctl == 0L);

		double mean_ctl = Rcpp::mean(ctl_data);
		double var_ctl = Rcpp::var(ctl_data);

		// Then loop over reductions:
		for(int r=0; r<red.length(); r++){
			Rcpp::IntegerVector tx_data;
			tx_data = Rcpp::rnbinom_mu(N_tx, k_tx, mu*red[r]);

			// Avoid possible segfaults with int vs long long:
			long long sum_tx = 0L;
			for(int i=0; i<N_tx; i++){
				sum_tx += tx_data[i];
			}
			double mean_tx = Rcpp::mean(tx_data);
			double var_tx = Rcpp::var(tx_data);
			double obsred = mean_tx / mean_ctl;

			Rcpp::NumericVector ks;
			if(useml){
				ks = estimate_k_ml(ctl_data, mean_ctl, var_ctl, tx_data, mean_tx, var_tx, 0.0, false);
			}else{
				ks = estimate_k(mean_ctl, var_ctl, mean_tx, var_tx, 0.0, false);
			}
			// Test:
			if(!R_finite(ks[0])){
				Rcpp::stop("Non-finite k1 generated");
			}
			if(ks[0] <= 0.0){
				Rcpp::stop("k1 below 0.0 generated");
			}
			if(!R_finite(ks[1])){
				Rcpp::stop("Non-finite k2 generated");
			}
			if(ks[1] <= 0.0){
				Rcpp::stop("k2 below 0.0 generated");
			}

			double estk_ctl = ks[0];
			double estk_tx = ks[1];

			// For cheating:
			//estk_ctl = k_ctl;
			//estk_tx = k_tx;

			// Then loop over thresholds:
			for(int t=0; t<thresholds.nrow(); t++){

				// This stays the same for all rows within this threshold:
				int tind = output_ThresholdIndex[row] - 1L;
				double th1 = thresholds(tind,0L);
				double th2 = thresholds(tind,1L);

				// TODO: Method should not need to be set but do as a debug check?

				// Apply each test and save the results:
				output_Method[row] = "BNB";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				bnb_pval(sum_ctl, N_ctl, estk_ctl, mean_ctl, var_ctl, sum_tx, N_tx, estk_tx, mean_tx, var_tx, 0.0, 1.0, th1, th2, conjugate_priors_db, delta, beta_iters, approx, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = get_type_pv(1.0-obsred, output_Stat1[row], output_Stat2[row], tail, th1, th2);
				row++;

				output_Method[row] = "Gamma";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				levecke_u_ci(mean_ctl, mean_tx, var_ctl, var_tx, N_ctl, N_tx, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_tx==0) ? "100%red" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "WAAVP";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				waavp_u_ci(mean_ctl, mean_tx, var_ctl, var_tx, N_ctl, N_tx, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_tx==0) ? "100%red" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "Asymptotic";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				mle_u_ci(mean_ctl, mean_tx, var_ctl, var_tx, N_ctl, N_tx, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_tx==0) ? "100%red" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				Rcpp::checkUserInterrupt();

			}
		}
	}

	return output;
}


Rcpp::DataFrame efficacy_frequencies_paired(int iters, Rcpp::NumericVector red, int N, double mu, double k_pre, double k_post, double k_c, Rcpp::NumericMatrix thresholds, Rcpp::NumericVector conjugate_priors, int delta, int beta_iters, int approx, double tail, bool useml, Rcpp::NumericVector dobson_cl, Rcpp::NumericVector dobson_priors)
{

	// Note: k_pre and k_post are the correlated k - I think - check!!!

	// Check that the thresholds have 2 columns and >0 rows:
	if(thresholds.ncol() != 2){
		Rcpp::stop("There must be exactly 2 columns");
	}
	if(thresholds.nrow() < 1){
		Rcpp::stop("There must be 1 or more rows");
	}
	for(int i=0; i<thresholds.nrow(); i++){
		if(thresholds(i,1) < thresholds(i,0)){
			Rcpp::stop("Second thresholds must be higher than the first");
		}
		if(thresholds(i,0) < 0 || thresholds(i,0) > 1){
			Rcpp::stop("Thresholds must be between 0 and 1");
		}
		if(thresholds(i,1) < 0 || thresholds(i,1) > 1){
			Rcpp::stop("Thresholds must be between 0 and 1");
		}
	}

	// Check that k_c is greater than pre and post:
	if(k_c < k_pre || k_c < k_post){
		Rcpp::stop("k_c must be larger than k_pre and k_post");
	}

	// TODO: tidy up check and error:
	if((double(N)*double(mu)) > (0.1 * INT_MAX)){
		Rcpp::stop("Possible int overflow");
	}

	// Use the expand.grid R function:
	Rcpp::Function expGrid("expand.grid");

	Rcpp::StringVector methods = { "BNB", "Gamma", "WAAVP", "Asymptotic", "Binomial" };
	Rcpp::IntegerVector tseq = Rcpp::seq(1L, thresholds.nrow());
	Rcpp::DoubleVector td = { NA_REAL };
	Rcpp::StringVector ts = { NA_STRING };
//	Rcpp::IntegerVector ti = { NA_INTEGER };
	Rcpp::IntegerVector iteration = Rcpp::seq(1L, iters);

	Rcpp::DataFrame output = expGrid(methods, red, tseq, td, td, iteration, td, td, td, ts, Rcpp::Named("stringsAsFactors")=false);
	output.names() = Rcpp::StringVector::create("Method", "Reduction", "ThresholdIndex", "ThresholdLower", "ThresholdUpper", "Iteration", "ObsReduction", "Stat1", "Stat2", "Typology");

	// Get references to list items:
	Rcpp::StringVector output_Method = output[0L];
	Rcpp::DoubleVector output_Reduction = output[1L];
	Rcpp::IntegerVector output_ThresholdIndex = output[2L];
	Rcpp::DoubleVector output_ThresholdLower = output[3L];
	Rcpp::DoubleVector output_ThresholdUpper = output[4L];
	Rcpp::IntegerVector output_Iteration = output[5L];
	Rcpp::DoubleVector output_ObsReduction = output[6L];
	Rcpp::DoubleVector output_Stat1 = output[7L];
	Rcpp::DoubleVector output_Stat2 = output[8L];
	Rcpp::StringVector output_Classification = output[9L];

	// Note:  Rcpp allows conversion from long long but not to double array directly
	double conjugate_priors_db[2] = { conjugate_priors[0], conjugate_priors[1] };
	double dobson_priors_db[2] = { dobson_priors[0], dobson_priors[1] };

	// To control the row index to write to:
	int row = 0L;

	// Precalculate some things:
	double muadj_pre = mu / k_pre;
	double muadj_post = mu / k_post;
	double b_pre = k_c - k_pre;
	double b_post = k_c - k_post;

	// First loop over correlation and pre-treatment datasets to simulate:
	for(int i=0; i<iters; i++){
		Rcpp::NumericVector gammas;
		gammas = Rcpp::rgamma(N, k_c, 1.0);

		Rcpp::IntegerVector pre_data(N);
		long long sum_pre = 0L;
		do {
			sum_pre = 0L;
			for(int i=0; i<N; i++){
				pre_data[i] = R::rpois(R::rbeta(k_pre, b_pre) * gammas[i] * muadj_pre);
				sum_pre += pre_data[i];
			}
		} while (sum_pre == 0L);

		double mean_pre = Rcpp::mean(pre_data);
		double var_pre = Rcpp::var(pre_data);

		// Then loop over reductions:
		for(int r=0; r<red.length(); r++){
			long long sum_post = 0L;
			Rcpp::IntegerVector post_data(N);
			for(int i=0; i<N; i++){
				post_data[i] = R::rpois(R::rbeta(k_post, b_post) * gammas[i] * muadj_post * red[r]);
				sum_post += post_data[i];
			}

			double mean_post = Rcpp::mean(post_data);
			double var_post = Rcpp::var(post_data);
			double obsred = mean_post / mean_pre;

			double cov = 0.0;
		    for(int i=0; i<N; i++){
				cov += (double(pre_data[i]) - mean_pre) * (double(post_data[i]) - mean_post);
			}
		    cov = cov / double(N - 1L);

			Rcpp::NumericVector ks;
			if(useml){
				ks = estimate_k_ml(pre_data, mean_pre, var_pre, post_data, mean_post, var_post, cov, true);
			}else{
				ks = estimate_k(mean_pre, var_pre, mean_post, var_post, cov, true);
			}
			// Test:
			if(!R_finite(ks[0])){
				Rcpp::stop("Non-finite k1 generated");
			}
			if(ks[0] <= 0.0){
				Rcpp::stop("k1 below 0.0 generated");
			}
			if(!R_finite(ks[1])){
				Rcpp::stop("Non-finite k2 generated");
			}
			if(ks[1] <= 0.0){
				Rcpp::stop("k2 below 0.0 generated");
			}
			double estk_pre = ks[0];
			double estk_post = ks[1];

			// For cheating:
			//estk_pre = k_pre / (1.0 - (std::sqrt(k_pre) * std::sqrt(k_post) / k_c));
			//estk_post = k_post / (1.0 - (std::sqrt(k_pre) * std::sqrt(k_post) / k_c));


			// Then loop over thresholds:
			for(int t=0; t<thresholds.nrow(); t++){

				// This stays the same for all rows within this threshold:
				int tind = output_ThresholdIndex[row] - 1L;
				double th1 = thresholds(tind,0L);
				double th2 = thresholds(tind,1L);


				// TODO: Method should not need to be set but do as a debug check?

				// Apply each test and save the results:

				output_Method[row] = "BNB";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				bnb_pval(sum_pre, N, estk_pre, mean_pre, var_pre, sum_post, N, estk_post, mean_post, var_post, cov, 1.0, th1, th2, conjugate_priors_db, delta, beta_iters, approx, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = get_type_pv(1.0-obsred, output_Stat1[row], output_Stat2[row], tail, th1, th2);
				row++;

				output_Method[row] = "Gamma";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				levecke_p_ci(mean_pre, mean_post, var_pre, var_post, cov, N, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post==0) ? "100%red" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "WAAVP";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				waavp_p_ci(mean_pre, mean_post, var_pre, var_post, cov, N, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post==0) ? "100%red" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "Asymptotic";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				mle_p_ci(mean_pre, mean_post, var_pre, var_post, cov, N, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post==0) ? "100%red" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "Binomial";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				dobson_ci(sum_pre, sum_post, dobson_cl[0], dobson_cl[1], dobson_priors_db, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post>sum_pre) ? "<0%red" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				Rcpp::checkUserInterrupt();

			}
		}
	}

	return output;
}


// TODO: power_matrix_unpaired


Rcpp::DataFrame power_matrix_paired(Rcpp::NumericVector Ns, Rcpp::NumericVector margin, Rcpp::NumericVector mu, double target, double k_pre, double k_post, double k_c, int iters, Rcpp::NumericVector conjugate_priors, int delta, int beta_iters, int approx, double tail, int useml)
{

	// Note: k_pre and k_post are the correlated k - I think - check!!!

	// TODO: constification

	// Check that all nim are >0:
	for(int i=0; i<margin.length(); ++i)
	{
		if(margin[i] <= 0.0)
		{
			Rcpp::stop("The non-inferiority margin must be strictly positive");
		}
	}

	// Check that k_c is greater than pre and post:
	if(k_c < k_pre || k_c < k_post){
		Rcpp::stop("k_c must be larger than k_pre and k_post");
	}

	// N must be small to large:
	std::sort(Ns.begin(), Ns.end());
	const int maxN = Rcpp::max(Ns);
	if((double(maxN)*Rcpp::max(mu)) > (0.1 * INT_MAX)){
		Rcpp::stop("Possible int overflow");
	}

	// Use the expand.grid R function:
	Rcpp::Function expGrid("expand.grid");

	Rcpp::DoubleVector td = { NA_REAL };
	Rcpp::StringVector ts = { NA_STRING };
	Rcpp::LogicalVector tl = Rcpp::LogicalVector::create( NA_LOGICAL );
//	Rcpp::IntegerVector ti = { NA_INTEGER };
	Rcpp::IntegerVector iteration = Rcpp::seq(1L, iters);
	Rcpp::StringVector hs = { "Efficacy", "Resistance" };

	// The ordering is important here:  hypothesis must change fastest (use same dataset), then N (same dataset subselected to maxN), then margin/reduction (use same pre-tx dataset), then iteration (resimulate), then mu (recalculate muadj etc)
	Rcpp::DataFrame output = expGrid(hs, Ns, target, margin, td, iteration, mu, td, td, td, ts, tl, Rcpp::Named("stringsAsFactors")=false);
	output.names() = Rcpp::StringVector::create("Hypothesis", "N", "Target", "Margin", "Reduction", "Iteration", "Mean", "ObsReduction", "Stat1", "Stat2", "Typology", "RejectH0");

	// Get references to list items:
	Rcpp::StringVector output_Hypothesis = output[0L];
	Rcpp::IntegerVector output_Ns = output[1L];
	Rcpp::DoubleVector output_Target = output[2L];
	Rcpp::DoubleVector output_Margin = output[3L];
	Rcpp::DoubleVector output_Reduction = output[4L];
	Rcpp::IntegerVector output_Mean = output[5L];
	Rcpp::IntegerVector output_Iteration = output[6L];
	Rcpp::DoubleVector output_ObsReduction = output[7L];
	Rcpp::DoubleVector output_Stat1 = output[8L];
	Rcpp::DoubleVector output_Stat2 = output[9L];
	Rcpp::StringVector output_Typology = output[10L];
	Rcpp::LogicalVector output_RejectH0 = output[11L];

	// Note:  Rcpp allows conversion from long long but not to double array directly
	double conjugate_priors_db[2] = { conjugate_priors[0], conjugate_priors[1] };

	// Precalculate some things:
	const double b_pre = k_c - k_pre;
	const double b_post = k_c - k_post;
	const double target_red = 1.0 - target;

	// For use with useml==0L (will be overwritten for other useml):
	Rcpp::DoubleVector ks = { k_pre, k_post };

	// To control the row index to write to:
	int row = 0L;

	// First loop over mus:
	for(int u=0; u<mu.length(); ++u)
	{

		// Precalculate some things for this mu:
		const double muadj_pre = mu[u] / k_pre;
		const double muadj_post = mu[u] / k_post;

		// Then loop over iterations to simulate one pre-treatment dataset common to all margins and N:
		for(int it=0; it<iters; ++it){

			Rcpp::checkUserInterrupt();

			Rcpp::NumericVector gammas;
			gammas = Rcpp::rgamma(maxN, k_c, 1.0);

			Rcpp::IntegerVector pre_full(maxN);
			do {
				for(int i=0; i<maxN; i++){
					pre_full[i] = R::rpois(R::rbeta(k_pre, b_pre) * gammas[i] * muadj_pre);
				}
			} while (pre_full[0L] == 0L);  // Hack to make sure all possible sample sizes have sum > 0

			// Then create a single post-treatment dataset for efficacious reductions (fixed over margins):
			Rcpp::IntegerVector post_full_sus(maxN);
			for(int i=0; i<maxN; i++){
				post_full_sus[i] = R::rpois(R::rbeta(k_post, b_post) * gammas[i] * muadj_post * target_red);
			}
			// TODO: use accumulating algorithm for calculating mean/var/sum eficiently
			// TODO: loop over sample sizes here to calculate pre-treatment (and post-treatment efficacious) means and vars once
			// double mean_pre = Rcpp::mean(pre_data);
			// double var_pre = Rcpp::var(pre_data);

			// Then loop over margins (and therefore resistant reductions):
			for(int m=0; m<margin.length(); ++m){

				const double red = target_red + margin[m];

				// Then create the post-treatment data for reduced reductions:
				Rcpp::IntegerVector post_full_res(maxN);
				for(int i=0; i<maxN; i++){
					post_full_res[i] = R::rpois(R::rbeta(k_post, b_post) * gammas[i] * muadj_post * red);
				}
				// TODO: use accumulating algorithm for calculating mean/var/sum eficiently
				// TODO: loop over sample sizes here to calculate pre-treatment (and post-treatment efficacious) means and vars once
				// double mean_pre = Rcpp::mean(pre_data);
				// double var_pre = Rcpp::var(pre_data);

				// Then loop over sample sizes:
				for(int n=0; n<Ns.length(); ++n){

					/*
					std::random_shuffle(pre_full.begin(), pre_full.end());
					std::random_shuffle(post_full_sus.begin(), post_full_sus.end());
					std::random_shuffle(post_full_res.begin(), post_full_res.end());
					*/

					const int N = Ns[n];

					// Indices to use for this sample size:
					// Rcpp::IntegerVector idx = Rcpp::seq(0L, N-1L);
					// TODO: work out how to subset vectors rather than copying

					Rcpp::IntegerVector pre_data(N);
					Rcpp::IntegerVector post_data(N);
					std::copy(pre_full.begin(), pre_full.begin()+N, pre_data.begin());

					// TODO: use once-calculated values here (or is LLVM smart enough to do this for me?!?):
					double mean_pre = Rcpp::mean(pre_data);
					double var_pre = Rcpp::var(pre_data);

					// Then loop over hypotheses:
					for(int h=0L; h<2L; ++h)
					{
						if(h==0L){
							// Hypothesis for efficacy comes first:
							std::copy(post_full_sus.begin(), post_full_sus.begin()+N, post_data.begin());
						}else{
							// Hypothesis for resistance comes second:
							std::copy(post_full_res.begin(), post_full_res.begin()+N, post_data.begin());
						}

						double mean_post = Rcpp::mean(post_data);
						double var_post = Rcpp::var(post_data);
						double obsred = mean_post / mean_pre;

						double cov = 0.0;

						long long sum_pre = 0L;
						long long sum_post = 0L;
					    for(int i=0; i<N; i++){
							cov += (double(pre_data[i]) - mean_pre) * (double(post_data[i]) - mean_post);
							sum_pre += pre_data[i];
							sum_post += post_data[i];
						}
					    cov = cov / double(N - 1L);

						// useml is either 0 (use known k), 1 (use faster approximation to estimate) or 2 (use actual ml to estimate)
						if(useml==0L){
							// Note: do nothing here for useml==0L
						}else if(useml==1L){
							ks = estimate_k(mean_pre, var_pre, mean_post, var_post, cov, true);
						}else if(useml==2L){
							ks = estimate_k_ml(pre_data, mean_pre, var_pre, post_data, mean_post, var_post, cov, true);
						}else{
							Rcpp::stop("Unrecognised useml value");
						}
						double estk_pre = ks[0];
						double estk_post = ks[1];

						const double th1 = target - margin[m];
						const double th2 = target;

						bnb_pval(sum_pre, N, estk_pre, mean_pre, var_pre, sum_post, N, estk_post, mean_post, var_post, cov, 1.0, th1, th2, conjugate_priors_db, delta, beta_iters, approx, &output_Stat1[row], &output_Stat2[row]);
						output_Typology[row] = get_type_pv(1.0-obsred, output_Stat1[row], output_Stat2[row], tail, th1, th2);

						// Record the reduction etc:
						output_ObsReduction[row] = obsred;

						// Record if the relevant hypothesis has been rejected:
						if(h==0L){
							// Hypothesis for efficacy comes first:
							output_RejectH0[row] = output_Stat1[row] <= tail;
							output_Reduction[row] = red;
						}else{
							// Hypothesis for resistance comes second:
							output_RejectH0[row] = output_Stat2[row] <= tail;
							output_Reduction[row] = target_red;
						}

						// Increment the row number:
						row++;
					}
				}
			}
		}
	}

	return output;
}



RCPP_MODULE(rcpp_module){

	using namespace Rcpp;

	function("estimate_k_ml", &estimate_k_ml);
	function("RCPP_efficacy_analysis", &efficacy_analysis);
	function("RCPP_typology_analysis", &typology_analysis);
	function("RCPP_efficacy_frequencies_paired", &efficacy_frequencies_paired);
	function("RCPP_efficacy_frequencies_unpaired", &efficacy_frequencies_unpaired);
	function("RCPP_power_matrix_paired", &power_matrix_paired);

}


// [[Rcpp::export]]
Rcpp::DataFrame fecrt_sim_paired_LP(int iters, Rcpp::NumericVector red, int N, double mu, double cv_pre, double cv_post, double ncor, Rcpp::NumericMatrix thresholds, Rcpp::NumericVector conjugate_priors, int delta, int beta_iters, int approx, double tail, bool useml, Rcpp::NumericVector dobson_cl, Rcpp::NumericVector dobson_priors, bool uselp)
{

	// Check that the thresholds have 2 columns and >0 rows:
	if(thresholds.ncol() != 2){
		Rcpp::stop("There must be exactly 2 columns");
	}
	if(thresholds.nrow() < 1){
		Rcpp::stop("There must be 1 or more rows");
	}
	for(int i=0; i<thresholds.nrow(); i++){
		if(thresholds(i,1) < thresholds(i,0)){
			Rcpp::stop("Second thresholds must be higher than the first");
		}
		if(thresholds(i,0) < 0 || thresholds(i,0) > 1){
			Rcpp::stop("Thresholds must be between 0 and 1");
		}
		if(thresholds(i,1) < 0 || thresholds(i,1) > 1){
			Rcpp::stop("Thresholds must be between 0 and 1");
		}
	}

	// TODO: tidy up check and error:
	if((double(N)*double(mu)) > (0.1 * INT_MAX)){
		Rcpp::stop("Possible int overflow");
	}

	// Use the expand.grid R function:
	Rcpp::Function expGrid("expand.grid");

	Rcpp::StringVector methods = { "BNB", "Levecke", "WAAVP", "MLE", "Dobson" };
	Rcpp::IntegerVector tseq = Rcpp::seq(1L, thresholds.nrow());
	Rcpp::DoubleVector td = { NA_REAL };
	Rcpp::StringVector ts = { NA_STRING };
	//	Rcpp::IntegerVector ti = { NA_INTEGER };
	Rcpp::IntegerVector iteration = Rcpp::seq(1L, iters);

	Rcpp::DataFrame output = expGrid(methods, red, tseq, td, td, iteration, td, td, td, ts, Rcpp::Named("stringsAsFactors")=false);
	output.names() = Rcpp::StringVector::create("Method", "Reduction", "ThresholdIndex", "ThresholdLower", "ThresholdUpper", "Iteration", "ObsReduction", "Stat1", "Stat2", "Classification");

	// Get references to list items:
	Rcpp::StringVector output_Method = output[0L];
	Rcpp::DoubleVector output_Reduction = output[1L];
	Rcpp::IntegerVector output_ThresholdIndex = output[2L];
	Rcpp::DoubleVector output_ThresholdLower = output[3L];
	Rcpp::DoubleVector output_ThresholdUpper = output[4L];
	Rcpp::IntegerVector output_Iteration = output[5L];
	Rcpp::DoubleVector output_ObsReduction = output[6L];
	Rcpp::DoubleVector output_Stat1 = output[7L];
	Rcpp::DoubleVector output_Stat2 = output[8L];
	Rcpp::StringVector output_Classification = output[9L];

	// Note:  Rcpp allows conversion from long long but not to double array directly
	double conjugate_priors_db[2] = { conjugate_priors[0], conjugate_priors[1] };
	double dobson_priors_db[2] = { dobson_priors[0], dobson_priors[1] };

	// To control the row index to write to:
	int row = 0L;

	// How to simulate:
	// distribution(const std::array<const double, s_dim> mu, const std::array<const double, s_dim> cv, const double lcor)
	const std::array<const double, 2L> mus { mu, mu };
	const std::array<const double, 2L> cvs { cv_pre, cv_post };

	distribution<dists::mvlnormpois> distn_lp(mus, cvs, ncor);
	distribution<dists::mvnbinom> distn_gp(mus, cvs, ncor);
	// Note: use draw rather than draw_mu and multiply post-tx by reds, then do a Poisson draw
	// as this is the same thing as repeated lognormal draws and cheaper

	// First loop over correlation and pre-treatment datasets to simulate:
	for(int i=0; i<iters; i++){

		Rcpp::NumericVector mus_1(N);
		Rcpp::NumericVector mus_2(N);
		if(uselp)
		{
			for(int i=0; i<N; ++i)
			{
				std::array<double, 2L> rv = distn_lp.draw_mu();
				mus_1[i] = rv[0L];
				mus_2[i] = rv[1L];
			}
		}else{
			for(int i=0; i<N; ++i)
			{
				std::array<double, 2L> rv = distn_gp.draw_mu();
				mus_1[i] = rv[0L];
				mus_2[i] = rv[1L];
			}
		}

		Rcpp::IntegerVector pre_data(N);
		long long sum_pre = 0L;
		do {
			sum_pre = 0L;
			for(int i=0; i<N; i++){
				pre_data[i] = R::rpois(mus_1[i]);
				sum_pre += pre_data[i];
			}
		} while (sum_pre == 0L);

		double mean_pre = Rcpp::mean(pre_data);
		double var_pre = Rcpp::var(pre_data);

		// Then loop over reductions:
		for(int r=0; r<red.length(); r++){
			long long sum_post = 0L;
			Rcpp::IntegerVector post_data(N);
			for(int i=0; i<N; i++){
				post_data[i] = R::rpois(mus_2[i] * red[r]);
				sum_post += post_data[i];
			}

			double mean_post = Rcpp::mean(post_data);
			double var_post = Rcpp::var(post_data);
			double obsred = mean_post / mean_pre;

			// We need both scaled_cov for estimate_k and cov for the estimation functions (including bnb if using a large sample approximation)
			double scaled_cov = 0.0;
			double cov = 0.0;
			for(int i=0; i<N; i++){
				scaled_cov += (double(pre_data[i])/mean_pre - 1.0) * (double(post_data[i])/mean_post - 1.0);
				cov += (double(pre_data[i]) - mean_pre) * (double(post_data[i]) - mean_post);
			}
			scaled_cov = scaled_cov / double(N - 1L);
			cov = cov / double(N - 1L);

			Rcpp::NumericVector ks;
			if(useml){
				ks = estimate_k_ml(pre_data, mean_pre, var_pre, post_data, mean_post, var_post, cov, true);
			}else{
				ks = estimate_k(mean_pre, var_pre, mean_post, var_post, cov, true);
			}
			// Test:
			if(!R_finite(ks[0])){
				Rcpp::stop("Non-finite k1 generated");
			}
			if(ks[0] <= 0.0){
				Rcpp::stop("k1 below 0.0 generated");
			}
			if(!R_finite(ks[1])){
				Rcpp::stop("Non-finite k2 generated");
			}
			if(ks[1] <= 0.0){
				Rcpp::stop("k2 below 0.0 generated");
			}
			double estk_pre = ks[0];
			double estk_post = ks[1];

			// For cheating:
			//estk_pre = k_pre / (1.0 - (std::sqrt(k_pre) * std::sqrt(k_post) / k_c));
			//estk_post = k_post / (1.0 - (std::sqrt(k_pre) * std::sqrt(k_post) / k_c));
			estk_pre = 1.0 / (cv_pre*cv_pre);
			estk_post = 1.0 / (cv_post*cv_post);

			/*
			Rcpp::NumericVector uks = estimate_k(mean_pre, var_pre, mean_post, var_post, 0.0, false);
			double estk_pre_u = uks[0];
			double estk_post_u = uks[1];
			 */

			//			Rcpp::Rcout << cov << " - " << estk_pre << " - " << estk_pre_u << " - " << estk_post << " - " << estk_post_u << "\n";

			// Then loop over thresholds:
			for(int t=0; t<thresholds.nrow(); t++){

				// This stays the same for all rows within this threshold:
				int tind = output_ThresholdIndex[row] - 1L;
				double th1 = thresholds(tind,0L);
				double th2 = thresholds(tind,1L);


				// TODO: Method should not need to be set but do as a debug check?

				// Apply each test and save the results:

				output_Method[row] = "BNB";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				int fl = bnb_pval(sum_pre, N, estk_pre, mean_pre, var_pre, sum_post, N, estk_post, mean_post, var_post, cov, 1.0, th1, th2, conjugate_priors_db, delta, beta_iters, approx, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = get_type_pv(1.0-obsred, output_Stat1[row], output_Stat2[row], tail, th1, th2);
				row++;

				output_Method[row] = "Levecke";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				levecke_p_ci(mean_pre, mean_post, var_pre, var_post, cov, N, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post==0) ? "SumPost0" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "WAAVP";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				waavp_p_ci(mean_pre, mean_post, var_pre, var_post, cov, N, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post==0) ? "SumPost0" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "MLE";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				mle_p_ci(mean_pre, mean_post, var_pre, var_post, cov, N, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post==0) ? "SumPost0" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "Dobson";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				dobson_ci(sum_pre, sum_post, dobson_cl[0], dobson_cl[1], dobson_priors_db, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post>sum_pre) ? "post>pre" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				Rcpp::checkUserInterrupt();

			}
		}
	}

	return output;
}
