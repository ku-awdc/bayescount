// This needs to be RcppDist due to depdendencies on distribution.h
#include <RcppDist.h>

Rcpp::NumericMatrix draw_lambda(const int n, const Rcpp::NumericVector& mu, const Rcpp::NumericVector& cv, const double lcor, const std::string& dist);
Rcpp::IntegerMatrix draw_count(const int n, const Rcpp::NumericVector& mu, const Rcpp::NumericVector& cv, const double lcor, const std::string& dist);

Rcpp::NumericVector estimate_fecrt(Rcpp::IntegerVector& pre, Rcpp::IntegerVector& post, const bool paired, const std::string& k_type,
                                    const double mean_ratio, const double H0_1, const double H0_2, const double tail,
                                    const Rcpp::NumericVector& conjugate_priors, const std::string& delta, const int beta_iters,
                                    const std::string& approx, const Rcpp::NumericVector& dobson_priors, const double true_effk_pre,
                                    const double trye_effk_post);

Rcpp::NumericVector summarise_fecrt(Rcpp::IntegerVector& pre, Rcpp::IntegerVector& post, const bool paired, const std::string& k_type);
