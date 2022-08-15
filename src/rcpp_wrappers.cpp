#include "rcpp_wrappers.h"

// Note: we can't include Rcpp.h directly as distribution includes RcppDist.h
#include <array>
#include <algorithm>  // For std::max(size_t, size_t)

#include "distribution.h"
#include "simulation.h"
#include "estimator.h"
#include "enums.h"

#include "fecrt_summariser.h"
#include "fecrt_estimator.h"


Rcpp::NumericMatrix draw_lambda(const int n, const Rcpp::NumericVector& mu, const Rcpp::NumericVector& cv, const double lcor, const std::string& dist)
{
  Rcpp::NumericMatrix rv;
  if(dist=="mvlnormpois")
  {
    rv = draw_lambda_template<dists::mvlnormpois>(n, mu, cv, lcor);
  }
  else if(dist=="lnormpois")
  {
    rv = draw_lambda_template<dists::lnormpois>(n, mu, cv, lcor);
  }
  else if(dist=="nbinom")
  {
    rv = draw_lambda_template<dists::nbinom>(n, mu, cv, lcor);
  }
  else if(dist=="mvnbinom")
  {
    rv = draw_lambda_template<dists::mvnbinom>(n, mu, cv, lcor);
  }
  else if(dist=="poisson")
  {
    rv = draw_lambda_template<dists::poisson>(n, mu, cv, lcor);
  }
  else
  {
    Rcpp::stop("Unrecognised dist");
  }

  return rv;
}


Rcpp::IntegerMatrix draw_count(const int n, const Rcpp::NumericVector& mu, const Rcpp::NumericVector& cv, const double lcor, const std::string& dist)
{
  Rcpp::IntegerMatrix rv;
  if(dist=="mvlnormpois")
  {
    rv = draw_count_template<dists::mvlnormpois>(n, mu, cv, lcor);
  }
  else if(dist=="lnormpois")
  {
    rv = draw_count_template<dists::lnormpois>(n, mu, cv, lcor);
  }
  else if(dist=="nbinom")
  {
    rv = draw_count_template<dists::nbinom>(n, mu, cv, lcor);
  }
  else if(dist=="mvnbinom")
  {
    rv = draw_count_template<dists::mvnbinom>(n, mu, cv, lcor);
  }
  else if(dist=="poisson")
  {
    rv = draw_count_template<dists::poisson>(n, mu, cv, lcor);
  }
  else
  {
    Rcpp::stop("Unrecognised dist");
  }

  return rv;
}

Rcpp::NumericVector summarise_fecrt(Rcpp::IntegerVector& pre, Rcpp::IntegerVector& post, const bool paired, const std::string& k_type)
{

	int maxN = 10L;
	double true_effk_pre = 0.0;
	double true_effk_post = 0.0;

	// template<bool t_paired, ktypes t_ktype, typename t_cont_type, containers t_container>
	// fecrt_summarise(const size_t maxN, const double true_effk_pre, const double true_effk_post)

	Rcpp::NumericVector rv;

	if (paired)
	{
		fecrt_summariser<true, ktypes::flex, Rcpp::IntegerVector, containers::rcppvector>
		estr(maxN, k_type, true_effk_pre, true_effk_post);
		estr.push_multiple(pre, post);
		rv = estr.summarise();
	} else {
		fecrt_summariser<false, ktypes::flex, Rcpp::IntegerVector, containers::rcppvector>
		estr(maxN, k_type, true_effk_pre, true_effk_post);
		estr.push_multiple(pre, post);
		rv = estr.summarise();
	}

	return rv;
}


Rcpp::NumericVector estimate_fecrt(Rcpp::IntegerVector& pre, Rcpp::IntegerVector& post, const bool paired, const std::string& k_type,
                                   const double mean_ratio, const double H0_1, const double H0_2, const double tail,
                                   const Rcpp::NumericVector& conjugate_priors, const std::string& delta, const int beta_iters,
                                   const std::string& approx, const Rcpp::NumericVector& dobson_priors,
                                   const double true_effk_pre, const double true_effk_post)
{

	if(conjugate_priors.length()!=2L) Rcpp::stop("Invalid length conjugate_priors");
	const std::array<double, 2L> conjugate_priors_arr = { conjugate_priors[0L], conjugate_priors[1L] };

	if(dobson_priors.length()!=2L) Rcpp::stop("Invalid length dobson_priors");
	const std::array<double, 2L> dobson_priors_arr = { dobson_priors[0L], dobson_priors[1L] };

	optswitch deltaenum;
	if(delta=="never")
	{
		deltaenum = optswitch::never;
	}
	else if(delta=="sometimes")
	{
		deltaenum = optswitch::sometimes;
	}
	else if(delta=="always")
	{
		deltaenum = optswitch::always;
	}
	else
	{
		Rcpp::stop("Unhandled optswitch for delta");
	}

	optswitch approxenum;
	if(approx=="never")
	{
		approxenum = optswitch::never;
	}
	else if(approx=="sometimes")
	{
		approxenum = optswitch::sometimes;
	}
	else if(approx=="always")
	{
		approxenum = optswitch::always;
	}
	else
	{
		Rcpp::stop("Unhandled optswitch for approx");
	}

	int maxN = 0L;

	if (paired)
	{
		maxN = pre.size();
		if(post.size() != maxN)
		{
			Rcpp::stop("Unequal pre and post size for paired data");
		}
	}
	else
	{
		maxN = std::max(pre.size(), post.size());
	}

	Rcpp::NumericVector rv;

	// template<bool t_all_methods, bool t_paired, ktypes t_ktype, typename t_cont_type, containers t_container>
	// class fecrt_estimator(const double mean_ratio, const double H0_1, const double H0_2,
	// const double tail, const std::array<double, 2L> conjugate_priors, const optswitch delta,
	// const int beta_iters, const optswitch approx, const std::array<double, 2L> dobson_priors,
	// const size_t maxN, const std::string& k_type, const double true_effk_pre, const double true_effk_post)

	if (paired)
	{
		fecrt_estimator<true, true, ktypes::flex, Rcpp::IntegerVector, containers::rcppvector>
			estimator(mean_ratio, H0_1, H0_2, tail, conjugate_priors_arr, deltaenum, beta_iters, approxenum,
    	dobson_priors_arr, maxN, k_type, true_effk_pre, true_effk_post);
		rv = estimator.push_multiple(pre, post);
	}

	return rv;
}

/*
Rcpp::NumericVector estimate_fecrt(Rcpp::IntegerVector& pre, Rcpp::IntegerVector& post, const bool paired, const std::string& k_type,
                                    const double mean_ratio, const double H0_1, const double H0_2, const double tail,
                                    const Rcpp::NumericVector& conjugate_priors, const std::string& delta, const int beta_iters,
                                    const std::string& approx, const Rcpp::NumericVector& dobson_priors, const double true_effk_pre,
                                    const double true_effk_post)
{
  constexpr size_t s_rvlen_all = 11L;
  constexpr size_t s_rvlen_fix = 3L;

  if(conjugate_priors.length()!=2L) Rcpp::stop("Invalid length conjugate_priors");
  const std::array<double, 2L> conjugate_priors_arr = { conjugate_priors[0L], conjugate_priors[1L] };

  if(dobson_priors.length()!=2L) Rcpp::stop("Invalid length dobson_priors");
  const std::array<double, 2L> dobson_priors_arr = { dobson_priors[0L], dobson_priors[1L] };

  optswitch deltaenum;
  if(delta=="never")
  {
    deltaenum = optswitch::never;
  }
  else if(delta=="sometimes")
  {
    deltaenum = optswitch::sometimes;
  }
  else if(delta=="always")
  {
    deltaenum = optswitch::always;
  }
  else
  {
    Rcpp::stop("Unhandled optswitch for delta");
  }

  optswitch approxenum;
  if(approx=="never")
  {
    approxenum = optswitch::never;
  }
  else if(approx=="sometimes")
  {
    approxenum = optswitch::sometimes;
  }
  else if(approx=="always")
  {
    approxenum = optswitch::always;
  }
  else
  {
    Rcpp::stop("Unhandled optswitch for approx");
  }

	/ *
  template<bool t_paired, bool t_all_methods, ktypes t_ktype, typename t_cont_type, containers t_container, size_t t_rvlen>
	estimator(const size_t maxN, const double mean_ratio, const double H0_1, const double H0_2,
		const double tail, const std::array<double, 2L> conjugate_priors, const optswitch delta,
		const int beta_iters, const optswitch approx, const std::array<double, 2L> dobson_priors,
		const double true_effk_pre, const double true_effk_post)
	 * /

  Rcpp::NumericVector rv;

  size_t maxN = pre.length();
  if(paired)
  {
    if(maxN != post.length()) Rcpp::stop("Unequal length pre and post treatment data for paired analysis");

    if(k_type=="fix")
    {
      estimator<true, false, ktypes::fix, Rcpp::IntegerVector, containers::rcppvector, s_rvlen_fix>
        estr(maxN, mean_ratio, H0_1, H0_2, tail, conjugate_priors_arr, deltaenum, beta_iters,
              approxenum, dobson_priors_arr, true_effk_pre, true_effk_post);

      estr.push_multiple(pre, post);
      rv = estr.estimate();
    }
    else if(k_type=="mm")
    {
      estimator<true, true, ktypes::mm, Rcpp::IntegerVector, containers::rcppvector, s_rvlen_all>
        estr(maxN, mean_ratio, H0_1, H0_2, tail, conjugate_priors_arr, deltaenum, beta_iters,
              approxenum, dobson_priors_arr, true_effk_pre, true_effk_post);

      estr.push_multiple(pre, post);
      rv = estr.estimate();
    }
    else if(k_type=="ql")
    {
      estimator<true, true, ktypes::ql, Rcpp::IntegerVector, containers::rcppvector, s_rvlen_all>
        estr(maxN, mean_ratio, H0_1, H0_2, tail, conjugate_priors_arr, deltaenum, beta_iters,
              approxenum, dobson_priors_arr, true_effk_pre, true_effk_post);

      estr.push_multiple(pre, post);
      rv = estr.estimate();
    }
    else if(k_type=="ml")
    {
      estimator<true, true, ktypes::ml, Rcpp::IntegerVector, containers::rcppvector, s_rvlen_all>
        estr(maxN, mean_ratio, H0_1, H0_2, tail, conjugate_priors_arr, deltaenum, beta_iters,
              approxenum, dobson_priors_arr, true_effk_pre, true_effk_post);

      estr.push_multiple(pre, post);
      rv = estr.estimate();
    }
    else
    {
      Rcpp::stop("Unhandled k_type");
    }

  }
  else
  {
    maxN = std::max<size_t>(maxN, post.length());

    if(k_type=="fix")
    {
      estimator<false, false, ktypes::fix, Rcpp::IntegerVector, containers::rcppvector, s_rvlen_fix>
        estr(maxN, mean_ratio, H0_1, H0_2, tail, conjugate_priors_arr, deltaenum, beta_iters,
              approxenum, dobson_priors_arr, true_effk_pre, true_effk_post);

      estr.push_multiple(pre, post);
      rv = estr.estimate();
    }
    else if(k_type=="mm")
    {
      estimator<false, true, ktypes::mm, Rcpp::IntegerVector, containers::rcppvector, s_rvlen_all>
        estr(maxN, mean_ratio, H0_1, H0_2, tail, conjugate_priors_arr, deltaenum, beta_iters,
              approxenum, dobson_priors_arr, true_effk_pre, true_effk_post);

      estr.push_multiple(pre, post);
      rv = estr.estimate();
    }
    else if(k_type=="ql")
    {
      estimator<false, true, ktypes::ql, Rcpp::IntegerVector, containers::rcppvector, s_rvlen_all>
        estr(maxN, mean_ratio, H0_1, H0_2, tail, conjugate_priors_arr, deltaenum, beta_iters,
              approxenum, dobson_priors_arr, true_effk_pre, true_effk_post);

      estr.push_multiple(pre, post);
      rv = estr.estimate();
    }
    else if(k_type=="ml")
    {
      estimator<false, true, ktypes::ml, Rcpp::IntegerVector, containers::rcppvector, s_rvlen_all>
        estr(maxN, mean_ratio, H0_1, H0_2, tail, conjugate_priors_arr, deltaenum, beta_iters,
              approxenum, dobson_priors_arr, true_effk_pre, true_effk_post);

      estr.push_multiple(pre, post);
      rv = estr.estimate();
    }
    else
    {
      Rcpp::stop("Unhandled k_type");
    }

  }

  return rv;

}
*/
