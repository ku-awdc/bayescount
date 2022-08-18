#ifndef FECRT_ESTIMATOR_H
#define FECRT_ESTIMATOR_H

// This is required by distribution.h and precludes Rcpp.h:
#include <RcppDist.h>
#include <array>
#include <type_traits>
#include <typeinfo>

#include "enums.h"
#include "fecrt_summariser.h"

// The estimator class just produces CI and p-values for the different methods from summary stats
// But it contains a summariser that produces means/vars/ks (and potentially stores data)

// Note: only the first of these template parameters is used here, the rest are passed down to
// the summariser class
template<bool t_all_methods, bool t_paired, ktypes t_ktype, typename t_cont_type, containers t_container>
class fecrt_estimator
{
private:

	fecrt_summariser<t_paired, t_ktype, t_cont_type, t_container> m_summariser;

  const double m_mean_ratio;
  const double m_H0_1;
  const double m_H0_2;
  const double m_tail;
  const double m_lci;
  const double m_uci;
  const std::array<double, 2L> m_conjugate_priors;
  const optswitch m_delta;
  const int m_beta_iters;
  const optswitch m_approx;
  const std::array<double, 2L> m_dobson_priors;

  static constexpr size_t s_summary_len = t_paired ? 12L : 8L;
  // double(Ns), mus, vars, ks (adjusted), if(t_paired): cov, cor, ks (unadjusted)

  static constexpr size_t s_estimate_len = t_all_methods ? (t_paired ? 12L : 10L) : 4L;
	// mus, BNB pvals, if(t_all_methods) MLE, WAAVP, Levecke, if(t_paired) Dobson

  // TODO: overload paired vs not paired on arguments or templates?

  /*
  std::array<double, 2L> estimate_mle(const double sum_pre, const double sum_post,
                                      const double var_pre, const double var_post) const
  {
  	// TODO: change from Dobson

  	std::array<double, 2L> rv;

    if(var_pre <= 0.0 || var_post <= 0.0)
    {
      rv[0L] = NA_REAL;
    	rv[1L] = NA_REAL;
    }
    else if constexpr(t_paired)
    {
      if(sum_post > sum_pre)
      {
      	rv[0L] = NA_REAL;
      	rv[1L] = NA_REAL;
      }
      else
      {
      	const double postsumd = static_cast<double>(sum_post);
      	const double sumdiffd = static_cast<double>(sum_pre - sum_post);
      	rv[0L] = 1.0 - R::qbeta(m_lci, postsumd+m_dobson_priors[0L], sumdiffd+m_dobson_priors[1L], 0L, 0L);
      	rv[1L] = 1.0 - R::qbeta(m_uci, postsumd+m_dobson_priors[0L], sumdiffd+m_dobson_priors[1L], 0L, 0L);
      }
    }
    else
    {
    	rv[0L] = NA_REAL;
    	rv[1L] = NA_REAL;
    }
  }
  */

  std::array<double, 2L> estimate_waavp(const double N_pre, const double N_post,
                                        const double mean_pre, const double mean_post,
                                        const double var_pre, const double var_post,
                                        const double cov) const
  {
  	std::array<double, 2L> rv;

    if(var_pre <= 0.0 || var_post <= 0.0)
    {
      rv[0L] = NA_REAL;
    	rv[1L] = NA_REAL;
    }
    else if constexpr(t_paired)
    {
  	// Method B of Lyndal-Murphy, M., Swain, a J., & Pepper, P. M. (2014). Methods to determine resistance to anthelmintics when continuing larval development occurs. Veterinary Parasitology, 199(3–4), 191–200. https://doi.org/10.1016/j.vetpar.2013.11.002

    	const double df = N_pre - 1.0;
    	// Signature of qt is:  double  qt(double, double, int, int);
    	const double tval = R::qt(1.0 - m_tail, df, 1L, 0L);

    	const double Nd = N_pre;  // Will be equal to N_post
    	// const double varred = var1 / (Nd * mu1 * mu1) + var2 / (Nd * mu2 * mu2) - 2.0 * cov12 / (Nd * mu1 * mu2);
    	const double varred = var_pre / (Nd * mean_pre * mean_pre) + var_post / (Nd * mean_post * mean_post) - 2.0 * cov / (Nd * mean_pre * mean_post);

    	rv[1L] = 1.0 - (mean_post / mean_pre * std::exp(-tval * std::sqrt(varred)));
    	rv[0L] = 1.0 - (mean_post / mean_pre * std::exp(tval * std::sqrt(varred)));
    }
    else
    {
    	// Method A of Lyndal-Murphy, M., Swain, a J., & Pepper, P. M. (2014). Methods to determine resistance to anthelmintics when continuing larval development occurs. Veterinary Parasitology, 199(3–4), 191–200. https://doi.org/10.1016/j.vetpar.2013.11.002

    	const double df = N_pre + N_post - 2.0;
    	// Signature of qt is:  double  qt(double, double, int, int);
    	const double tval = R::qt(1.0 - m_tail, df, 1L, 0L);

    	// const double varred = var1 / (N1d * mu1 * mu1) + var2 / (N2d * mu2 * mu2);
    	const double varred = var_pre / (N_pre * mean_pre * mean_pre) + var_post / (N_post * mean_post * mean_post);

    	rv[1L] = 1.0 - (mean_post / mean_pre * std::exp(-tval * std::sqrt(varred)));
    	rv[0L] = 1.0 - (mean_post / mean_pre * std::exp(tval * std::sqrt(varred)));
    }

    return rv;
  }

  /*
  void estimate_levecke(const double N_pre, const double N_post,
   const double mean_pre, const double mean_post,
   const double var_pre, const double var_post,
   const double cov, double& lci, double& uci) const
  {
   // void levecke_u_ci(double mu1, double mu2, double var1, double var2, int N1, int N2, double tail, double *ci_l, double *ci_u);
   // void levecke_p_ci(double mu1, double mu2, double var1, double var2, double cov12, int N, double tail, double *ci_l, double *ci_u);

    if(var_pre <= 0.0 || var_post <= 0.0)
    {
      uci = NA_REAL;
      lci = NA_REAL;
    }
    else if constexpr(t_paired)
    {
      if(sum_post > sum_pre)
      {
        lci = NA_REAL;
        uci = NA_REAL;
      }
      else
      {
      	const double postsumd = static_cast<double>(sum_post);
      	const double sumdiffd = static_cast<double>(sum_pre - sum_post);
      	lci = 1.0 - R::qbeta(m_lci, postsumd+m_dobson_priors[0L], sumdiffd+m_dobson_priors[1L], 0L, 0L);
      	uci = 1.0 - R::qbeta(m_uci, postsumd+m_dobson_priors[0L], sumdiffd+m_dobson_priors[1L], 0L, 0L);
      }
    }
    else
    {
      lci = NA_REAL;
      uci = NA_REAL;
    }
  }
  */

  std::array<double, 2L> estimate_dobson(const double sum_pre, const double sum_post,
                                         const double var_pre, const double var_post) const
  {
  	std::array<double, 2L> rv;

  	if(var_pre <= 0.0 || var_post <= 0.0)
  	{
  		rv[0L] = NA_REAL;
  		rv[1L] = NA_REAL;
  	}
  	else if constexpr(t_paired)
  	{
  		if(sum_post > sum_pre)
  		{
  			rv[0L] = NA_REAL;
  			rv[1L] = NA_REAL;
  		}
  		else
  		{
  			const double postsumd = static_cast<double>(sum_post);
  			const double sumdiffd = static_cast<double>(sum_pre - sum_post);
  			rv[0L] = 1.0 - R::qbeta(m_lci, postsumd+m_dobson_priors[0L], sumdiffd+m_dobson_priors[1L], 0L, 0L);
  			rv[1L] = 1.0 - R::qbeta(m_uci, postsumd+m_dobson_priors[0L], sumdiffd+m_dobson_priors[1L], 0L, 0L);
  		}
  	}
  	else
  	{
  		rv[0L] = NA_REAL;
  		rv[1L] = NA_REAL;
  	}
  }


public:
	fecrt_estimator(const double mean_ratio, const double H0_1, const double H0_2,
    const double tail, const std::array<double, 2L> conjugate_priors, const optswitch delta,
    const int beta_iters, const optswitch approx, const std::array<double, 2L> dobson_priors,
    // Arguments passed stright to summariser:
    const size_t maxN, const std::string& k_type, const double true_effk_pre, const double true_effk_post) :
    m_mean_ratio(mean_ratio), m_H0_1(H0_1), m_H0_2(H0_2),
    m_tail(tail), m_lci(tail), m_uci(1.0 - tail), m_conjugate_priors(conjugate_priors), m_delta(delta),
    m_beta_iters(beta_iters), m_approx(approx), m_dobson_priors(dobson_priors),
    m_summariser(maxN, k_type, true_effk_pre, true_effk_post)
  {

  }

  void push_pre(const int obs)
  {
  	m_summariser.push_pre(obs);
  }

  void push_post(const int obs)
  {
  	m_summariser.push_post(obs);
  }

  void push_pair(const int pre, const int post)
  {
  	m_summariser.push_pair(pre, post);
  }

  Rcpp::NumericVector push_multiple(const Rcpp::IntegerVector& pre, const Rcpp::IntegerVector& post,
                                    const bool return_estimate = true)
  {
  	{
  		Rcpp::NumericVector summaries = m_summariser.push_multiple(pre, post);
  		if (!return_estimate)
  		{
  			return summaries;
  		}
  	}

  	std::array<double, s_estimate_len> ests = estimate();
  	Rcpp::NumericVector rv(s_estimate_len);
  	for(size_t i=0L; i<s_estimate_len; ++i)
  	{
  		rv[i] = ests[i];
  	}

  	return rv;
  }

	std::array<double, s_estimate_len> estimate()
	{
		std::array<double, s_summary_len> summaries = m_summariser.summarise();
		std::array<double, s_estimate_len> rv;

		// Copy means:
		rv[0L] = summaries[0L];
		rv[1L] = summaries[1L];

		// Do BNB:

		if constexpr(t_all_methods)
		{
			// Do MLE:
			// Do WAAVP:
			// Do Levecke:

			if constexpr(t_paired)
			{
				// Do Dobson:
			}
		}

		return rv;
	}

	/*
	void estimate(double& efficacy, double& p1, double& p2)
  {
  	std::array<double, 6L> sv = m_summariser.summarise_adjusted();


    // Update internally stored variance:
    if constexpr(t_paired)
    {
      var_pre = m_varnum_pre / static_cast<double>(m_N-1L);
      var_post = m_varnum_post / static_cast<double>(m_N-1L);
      m_cov = m_covnum / static_cast<double>(m_N-1L);
    }
    else
    {
      var_pre = m_varnum_pre / static_cast<double>(m_Npre-1L);
      var_post = m_varnum_post / static_cast<double>(m_Npost-1L);
    }

    // Estimate k if necessary:
    if constexpr(t_ktype == ktypes::fix)
    {
      m_effk_pre = m_true_effk_pre;
      m_effk_post = m_true_effk_post;
    }
    else if constexpr(t_ktype == ktypes::mm)
    {
      m_effk_pre = m_true_effk_pre;
      m_effk_post = m_true_effk_post;
    }
    else
    {
      Rcpp::stop("Unhandled t_ktype in estimate()");
    }

    // First value is always efficacy:
    efficacy = 1.0 - m_mean_post / m_mean_pre;

    // TODO: re-write bnb_pval as internal method
    int delta = static_cast<int>(m_delta);
    int approx = static_cast<int>(m_approx);
    double conjugate_priors[2] = { m_conjugate_priors[0L], m_conjugate_priors[1L] };
    double p_1 = 0.0;
    double p_2 = 0.0;

    int n_pre;
    int n_post;
    if constexpr(t_paired)
    {
      n_pre = m_Npre;
      n_post = m_Npost;
    }
    else
    {
      n_pre = m_N;
      n_post = m_N;
    }

    int bnbfailed = bnb_pval(sum_pre, n_pre, m_effk_pre, m_mean_pre, var_pre, sum_post, n_post, m_effk_post, m_mean_post, var_post, m_cov,
                            m_mean_ratio, m_H0_1, m_H0_2, conjugate_priors, delta, m_beta_iters, approx, &p_1, &p_2);
    p1 = p_1;
    p2 = p_2;

  }
	 */

	/*
  void estimate(double& efficacy, double& bnb_p1, double& bnb_p2, double& mle_lci, double& mle_uci, double& waavp_lci, double& waavp_uci,
                double& levecke_lci, double& levecke_uci, double& dobson_lci, double& dobson_uci)
  {
    if constexpr(!t_all_methods)
    {
      Rcpp::stop("Unable to call full estimate function when !t_all_methods");
    }

    // Delegate var/cov/k estimation and Denwood method:
    estimate(efficacy, bnb_p1, bnb_p2);

    // And then do the other methods:
    estimate_mle(mle_lci, mle_uci);
    estimate_waavp(waavp_lci, waavp_uci);
    estimate_levecke(levecke_lci, levecke_uci);
    estimate_dobson(dobson_lci, dobson_uci);
  }

  Rcpp::NumericVector estimate()
  {
    if(m_N+m_Npre < 0L )
    {
      Rcpp::stop("Can't estimate: no pre-treatment data!");
    }
    if(m_N+m_Npost == 0L)
    {
      Rcpp::stop("Can't estimate: no post-treatment data!");
    }

    // To store output (append effk/variance/covariance for R):
    Rcpp::NumericVector rv(t_rvlen+5L);

    // Get esimtates:
    if constexpr(t_all_methods)
    {
      estimate(rv[0L], rv[1L], rv[2L], rv[3L], rv[4L], rv[5L], rv[6L], rv[7L], rv[8L], rv[9L], rv[10L]);
      static_assert(t_rvlen==11L, "Logic error: incorrect argument length for estimate(...)");

      rv.attr("names") = Rcpp::CharacterVector::create("efficacy", "p1", "p2", "mle_lci", "mle_uci", "waavp_lci", "waavp_uci",
                                          "levecke_lci", "levecke_uci", "dobson_lci", "dobson_uci", "effk1", "effk2", "var1", "var2", "cov");
    }
    else
    {
      estimate(rv[0L], rv[1L], rv[2L]);

      rv.attr("names") = Rcpp::CharacterVector::create("efficacy", "p1", "p2", "effk1", "effk2", "var1", "var2", "cov");
    }

    rv[t_rvlen] = m_effk_pre;
    rv[t_rvlen+1L] = m_effk_post;
    rv[t_rvlen+2L] = var_pre;
    rv[t_rvlen+3L] = var_post;
    rv[t_rvlen+4L] = m_cov;

    return rv;
  }
	*/

};


#endif // FECRT_ESTIMATOR_H
