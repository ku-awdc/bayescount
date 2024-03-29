#ifndef ESTIMATOR_H
#define ESTIMATOR_H

// This is required by distribution.h and precludes Rcpp.h:
#include <RcppDist.h>
#include <array>
#include <type_traits>
#include <typeinfo>

#include "enums.h"
#include "fecrt.h"
//#include "data_storage.h"

// TODO: an estimater class that sets up the estimation method and potentially stores pre and post internally
// Then this gets passed to simulation, so that simulation doesn't need to store counts or do any estimating itself
// estimator class could have methods to add data 1 point at a time (running mean) or as an array or as a vector
// this means the newton rhapson function should be templated (or included as a method within estimator class??)

// TODO: fecrt_estimator does the estimation methods (including all BNB stuff)
//        fecrt_summariser does the count storage and mean/var/cov/k estimation


template<bool t_paired, bool t_all_methods, ktypes t_ktype, typename t_cont_type, containers t_container, size_t t_rvlen>
class estimator
{
private:
  const size_t m_maxN;
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
  // Should be pre-calculated as effective k for paired:
  const double m_true_effk_pre;
  const double m_true_effk_post;

  // Note: these are not actually needed if !s_store_count
  t_cont_type m_pre;
  t_cont_type m_post;

  unsigned long long m_sum_pre = 0L;
  unsigned long long m_sum_post = 0L;

  double m_mean_pre = 0.0;
  double m_mean_post = 0.0;
  double m_varnum_pre = 0.0;
  double m_varnum_post = 0.0;
  double m_covnum = 0.0;

  // TODO: nest std::conditional to get correct child class for fecrt_summariser (paired or not)?
  // Temporary solution - will be replaced with data_storage class:
  using type_paired = std::conditional_t<t_paired, size_t, const size_t>;
  using type_unpaired = std::conditional_t<t_paired, const size_t, size_t>;
  type_paired m_N = -1L;
  type_unpaired m_Npre = -1L;
  type_unpaired m_Npost = -1L;

  double m_var_pre = NA_REAL;
  double m_var_post = NA_REAL;
  double m_cov = NA_REAL;
  double m_effk_pre = NA_REAL;
  double m_effk_post = NA_REAL;

  static constexpr bool s_store_count = (t_ktype==ktypes::ql || t_ktype==ktypes::ml);

  void add_pre(const int obs)
  {
    if constexpr(t_paired)
    {
      Rcpp::stop("Logic error: add_pre called for a paired analysis");
    }

    // If needed, store the data:
    if constexpr(s_store_count)
    {
      if(m_Npre >= m_maxN)
      {
        Rcpp::stop("Max storage capacity exceeded for pre-treatment data");
      }
      m_pre[m_Npre] = obs;
    }

    m_Npre++;

		double delta = obs - m_mean_pre;
    m_mean_pre += delta / static_cast<double>(m_Npre);
		m_varnum_pre += delta * (obs - m_mean_pre);
    // Variance = m_varnum_pre / (m_Npre-1.0)

    m_sum_pre += obs;

  }

  void add_post(const int obs)
  {
    if constexpr(t_paired)
    {
      Rcpp::stop("Logic error: add_post called for a paired analysis");
    }

    // If needed, store the data:
    if constexpr(s_store_count)
    {
      if(m_Npost >= m_maxN)
      {
        Rcpp::stop("Max storage capacity exceeded for post-treatment data");
      }
      m_post[m_Npost] = obs;
    }

    m_Npost++;

		double delta = obs - m_mean_post;
    m_mean_post += delta / static_cast<double>(m_Npost);
		m_varnum_post += delta * (obs - m_mean_post);
    // Variance = m_varnum_post / (m_Npost-1.0)

    m_sum_post += obs;

  }

  void add_pair(const int pre, const int post)
  {
    if constexpr(!t_paired)
    {
      add_pre(pre);
      add_post(post);
    }
    else
    {
      // If needed, store the data:
      if constexpr(s_store_count)
      {
        if(m_N > m_maxN)
        {
          Rcpp::stop("Max storage capacity exceeded");
        }
        m_pre[m_N] = pre;
        m_post[m_N] = post;
      }

      m_N++;

  		double delta1 = pre - m_mean_pre;
  		double delta2 = post - m_mean_post;

      m_mean_pre += delta1 / static_cast<double>(m_N);
  		m_varnum_pre += delta1 * (pre - m_mean_pre);
      // Variance = m_varnum_pre / (m_N-1.0)

      m_mean_post += delta2 / static_cast<double>(m_N);
  		m_varnum_post += delta2 * (post - m_mean_post);
      // Variance = m_varnum_post / (m_N-1.0)

      const double iip1 = static_cast<double>(m_N-1L) / static_cast<double>(m_N);
      m_covnum += delta1 * delta2 * iip1;
      // Coariance = m_covnum / (m_N-1.0)

      m_sum_pre += pre;
      m_sum_post += post;

    }
  }

  void estimate_mle(double& lci, double& uci) const
  {
    if(m_var_pre <= 0.0 || m_var_post <= 0.0)
    {
      uci = NA_REAL;
      lci = NA_REAL;
    }
    else if constexpr(t_paired)
    {
      if(m_sum_post > m_sum_pre)
      {
        lci = NA_REAL;
        uci = NA_REAL;
      }
      else
      {
      	const double postsumd = static_cast<double>(m_sum_post);
      	const double sumdiffd = static_cast<double>(m_sum_pre - m_sum_post);
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

  void estimate_waavp(double& lci, double& uci) const
  {
    if(m_var_pre <= 0.0 || m_var_post <= 0.0)
    {
      uci = NA_REAL;
      lci = NA_REAL;
    }
    else if constexpr(t_paired)
    {
  	// Method B of Lyndal-Murphy, M., Swain, a J., & Pepper, P. M. (2014). Methods to determine resistance to anthelmintics when continuing larval development occurs. Veterinary Parasitology, 199(3–4), 191–200. https://doi.org/10.1016/j.vetpar.2013.11.002

    	const double df = static_cast<double>(m_N - 1L);
    	// Signature of qt is:  double  qt(double, double, int, int);
    	const double tval = R::qt(1.0 - m_tail, df, 1L, 0L);

    	const double Nd = static_cast<double>(m_N);
    	// const double varred = var1 / (Nd * mu1 * mu1) + var2 / (Nd * mu2 * mu2) - 2.0 * cov12 / (Nd * mu1 * mu2);
    	const double varred = m_var_pre / (Nd * m_mean_pre * m_mean_pre) + m_var_post / (Nd * m_mean_post * m_mean_post) - 2.0 * m_cov / (Nd * m_mean_pre * m_mean_post);

    	uci = 1.0 - (m_mean_post / m_mean_pre * std::exp(-tval * std::sqrt(varred)));
    	lci = 1.0 - (m_mean_post / m_mean_pre * std::exp(tval * std::sqrt(varred)));
    }
    else
    {
    	// Method A of Lyndal-Murphy, M., Swain, a J., & Pepper, P. M. (2014). Methods to determine resistance to anthelmintics when continuing larval development occurs. Veterinary Parasitology, 199(3–4), 191–200. https://doi.org/10.1016/j.vetpar.2013.11.002
	
    	const double df = static_cast<double>(m_Npre + m_Npost - 2L);
    	// Signature of qt is:  double  qt(double, double, int, int);
    	const double tval = R::qt(1.0 - m_tail, df, 1L, 0L);
	
    	// const double varred = var1 / (N1d * mu1 * mu1) + var2 / (N2d * mu2 * mu2);
    	const double varred = m_var_pre / (static_cast<double>(m_Npre) * m_mean_pre * m_mean_pre) + m_var_post / (static_cast<double>(m_Npost) * m_mean_post * m_mean_post);
	
    	uci = 1.0 - (m_mean_post / m_mean_pre * std::exp(-tval * std::sqrt(varred)));
    	lci = 1.0 - (m_mean_post / m_mean_pre * std::exp(tval * std::sqrt(varred)));
    }
  }
  
  void estimate_levecke(double& lci, double& uci) const
  {
    if(m_var_pre <= 0.0 || m_var_post <= 0.0)
    {
      uci = NA_REAL;
      lci = NA_REAL;
    }
    else if constexpr(t_paired)
    {
      if(m_sum_post > m_sum_pre)
      {
        lci = NA_REAL;
        uci = NA_REAL;
      }
      else
      {
      	const double postsumd = static_cast<double>(m_sum_post);
      	const double sumdiffd = static_cast<double>(m_sum_pre - m_sum_post);
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

  void estimate_dobson(double& lci, double& uci) const
  {
    if constexpr(t_paired)
    {
      if(m_sum_post > m_sum_pre)
      {
        lci = NA_REAL;
        uci = NA_REAL;
      }
      else
      {
      	const double postsumd = static_cast<double>(m_sum_post);
      	const double sumdiffd = static_cast<double>(m_sum_pre - m_sum_post);
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


public:
  estimator(const size_t maxN, const double mean_ratio, const double H0_1, const double H0_2,
    const double tail, const std::array<double, 2L> conjugate_priors, const optswitch delta,
    const int beta_iters, const optswitch approx, const std::array<double, 2L> dobson_priors,
    const double true_effk_pre, const double true_effk_post) :
    m_maxN(maxN), m_mean_ratio(mean_ratio), m_H0_1(H0_1), m_H0_2(H0_2),
    m_tail(tail), m_lci(tail), m_uci(1.0 - tail), m_conjugate_priors(conjugate_priors), m_delta(delta),
    m_beta_iters(beta_iters), m_approx(approx), m_dobson_priors(dobson_priors),
    m_true_effk_pre(true_effk_pre), m_true_effk_post(true_effk_post)
  {
    // If k is fixed or mm then we don't need to store data:
    if constexpr(!s_store_count)
    {
      // Check that we are not getting all methods if fixing k:
      static_assert(!(t_all_methods && t_ktype == ktypes::fix), "Cannot fix k for all analysis methods");
    }
    else if constexpr(t_container==containers::stdarray)
    {
      // t_cont_type should be an array
      if(maxN >= m_pre.size() || maxN >= m_post.size())
      {
        Rcpp::stop("Specified maxN is larger than the array size");
      }
    }
    else if constexpr(t_container==containers::stdvector)
    {
      // t_cont_type should be a std::vector
      m_pre.resize(m_maxN);
      m_post.resize(m_maxN);
    }
    else if constexpr(t_container==containers::rcppvector)
    {
      // t_cont_type should be an Rcpp::IntegerVector
      Rcpp::IntegerVector pre(m_maxN);
      m_pre = pre;
      Rcpp::IntegerVector post(m_maxN);
      m_post = post;
    }
    else
    {
      Rcpp::stop("Unrecognised value for t_resizable");
    }

    if constexpr(t_all_methods)
    {
      // All methods:  Denwood, MLE, WAAVP, Levecke, Dobson
      static_assert(t_rvlen == (1L + 2L * 5L), "t_rvlen is an unexpected size for all methods");
    }
    else
    {
      static_assert(t_rvlen == (1L + 2L), "t_rvlen is an unexpected size for Denwood method");
    }

    reset();

  }

  void reset()
  {
    m_sum_pre = 0L;
    m_sum_post = 0L;
    m_mean_pre = 0.0;
    m_mean_post = 0.0;
    m_varnum_pre = 0.0;
    m_varnum_post = 0.0;
    m_covnum = 0.0;
    
    m_var_pre = NA_REAL;
    m_var_post = NA_REAL;
    m_cov = NA_REAL;
    m_effk_pre = NA_REAL;
    m_effk_post = NA_REAL;

    if constexpr(t_paired)
    {
      m_N = 0L;
    }
    else
    {
      m_Npre = 0L;
      m_Npost = 0L;
    }

  }

  void push_pre(const int obs)
  {
    if constexpr(t_paired)
    {
      Rcpp::stop("Attempt to push pre data alone to a paired analysis!");
    }
    add_pre(obs);
  }

  void push_post(const int obs)
  {
    if constexpr(t_paired)
    {
      Rcpp::stop("Attempt to push post data alone to a paired analysis!");
    }
    add_post(obs);
  }

  void push_pair(const int pre, const int post)
  {
    add_pair(pre, post);
  }

  // This will only be used from R so we can live with a deep copy:
  void push_multiple(const Rcpp::IntegerVector& pre, const Rcpp::IntegerVector& post)
  {
    if constexpr(t_paired)
    {
      if(pre.size() != post.size())
      {
        Rcpp::stop("Attempt to push unequal length pre- and post-data to a paired analysis!");
      }
      for(size_t i=0L; i<pre.size(); ++i)
      {
        add_pair(pre[i], post[i]);
      }
    }
    else
    {
      for(size_t i=0L; i<pre.size(); ++i)
      {
        add_pre(pre[i]);
      }
      for(size_t i=0L; i<post.size(); ++i)
      {
        add_post(post[i]);
      }
    }
  }

  void estimate(double& efficacy, double& p1, double& p2)
  {
    // Update internally stored variance:
    if constexpr(t_paired)
    {
      m_var_pre = m_varnum_pre / static_cast<double>(m_N-1L);
      m_var_post = m_varnum_post / static_cast<double>(m_N-1L);
      m_cov = m_covnum / static_cast<double>(m_N-1L);
    }
    else
    {
      m_var_pre = m_varnum_pre / static_cast<double>(m_Npre-1L);
      m_var_post = m_varnum_post / static_cast<double>(m_Npost-1L);
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

    int bnbfailed = bnb_pval(m_sum_pre, n_pre, m_effk_pre, m_mean_pre, m_var_pre, m_sum_post, n_post, m_effk_post, m_mean_post, m_var_post, m_cov,
                            m_mean_ratio, m_H0_1, m_H0_2, conjugate_priors, delta, m_beta_iters, approx, &p_1, &p_2);
    p1 = p_1;
    p2 = p_2;

  }

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
    rv[t_rvlen+2L] = m_var_pre;
    rv[t_rvlen+3L] = m_var_post;
    rv[t_rvlen+4L] = m_cov;

    return rv;
  }

};


#endif // ESTIMATOR_H
