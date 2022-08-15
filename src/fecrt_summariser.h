#ifndef FECRT_SUMMARISER_H
#define FECRT_SUMMARISER_H

// This is required by distribution.h and precludes Rcpp.h:
#include <RcppDist.h>
#include <array>
#include <type_traits>
#include <typeinfo>

#include "enums.h"

// fecrt_summarise does the count storage and mean/var/cov/k estimation
// It needs to know about paired, ktypes and containers but not distributions
// It also needs to know the true k, but only uses it for fixk ktype

template<bool t_paired, ktypes t_ktype, typename t_cont_type, containers t_container>
class fecrt_summariser
{
private:
  const size_t m_maxN;
	const std::string m_ktype = "templated";

  // Note: not actually needed unless t_ktype is fix
  // Also: should be pre-calculated as effective k for paired:
  const double m_true_effk_pre = 0.0;
  const double m_true_effk_post = 0.0;

  // Note: these are not actually needed if !s_store_count (defined below):
  t_cont_type m_pre;
  t_cont_type m_post;

  // TODO: get rid of the sums as we can cast back from mean*N (overflow safer?)
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
  // TODO: std::optional for these? https://en.cppreference.com/w/cpp/utility/optional
  type_paired m_N = 0L;
  type_unpaired m_Npre = 0L;
  type_unpaired m_Npost = 0L;

  static constexpr bool s_store_count = (t_ktype==ktypes::flex || t_ktype==ktypes::ql || t_ktype==ktypes::ml);
  static constexpr size_t s_summary_len = t_paired ? 12L : 8L;
  // double(Ns), mus, vars, ks (adjusted), if(t_paired): cov, cor, ks (unadjusted)

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

    const double dobs = static_cast<double>(obs);
		const double delta = obs - m_mean_pre;
    m_mean_pre += delta / static_cast<double>(m_Npre);
		m_varnum_pre += delta * (dobs - m_mean_pre);
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

    const double dobs = static_cast<double>(obs);
    const double delta = obs - m_mean_post;
    m_mean_post += delta / static_cast<double>(m_Npost);
		m_varnum_post += delta * (dobs - m_mean_post);
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

      const double dN = static_cast<double>(m_N);
      const double dpre = static_cast<double>(pre);
      const double dpost = static_cast<double>(post);
      const double delta1 = dpre - m_mean_pre;
  		const double delta2 = dpost - m_mean_post;

      m_mean_pre += delta1 / dN;
  		m_varnum_pre += delta1 * (dpre - m_mean_pre);
      // Variance = m_varnum_pre / (m_N-1.0)

      m_mean_post += delta2 / dN;
  		m_varnum_post += delta2 * (dpost - m_mean_post);
      // Variance = m_varnum_post / (m_N-1.0)

      const double iip1 = (dN-1.0) / dN;
      // const double iip1 = static_cast<double>(m_N-1L) / static_cast<double>(m_N);
      m_covnum += delta1 * delta2 * iip1;
      // Coariance = m_covnum / (m_N-1.0)

      m_sum_pre += pre;
      m_sum_post += post;

    }
  }

  //  This method is based on a function extracted from R-3.6.1/src/library/stats/src/optimize.c
  double find_theta(const t_cont_type& data, const size_t N, const double mu, const double min=0.001, const double max=20.0, const double tol=0.01)
  {
  	// Note: N is needed as an argument as we may not be using the whole array

  	const double ax = min;
  	const double bx = max;

  	/*  c is the squared inverse of the golden ratio */
  	const double c = (3. - sqrt(5.)) * .5;

  	/* Local variables */
  	double a, b, d, e, p, q, r, u, v, w, x;
  	double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

  	/*  eps is approximately the square root of the relative machine precision. */
  	eps = DBL_EPSILON;
  	tol1 = eps + 1.;/* the smallest 1.000... > 1 */
  	eps = std::sqrt(eps);

  	a = ax;
  	b = bx;
  	v = a + c * (b - a);
  	w = v;
  	x = v;

  	d = 0.;/* -Wall */
  	e = 0.;
  	//    fx = (*f)(x, info);
  	{}
  	fx = 0.0;
  	for (size_t i=0L; i<N; ++i)
  	{
  		fx -= dnbinom_mu(data[i], x, mu, true);
  	}
  	/*
  	Rcpp::Rcout << fx;
  	{
  		Rcpp::NumericVector probs = Rcpp::dnbinom_mu(data, x, mu, true);
  		fx = - Rcpp::sum(probs);
  	}
  	Rcpp::Rcout << " - " << fx << std::endl;
  	*/

  	fv = fx;
  	fw = fx;
  	tol3 = tol / 3.;

  	/*  main loop starts here ----------------------------------- */

  	for(;;) {

  		xm = (a + b) * .5;
  		tol1 = eps * fabs(x) + tol3;
  		t2 = tol1 * 2.;

  		/* check stopping criterion */

  		if (fabs(x - xm) <= t2 - (b - a) * .5) break;
  		p = 0.;
  		q = 0.;
  		r = 0.;
  		if (fabs(e) > tol1) { /* fit parabola */

  		r = (x - w) * (fx - fv);
  			q = (x - v) * (fx - fw);
  			p = (x - v) * q - (x - w) * r;
  			q = (q - r) * 2.;
  			if (q > 0.) p = -p; else q = -q;
  			r = e;
  			e = d;
  		}

  		if (fabs(p) >= fabs(q * .5 * r) ||
        p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

  		if (x < xm) e = b - x; else e = a - x;
  		d = c * e;
  		}
  		else { /* a parabolic-interpolation step */

  		d = p / q;
  			u = x + d;

  			/* f must not be evaluated too close to ax or bx */

  			if (u - a < t2 || b - u < t2) {
  				d = tol1;
  				if (x >= xm) d = -d;
  			}
  		}

  		/* f must not be evaluated too close to x */

  		if (fabs(d) >= tol1)
  			u = x + d;
  		else if (d > 0.)
  			u = x + tol1;
  		else
  			u = x - tol1;

  		// fu = (*f)(u, info);
  		fu = 0.0;
  		for (size_t i=0L; i<N; ++i)
  		{
  			fu -= dnbinom_mu(data[i], u, mu, true);
  		}
  		/*
  		{
  			Rcpp::NumericVector probs = Rcpp::dnbinom_mu(data, u, mu, true);
  			fu = - Rcpp::sum(probs);
  		}
  		*/

  		/*  update  a, b, v, w, and x */

  		if (fu <= fx) {
  			if (u < x) b = x; else a = x;
  			v = w;    w = x;   x = u;
  			fv = fw; fw = fx; fx = fu;
  		} else {
  			if (u < x) a = u; else b = u;
  			if (fu <= fw || w == x) {
  				v = w; fv = fw;
  				w = u; fw = fu;
  			} else if (fu <= fv || v == x || v == w) {
  				v = u; fv = fu;
  			}
  		}
  	}
  	/* end of main loop */

  	return x;
  }


  void setup()
  {
  	// If k is fixed or mm then we don't need to store data:
  	if constexpr(!s_store_count)
  	{
  		// Check that we are not getting all methods if fixing k:
  		static_assert(t_ktype == ktypes::fix, "ktype is not fixed for not storing counts");
  	}
  	else if constexpr(t_container==containers::stdarray)
  	{
  		// t_cont_type should be an array
  		if(m_maxN >= m_pre.size() || m_maxN >= m_post.size())
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

  	reset();

  }


public:
	// Probably not going to be used:
	fecrt_summariser(const size_t maxN, const std::string& k_type, const double true_effk_pre, const double true_effk_post) :
    m_maxN(maxN), m_ktype(k_type), m_true_effk_pre(true_effk_pre), m_true_effk_post(true_effk_post)
  {
		setup();
  }

	// Used for interface from R:
	fecrt_summariser(const size_t maxN, const std::string& k_type) :
		m_maxN(maxN), m_ktype(k_type)
	{
		setup();
	}

	// Used for interfact from C++ for fixk:
	fecrt_summariser(const size_t maxN, const double true_effk_pre, const double true_effk_post) :
		m_maxN(maxN)
	{
		setup();
	}

	// Used for interfact from C++ for not fixk:
	fecrt_summariser(const size_t maxN) :
		m_maxN(maxN)
	{
		setup();
	}

	// TODO: remove reset function
	void reset()
  {
		/*
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
		 */

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
  Rcpp::NumericVector push_multiple(const Rcpp::IntegerVector& pre, const Rcpp::IntegerVector& post)
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

    std::array<double, s_summary_len> rva = summarise();
    Rcpp::NumericVector rv(s_summary_len);
    for (size_t i=0L; i<s_summary_len; ++i)
   	{
    	rv[i] = rva[i];
   	}

		return rv;
  }

	std::array<double, s_summary_len> summarise()
	{

		// Return values:
		// double(Ns), mus, vars, ks (adjusted), if(t_paired): cov, cor, ks (unadjusted)

		std::array<double, s_summary_len> rv;

		// Get Ns:
		if constexpr(t_paired)
		{
			rv[0L] = static_cast<double>(m_N);
			rv[1L] = static_cast<double>(m_N);
		}
		else
		{
			rv[0L] = static_cast<double>(m_Npre);
			rv[1L] = static_cast<double>(m_Npost);
		}

		// Get means:
		rv[2L] = m_mean_pre;
		rv[3L] = m_mean_post;

		// Get variances:
		rv[4L] = m_varnum_pre / (rv[0L] - 1.0);
		rv[5L] = m_varnum_post / (rv[1L] - 1.0);

		// Estimate k if necessary:
		if constexpr(t_ktype == ktypes::fix)
		{
			rv[6L] = m_true_effk_pre;
			rv[7L] = m_true_effk_post;
		}
		else if constexpr(t_ktype == ktypes::ml)
		{
			rv[6L] = find_theta(m_pre, t_paired ? m_N : m_Npre, m_mean_pre);
			rv[7L] = find_theta(m_post, t_paired ? m_N : m_Npost, m_mean_post);
		}
		else if constexpr(t_ktype == ktypes::flex)
		{
			if (m_ktype == "fix")
			{
				rv[6L] = m_true_effk_pre;
				rv[7L] = m_true_effk_post;
			}
			else if (m_ktype == "ml")
			{
				rv[6L] = find_theta(m_pre, t_paired ? m_N : m_Npre, m_mean_pre);
				rv[7L] = find_theta(m_post, t_paired ? m_N : m_Npost, m_mean_post);
			}
			else
			{
				Rcpp::stop("Unhandled m_ktype in summarise()");
			}
		}
		else
		{
			// TODO: add ql and ml
			Rcpp::stop("Unhandled t_ktype in summarise()");
		}

		// Adjust k for paired:
		if constexpr(t_paired)
		{
			// Copy unadjusted values
		}

		// Get covariance and adjust k:
		if constexpr(t_paired)
		{
			rv[8L] = m_covnum / (rv[0L] - 1.0);
			rv[9L] = rv[8L] / (std::sqrt(rv[4L]) * std::sqrt(rv[5L]));

			// Copy unadjusted ks:
			rv[10L] = rv[6l];
			rv[11L] = rv[7l];

			// Adjust ks:
			rv[6L] = NA_REAL;
			rv[7L] = NA_REAL;
			// TODO: adjust properly

		}

		return rv;
	}

	// TODO: template this for ktype and s_summary_len???
	std::array<double, 6L> summarise_adjusted()
	{
		std::array<double, s_summary_len> sv = summarise();

		// Return values are always: mean1, mean2, var1, var2, adj_k1, adj_k2
		std::array<double, 6L> rv;
		for(size_t i=0L; i<6L; ++i)
		{
			rv[i] = sv[i];
		}

		bool adjust_k = false;
		if constexpr(t_paired)
		{
			if constexpr(t_ktype == ktypes::fix)
			{
				// Do nothing
			}
			else if constexpr(t_ktype == ktypes::flex)
			{
				// Runtime check is needed
				if (m_ktype != "fix")
				{
					adjust_k = true;
				}
			}
			else
			{
				adjust_k = true;
			}
		}

		if (adjust_k)
		{
			// TODO: adjust k estimates
		}

		return rv;
	}

};


#endif // FECRT_SUMMARISER_H
