#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <array>

// Note: we can't include Rcpp as we need:
#include <RcppDist.h>

#include "enums.h"

template<dists dist>
class distribution;


template<>
class distribution<dists::lnormpois>
{
private:
  static constexpr size_t s_dim = 2L;
  const std::array<const double, s_dim> m_mu;
  const std::array<const double, s_dim> m_cv;
  std::array<double, s_dim> m_lmu;
  std::array<double, s_dim> m_lsd;

public:
  distribution(const std::array<const double, s_dim> mu, const std::array<const double, s_dim> cv, const double lcor) :
    m_mu(mu), m_cv(cv)
  {
    // Required for lcor to make sense:
    if(s_dim!=2L) Rcpp::stop("Invalid s_dim != 2L");
    // Note: s_dim is constexpr so this should get compiled away to nothing...

    if(abs(lcor) > 0.01) Rcpp::stop("Non-correlated distribution used where lcor != 0");

    for(size_t i=0L; i<s_dim; ++i)
    {
      const double lvar = log(pow(m_cv[i], 2.0) + 1.0);
      m_lmu[i] = log(m_mu[i]) - (lvar / 2.0);
      m_lsd[i] = sqrt(lvar);
    }
  }

  std::array<double, s_dim> draw_mu() const
  {
    std::array<double, s_dim> rv;
    for(size_t i=0L; i<s_dim; ++i)
    {
      rv[i] = R::rlnorm(m_lmu[i], m_lsd[i]);
    }
  	return rv;
  }

  std::array<int, s_dim> draw() const
  {
    std::array<int, s_dim> rv;
    for(size_t i=0L; i<s_dim; ++i)
    {
      const double lambda = R::rlnorm(m_lmu[i], m_lsd[i]);
      rv[i] = R::rpois(lambda);
    }
  	return rv;
  }

};


template<>
class distribution<dists::mvlnormpois>
{
private:
  static constexpr size_t s_dim = 2L;
  const std::array<const double, s_dim> m_mu;
  const std::array<const double, s_dim> m_cv;
  const double m_lcor;
  arma::vec m_lmu;
  arma::mat m_lsigma;

public:
  distribution(const std::array<const double, s_dim> mu, const std::array<const double, s_dim> cv, const double lcor) :
    m_mu(mu), m_cv(cv), m_lcor(lcor)
  {
    // Required for lcor to make sense:
    if(s_dim!=2L) Rcpp::stop("Invalid s_dim != 2L");
    // Note: s_dim is constexpr so this should get compiled away to nothing...

    if(abs(lcor) > 1.0) Rcpp::stop("Invalid abs(lcor) > 1.0");

    Rcpp::NumericVector lmu(s_dim);
    Rcpp::NumericMatrix lsigma(s_dim, s_dim);

    for(size_t i=0L; i<s_dim; ++i)
    {
      lsigma(i,i) = log(pow(m_cv[i], 2.0) + 1.0);
      lmu[i] = log(m_mu[i]) - (lsigma(i,i) / 2.0);
    }

    const double lcov = lcor * sqrt(lsigma(0L,0L)) * sqrt(lsigma(1L,1L));
    lsigma(0L,1L) = lcov;
    lsigma(1L,0L) = lcov;

    m_lmu = Rcpp::as<arma::vec>(lmu);
    m_lsigma = Rcpp::as<arma::mat>(lsigma);

  }

  std::array<double, s_dim> draw_mu() const
  {

    const arma::mat rvarr = rmvnorm(1L, m_lmu, m_lsigma);

    std::array<double, s_dim> rv;
    for(size_t i=0L; i<s_dim; ++i)
    {
      rv[i] = exp(rvarr[0L,i]);
    }
  	return rv;
  }

  std::array<int, s_dim> draw() const
  {

    std::array<double, s_dim> mus = draw_mu();

    std::array<int, s_dim> rv;
    for(size_t i=0L; i<s_dim; ++i)
    {
      rv[i] = R::rpois(mus[i]);
    }
  	return rv;
  }

};


template<>
class distribution<dists::nbinom>
{
private:
  static constexpr size_t s_dim = 2L;
  const std::array<const double, s_dim> m_mu;
  const std::array<const double, s_dim> m_cv;
  const std::array<const double, s_dim> m_k;
  const std::array<const double, s_dim> m_rate;
public:
  // Note: this assumes that s_dim==2L
  distribution(const std::array<const double, s_dim> mu, const std::array<const double, s_dim> cv, const double lcor) :
    m_mu(mu), m_cv(cv), m_k { pow(m_cv[0L], -2.0), pow(m_cv[1L], -2.0) }, m_rate { m_mu[0L]/m_k[0L], m_mu[1L]/m_k[1L] }
  {
    if(abs(lcor) > 0.01) Rcpp::stop("NB distribution used where lcor != 0");
  }

  std::array<double, s_dim> draw_mu() const
  {
    std::array<double, s_dim> rv;
    for(size_t i=0L; i<s_dim; ++i)
    {
      if(m_cv[i] <= 0.0)
      {
        rv[i] = m_mu[i];
      }
      else
      {
        rv[i] = R::rgamma(m_k[i], m_rate[i]);;
      }
    }
  	return rv;
  }

  std::array<int, s_dim> draw() const
  {
    std::array<int, s_dim> rv;
    for(size_t i=0L; i<s_dim; ++i)
    {
      if(m_cv[i] <= 0.0)
      {
        rv[i] = R::rpois(m_mu[i]);
      }
      else
      {
        rv[i] = rnbinom_mu(m_k[i], m_mu[i]);
      }
    }
  	return rv;
  }

};


template<>
class distribution<dists::mvnbinom>
{
private:
  static constexpr size_t s_dim = 2L;
  const std::array<const double, s_dim> m_mu;
  const std::array<const double, s_dim> m_cv;
  const double m_lcor;
  const std::array<const double, s_dim> m_k;
  const std::array<const double, s_dim> m_rate;
  arma::vec m_lmu;
  arma::mat m_lsigma;
  std::array<double, s_dim> m_lsd;

public:
  distribution(const std::array<const double, s_dim> mu, const std::array<const double, s_dim> cv, const double lcor) :
    m_mu(mu), m_cv(cv), m_lcor(lcor), m_k { pow(m_cv[0L], -2.0), pow(m_cv[1L], -2.0) }, m_rate { m_mu[0L]/m_k[0L], m_mu[1L]/m_k[1L] }
  {
    // Required for lcor to make sense:
    if(s_dim!=2L) Rcpp::stop("Invalid s_dim != 2L");
    // Note: s_dim is constexpr so this should get compiled away to nothing...

    if(abs(lcor) > 1.0) Rcpp::stop("Invalid abs(lcor) > 1.0");

    Rcpp::NumericVector lmu(s_dim);
    Rcpp::NumericMatrix lsigma(s_dim, s_dim);

    for(size_t i=0L; i<s_dim; ++i)
    {
      lsigma(i,i) = log(pow(m_cv[i], 2.0) + 1.0);
      lmu[i] = log(m_mu[i]) - (lsigma(i,i) / 2.0);
      m_lsd[i] = sqrt(lsigma(i,i));
    }

    const double lcov = lcor * sqrt(lsigma(0L,0L)) * sqrt(lsigma(1L,1L));
    lsigma(0L,1L) = lcov;
    lsigma(1L,0L) = lcov;

    m_lmu = Rcpp::as<arma::vec>(lmu);
    m_lsigma = Rcpp::as<arma::mat>(lsigma);

  }

  std::array<double, s_dim> draw_mu() const
  {

    const arma::mat rvarr = rmvnorm(1L, m_lmu, m_lsigma);

    std::array<double, s_dim> rv;
    for(size_t i=0L; i<s_dim; ++i)
    {
		  // Signature eg:  double pgamma (double x, double alph, double scale, int lower_tail, int log_p)
      const double pv = R::plnorm(exp(rvarr[0L,i]), m_lmu[i], m_lsd[i], true, true);
      rv[i] = R::qgamma(pv, m_k[i], m_rate[i], true, true);
    }
  	return rv;
  }

  std::array<int, s_dim> draw() const
  {
    
    const arma::mat rvarr = rmvnorm(1L, m_lmu, m_lsigma);

    std::array<int, s_dim> rv;
    for(size_t i=0L; i<s_dim; ++i)
    {
		  // Signature eg:  double pgamma (double x, double alph, double scale, int lower_tail, int log_p)
      const double pv = R::plnorm(exp(rvarr[0L,i]), m_lmu[i], m_lsd[i], true, true);
      const double lambda = R::qgamma(pv, m_k[i], m_rate[i], true, true);
      rv[i] = R::rpois(lambda);
    }
  	return rv;
  }

};


template<>
class distribution<dists::poisson>
{
private:
  static constexpr size_t s_dim = 2L;
  const std::array<const double, s_dim> m_mu;
public:
  distribution(const std::array<const double, s_dim> mu, const std::array<const double, s_dim> cv, const double lcor) :
    m_mu(mu)
  {
    if(abs(lcor) > 0.01) Rcpp::stop("NB distribution used where lcor != 0");
    for(size_t i=0L; i<s_dim; ++i)
    {
      if(cv[i] > 0.0) Rcpp::stop("Poisson distribution used where cv > 0.0");
    }
  }

  std::array<double, s_dim> draw_mu() const
  {
    std::array<double, s_dim> rv;
    for(size_t i=0L; i<s_dim; ++i)
    {
      rv[i] = m_mu[i];
    }
  	return rv;
  }

  std::array<int, s_dim> draw() const
  {
    std::array<int, s_dim> rv;
    for(size_t i=0L; i<s_dim; ++i)
    {
      rv[i] = R::rpois(m_mu[i]);
    }
  	return rv;
  }

};


template<dists dist>
Rcpp::NumericMatrix draw_lambda_template(const int n, const Rcpp::NumericVector& mu, const Rcpp::NumericVector& cv, const double lcor)
{
  constexpr size_t s_dim = 2L;

  if(mu.size() != s_dim) Rcpp::stop("Invalid mu size");
  if(cv.size() != s_dim) Rcpp::stop("Invalid cv size");
  if(n <= 0L) Rcpp::stop("Invalid n");

  const std::array<const double, s_dim> muarr = { mu[0L], mu[1L] };
  const std::array<const double, s_dim> cvarr = { cv[0L], cv[1L] };

  Rcpp::NumericMatrix rv(n, s_dim);
  distribution<dist> dist_obj(muarr, cvarr, lcor);

  for(size_t i=0L; i<n; ++i)
  {
    const std::array<double, s_dim> tt = dist_obj.draw_mu();
    for(size_t j=0L; j<s_dim; j++)
    {
      rv(i,j) = tt[j];
    }

  }

  return rv;
}


template<dists dist>
Rcpp::IntegerMatrix draw_count_template(const int n, const Rcpp::NumericVector& mu, const Rcpp::NumericVector& cv, const double lcor)
{
  constexpr size_t s_dim = 2L;

  if(mu.size() != s_dim) Rcpp::stop("Invalid mu size");
  if(cv.size() != s_dim) Rcpp::stop("Invalid cv size");
  if(n <= 0L) Rcpp::stop("Invalid n");

  const std::array<const double, s_dim> muarr = { mu[0L], mu[1L] };
  const std::array<const double, s_dim> cvarr = { cv[0L], cv[1L] };

  Rcpp::IntegerMatrix rv(n, s_dim);
  distribution<dist> dist_obj(muarr, cvarr, lcor);

  for(size_t i=0L; i<n; ++i)
  {
    const std::array<int, s_dim> tt = dist_obj.draw();
    for(size_t j=0L; j<s_dim; j++)
    {
      rv(i,j) = tt[j];
    }

  }

  return rv;
}

#endif // DISTRIBUTION_H
