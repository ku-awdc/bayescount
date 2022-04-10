#ifndef SIMULATION_H
#define SIMULATION_H

// This is required by distribution.h and precludes Rcpp.h:
#include <RcppDist.h>
#include <array>

#include "enums.h"
#include "distribution.h"
#include "estimator.h"

// TODO: an estimater class that sets up the estimation method and potentially stores pre and post internally
// Then this gets passed to simulation, so that simulation doesn't need to store counts or do any estimating itself
// estimator class could have methods to add data 1 point at a time (running mean) or as an array or as a vector
// this means the newton rhapson function should be templated (or included as a method within estimator class??)

template<dists t_dist, size_t t_maxNarr, analyses t_analysis, size_t t_rvlen>
class simulation
{
private:
  static constexpr size_t s_dim = 2L;
  
  const Rcpp::IntegerVector& m_Ns;
  const std::array<const double, s_dim> m_mu;
  const std::array<const double, s_dim> m_cv;
  const double m_lcor;
  const size_t m_maxN;
  const distribution<t_dist> m_distribution;

  // Note: these are not actually needed if analyses::*_fixcv, but I could cheat and use 
  std::array<int, t_maxNarr> m_pre;
  std::array<int, t_maxNarr> m_post;

public:
  simulation(const Rcpp::IntegerVector& Ns, const std::array<const double, s_dim> mu, const std::array<const double, s_dim> cv, const double lcor) :
    m_Ns(Ns), m_mu(mu), m_cv(cv), m_lcor(lcor), m_maxN(max(m_Ns)), m_distribution(m_mu, m_cv, m_lcor)
  {
    if(m_maxN > t_maxNarr) Rcpp::stop("Invalid max Ns > maxN");
    // TODO: check Ns is sorted
    
    if(t_analysis == analyses::bnb_fixcv || t_analysis == analyses::bnb_estcv)
    {
      if(t_rvlen != 2L) Rcpp::stop("Invalid t_rvlen for analysis = bnb_*");
    }
    // TODO: check that rvlen matches that implied by analyses
  }

  std::array<double, t_rvlen> run(Rcpp::NumericVector& efficacy, Rcpp::NumericVector& pval_1, Rcpp::NumericVector& pval_2)
  {
    size_t ndone = 0L;
    for(size_t i=0L; i<m_maxN; ++i)
    {
      const std::array<const int, s_dim> rv = m_distribution.draw();
      m_pre[i] = rv[0L];
      m_post[i] = rv[1L];
      
      if(i==m_Ns[ndone])
      {
        
        ndone++;
      }
    }
    
    std::array<double, t_rvlen> rv;
    return rv;
  }

};


#endif // SIMULATION_H
