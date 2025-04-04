#ifndef BAYESCOUNT_WAAVP_H
#define BAYESCOUNT_WAAVP_H

#include <array>

namespace bayescount
{
  
  inline std::array<double, 2> waavp_p_ci(double mu1, double mu2, double var1, double var2, double cov12, int N, double tail){

  	// Method B of Lyndal-Murphy, M., Swain, a J., & Pepper, P. M. (2014). Methods to determine resistance to anthelmintics when continuing larval development occurs. Veterinary Parasitology, 199(3–4), 191–200. https://doi.org/10.1016/j.vetpar.2013.11.002

  	int df = N -1;
  	// Signature of qt is:  double  qt(double, double, int, int);
  	double tval = R::qt(1.0 - tail, (double) df, 1, 0);

  	double Nd = (double) N;
  	double varred = var1 / (Nd * mu1 * mu1) + var2 / (Nd * mu2 * mu2) - 2.0 * cov12 / (Nd * mu1 * mu2);

    std::array<double, 2> rv;
  	rv[1] = 1 - (mu2 / mu1 * std::exp(-tval * std::sqrt(varred)));
  	rv[0] = 1 - (mu2 / mu1 * std::exp(tval * std::sqrt(varred)));
    return rv;

  }


  inline std::array<double, 2> waavp_u_ci(double mu1, double mu2, double var1, double var2, int N1, int N2, double tail){

  	// Method A of Lyndal-Murphy, M., Swain, a J., & Pepper, P. M. (2014). Methods to determine resistance to anthelmintics when continuing larval development occurs. Veterinary Parasitology, 199(3–4), 191–200. https://doi.org/10.1016/j.vetpar.2013.11.002

  	int df = N1 + N2 - 2;
  	// Signature of qt is:  double  qt(double, double, int, int);
  	double tval = R::qt(1.0 - tail, (double) df, 1, 0);

  	double N1d = (double) N1;
  	double N2d = (double) N2;
  	double varred = var1 / (N1d * mu1 * mu1) + var2 / (N2d * mu2 * mu2);

    std::array<double, 2> rv;
  	rv[1] = 1 - (mu2 / mu1 * std::exp(-tval * std::sqrt(varred)));
  	rv[0] = 1 - (mu2 / mu1 * std::exp(tval * std::sqrt(varred)));
    return rv;

  }

} // namespace bayescount

#endif // BAYESCOUNT_WAAVP_H
