#include "distribution.h"
#include "simulation.h"
#include "estimator.h"


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

// TODO: no need for this to be a class - a simple (non-templated) function wrapper is fine
class estimator_pair_unfix
{
private:
  static constexpr size_t s_rvlen = 9L;
  estimator<true, true, ktypes::ml, Rcpp::IntegerVector, containers::rcppvector, s_rvlen> m_estimator;
public:
  estimator_pair_unfix(const int maxN, const double k_pre, const double k_post) :
    m_estimator(maxN, k_pre, k_post)
	/*
	estimator(const size_t maxN, const double mean_ratio, const double H0_1, const double H0_2,
		const double lci, const double uci, const std::array<double, 2L> conjugate_priors, const optswitch delta,
		const int beta_iters, const optswitch approx, const std::array<double, 2L> dobson_priors,
		const double true_effk_pre, const double true_effk_post)
	 */
  {

  }

  void push_data(Rcpp::IntegerVector pre, Rcpp::IntegerVector post)
  {
    m_estimator.push_multiple(pre, post);
  }

  Rcpp::NumericVector estimate()
  {
    const std::array<double, s_rvlen> ests = m_estimator.estimate();
    Rcpp::NumericVector rv(s_rvlen);
    for(size_t i=0L; i<s_rvlen; ++i)
    {
      rv[i] = ests[i];
    }
    return rv;
  }
};

RCPP_EXPOSED_CLASS(estimator_pair_unfix)

RCPP_MODULE(bayescount_module){

	using namespace Rcpp;

  function("Rcpp_draw_lambda", &draw_lambda);
  function("Rcpp_draw_count", &draw_count);

	class_<estimator_pair_unfix>("Rcpp_estimator_pair_unfix")
		.constructor<const int, const double, const double>("Constructor")

    .method("push_data", &estimator_pair_unfix::push_data, "Push data")
    .method("estimate", &estimator_pair_unfix::estimate, "Produce estimates")

    /*
		.method("show", &Simulation::show, "The show method")
		.method("AddPatch", &Simulation::AddPatch, "Add a patch to the specified population")
		.method("Reset", &Simulation::Reset, "Reset all patches and interactions")
		.method("Update", &Simulation::Update, "Update with given number of time steps")

		.property("compartments", &Simulation::GetCompartments, "Get the total for each compartment")
		.property("states", &Simulation::GetStates, "Get the total for each state")
    */
	;


}
