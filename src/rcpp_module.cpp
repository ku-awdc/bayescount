// This needs to be RcppDist due to depdendencies on distribution.h
#include <RcppDist.h>

#include "rcpp_wrappers.h"

/*
class estimator_pair_unfix
{
private:
  static constexpr size_t s_rvlen = 9L;
  estimator<true, true, ktypes::ml, Rcpp::IntegerVector, containers::rcppvector, s_rvlen> m_estimator;
public:
  estimator_pair_unfix(const int maxN, const double k_pre, const double k_post) :
    m_estimator(maxN, k_pre, k_post)

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
*/

#include "../inst/include/bayescount/bnb_pval.h"


RCPP_MODULE(bayescount_module){

	using namespace Rcpp;

  function("Rcpp_draw_lambda", &draw_lambda);
  function("Rcpp_draw_count", &draw_count);
  function("Rcpp_estimate_fecrt", &estimate_fecrt);
  function("Rcpp_summarise_fecrt", &summarise_fecrt);
  function("Rcpp_bnb_pval_100", &bayescount::bnb_pval_100);

    /*
	class_<estimator_pair_unfix>("Rcpp_estimator_pair_unfix")
		.constructor<const int, const double, const double>("Constructor")

    .method("push_data", &estimator_pair_unfix::push_data, "Push data")
    .method("estimate", &estimator_pair_unfix::estimate, "Produce estimates")

		.method("show", &Simulation::show, "The show method")
		.method("AddPatch", &Simulation::AddPatch, "Add a patch to the specified population")
		.method("Reset", &Simulation::Reset, "Reset all patches and interactions")
		.method("Update", &Simulation::Update, "Update with given number of time steps")

		.property("compartments", &Simulation::GetCompartments, "Get the total for each compartment")
		.property("states", &Simulation::GetStates, "Get the total for each state")
	;
    */


}
