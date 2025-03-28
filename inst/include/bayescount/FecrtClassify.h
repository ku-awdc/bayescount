/*
  FECRT analysis and classification class
  Takes data (pre/post or treatment and control) one obs at a time
  Returns classifications and/or CI / p-values

  TODO: finish me!
*/

enum class FecrtRT
{
  noninferiority,
  inferiority,
  both
};

enum class FecrtType
{
  paired,
  unpaired
};

enum class FecrtMethod
{
  delta,
  bnb,
  delta_bnb
};

template<FecrtRT s_fclass>
struct FecrtRV;

template<>
struct FecrtRV<FecrtRT::noninferiority>
{

};

template<FecrtRT s_fclass, FecrtType s_type, FecrtMethod s_method>
class FecrtClassify
{
private:
  
  // Non-linear transformation function:
  double g_fun(double p, double r, double s, double t){
  	return(p*r*t / (p*r*t - p*s + s));
  }

  // First 4 derivatives of g:
  double gp1_fun(double p, double r, double s, double t){
  	return(r*s*t / std::pow((p-1.0)*s - p*r*t, 2));
  }
  double gp2_fun(double p, double r, double s, double t){
  	return(-2.0*r*s*t*(s - r*t) / std::pow((p-1.0)*s - p*r*t, 3));
  }
  double gp3_fun(double p, double r, double s, double t){
  	return(6.0*r*s*t* std::pow(s - r*t, 2) / std::pow(p*r*t - p*s + s, 4));
  }
  double gp4_fun(double p, double r, double s, double t){
  	return(24.0*r*s*t* std::pow(r*t - s, 3) / std::pow((p-1.0)*s - p*r*t, 5));
  }

  // Functions for expectations of powers 2-5 of a Beta distribution with parameters a(lpha), b(eta):
  double epf2(double a, double b){
  	return(((a+1.0)*a)  / ((a+b+1.0)*(a+b)));
  }
  double epf3(double a, double b){
  	return(((a+2.0)*(a+1.0)*a)  / ((a+b+2.0)*(a+b+1.0)*(a+b)));
  }
  double epf4(double a, double b){
  	return(((a+3.0)*(a+2.0)*(a+1.0)*a)  / ((a+b+3.0)*(a+b+2.0)*(a+b+1.0)*(a+b)));
  }
  double epf5(double a, double b){
  	return(((a+4.0)*(a+3.0)*(a+2.0)*(a+1.0)*a)  / ((a+b+4.0)*(a+b+3.0)*(a+b+2.0)*(a+b+1.0)*(a+b)));
  }

  // The binomial expansions for E((x - mu)^t) for t in 2:5:
  double expow2(double a, double b){
  	double m = a / (a+b);
  	return(epf2(a,b) - std::pow(m,2));
  }
  double expow3(double a, double b){
  	double m = a / (a+b);
  	return(epf3(a,b) - 3.0*m*epf2(a,b) + 2*std::pow(m,3));
  }
  double expow4(double a, double b){
  	double m = a / (a+b);
  	return(epf4(a,b) - 4.0*m*epf3(a,b) + 6.0*std::pow(m,2)*epf2(a,b) - 3.0*std::pow(m,4));
  }
  double expow5(double a, double b){
  	double m = a / (a+b);
  	return(epf5(a,b) - 5.0*m*epf4(a,b) + 10.0*std::pow(m,2)*epf3(a,b) - 10.0*std::pow(m,3)*epf2(a,b) + 4.0*std::pow(m,5));
  }


  // Delta method approximation to the change in mean (first 2 terms of Taylor series):
  double delta_mean(double p_mu, double p_var, double r, double s, double t){

  	// p_mu and p_var are calculated by calling function (for efficiency) as:
  	// double p_mu = alpha / (alpha + beta);
  	// double p_var = (alpha * beta) / (std::pow(alpha + beta, 2) * (alpha + beta + 1));

  	double rm = g_fun(p_mu, r, s, t) + 0.5 * gp2_fun(p_mu, r, s, t) * p_var;
  	return(rm);
  }

  // Delta method approximation to the change in variance (higher order Taylor series):
  double delta_var(double alpha, double beta, double p_mu, double p_var, double r, double s, double t){

  	// p_mu and p_var are calculated by calling function (for efficiency) as:
  	// double p_mu = alpha / (alpha + beta);
  	// double p_var = (alpha * beta) / (std::pow(alpha + beta, 2) * (alpha + beta + 1));

  	double g1 = gp1_fun(p_mu, r, s, t);
  	double g2 = gp2_fun(p_mu, r, s, t);
  	double g3 = gp3_fun(p_mu, r, s, t);
  	double g4 = gp4_fun(p_mu, r, s, t);

  	double rv = std::pow(g1,2) * p_var +
  				2.0 * g1 * g2/2.0 * expow3(alpha, beta) +
  				(std::pow(g2,2)/4.0 + 2.0*g1 * g3/6.0) * expow4(alpha, beta) +
  				(2.0*g1 * g4/24.0 + g2 * g3/6.0) * expow5(alpha, beta);
  	return(rv);
  }
  


public:

  void input(int const x, int const y);

  FecrtRV<s_fclass> classify() const noexcept
  {
    FecrtRV<s_fclass> rv;
    
    
  }


};
