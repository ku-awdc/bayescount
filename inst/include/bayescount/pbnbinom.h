/*
General class to analyse FECRT/ERR data using a number of methods, based on input summary statistics

Underlying code within the anonymous namespace (hypergeometric distributions) was adapted from SuppDists version 1.1-9
Original file is copyright of Bob Wheeler, licensed GPL>=2
*/

#ifndef BAYESCOUNT__PBNBINOM_H
#define BAYESCOUNT__PBNBINOM_H

#include <Rcpp.h>


namespace bayescount
{
  namespace {
    /*
    Underlying code within this anonymous namespace (hypergeometric distributions) was adapted
    from SuppDists version 1.1-9
    Original file is copyright of Bob Wheeler, licensed GPL>=2
    */


    enum hyperType {
    	classic,
    	IAi,
    	IAii,
    	IB,
    	IIA,
    	IIB,
    	IIIA,
    	IIIB,
    	IV,
    	noType
    };

    inline constexpr double LOG10 = 2.3025850929940456840179915;
    inline constexpr double MAXEXP = LOG10*DBL_MAX_10_EXP;	// Maximum argument for exp()
    inline constexpr double LOGSQRT2PI = 0.9189385332046727417803296;

    template<typename T>
    T maxm(T a, T b)
    {
      T rv = (((a)>(b))?(a):(b));
      return rv;
    }

    template<typename T>
    T minm(T a, T b)
    {
      T rv = (((a)<(b))?(a):(b));
      return rv;
    }

    template<typename T>
    T absm(T a)
    {
      T rv = (((a)<0)?-(a):(a));
      return rv;
    }

    template<typename T>
    T signm(T a)
    {
      T rv = (((a)>0)?1:(((a)<0)?-1:0));
      return rv;
    }

    template<typename T>
    T SQR(T a)
    {
      T rv = ((a)*(a));
      return rv;
    }


    /*
    	Natural logarithm of the gamma function.
    	     CACM 291 due to M.C. Pike and I.D. Hill
    	     Accurate to at least 10 places.
    */
    inline double loggamma(double x)	 {
    	const double  T1=1.0/12.0;
    	const double  T3=1.0/360.0;
    	const double  T5=1.0/1260.0;
    	const double  T7=1.0/1680.0;
    	const double  T9=1.0/1188.0;

    	double f;

    	if (x == 1.0 || x == 2.0) {
    		return 0.0;
    	}

    	if (x>=7.0) {
    		f=0.0;
    	}
    	else {
    		for (f=1.0;x<7.0;x+=1.0) {
    			f*=x;
    		}
    		f=-log(f);
    	}

    	double z=1.0/(x*x);
    	f+=(x-0.5)*log(x)-x+LOGSQRT2PI;
    	double t=(T9*z-T7)*z+T5;
    	f+=((t*z-T3)*z+T1)/x;

    	return (f);
    }




     	// Normal approximation to the hypergeometric distribution function due to Peizer
    	// See Ling, R.F. and Pratt, J.W. (1984) The accuracy of Peizer approximations
    	//  to the hypergeometric distribution, with comparisons to some other
    	//  approximations. JASA 79-385. 49-60.
    inline double PeizerHypergeometric(
    	int x,	  // Number of marked items in sample
    	int S,	  // Sample size
    	int n,	  // Total number of marked items
    	int N	  // Total number of items
    )
    {
    	const double oneSix=1.0/6.0;

    	double dn=(double)n;
    	double dm=(double)(N-n);
    	double dr=(double)S;
    	double ds=(double)(N-S);
    	double dN=(double)N;
    	double dnp=dn+oneSix;
    	double dmp=dm+oneSix;
    	double drp=dr+oneSix;
    	double dsp=ds+oneSix;
    	double dNp=dN-oneSix;
    	double A=(double)x+0.5;
    	double B=maxm(dn-A,0.5); // prevents B or C from going neg when x=n or x=S
    	double C=maxm(dr-A,0.5);
    	double D=(dm-dr)+A;
    	double Ap=A+oneSix+0.02/(A+0.5)+0.01/(dn+1.0)+0.01/(dr+1.0);
    	double Bp=B-oneSix+0.02/(B+0.5)+0.01/(dn+1.0)+0.01/(ds+1.0);
    	double Cp=C-oneSix+0.02/(C+0.5)+0.01/(dm+1.0)+0.01/(dr+1.0);
    	double Dp=D+oneSix+0.02/(D+0.5)+0.01/(dm+1.0)+0.01/(ds+1.0);

    	double L=A*log((A*dN)/(dn*dr))+B*log((B*dN)/(dn*ds))+C*log((C*dN)/(dm*dr))+D*log((D*dN)/(dm*ds));

    	double z=((Ap*Dp-Bp*Cp)/fabs(A*D-B*C))*sqrt(2.0*L*((dm*dn*dr*ds*dNp)/(dmp*dnp*drp*dsp*dN)));

    	return R::pnorm(z,0,1,true,false);

    }


    /*
    	The Gaussian hypergeometric function, usually denoted as
    	 F[a,b;c;x]
    */
    inline double GaussianHypergometricFcn(
    	double a,
    	double b,
    	double c,
    	double x
    )
    {
    	int const MAXITERATES=100;

    	if (c<0.0 && floor(c) == c)
    		return NA_REAL;

    	double sum=0.0;
    	double term=1.0;
    	int j=1;
    	double dj;
    	double djm1;
    	do {
    		dj=(double)j;
    		djm1=dj-1.0;
    		sum+=term;
    		term*=((a+djm1)*(b+djm1))/(c+djm1)*(x/dj);
    		j++;
    	} while(!((sum+term) == sum || j>MAXITERATES));

    	return sum;
    }


    	// Returns true if the double is an int
    inline bool isint(
    	double x
    )
    {
    	return x == floor(x);
    }



    inline bool checkHyperArgument(
    	int k,
    	double a, 				// Sample size
    	double m,      	// Total number of marked items
    	double N,       		// Total number of items
    	hyperType variety
    )
    {

    	switch (variety) {
    		case classic:
    			return (maxm(0,(int)(a+m-N))<=k && k<=minm((int)a,(int)m));
    			break;
    		case IAi:
    			return (0<=k && k<=(int)m);
    			break;
    		case IAii:
    			return (0<=k && k<=(int)a);
    			break;
    		case IB:		// Specified 1.0<N to avoid problems with small parameters
    			return (0<=k);
    			break;
    		case IIA:
    			return (0<=k && k<=(int)m);
    			break;
    		case IIB:
    			return (0<=k);
    			break;
    		case IIIA:
    			return (0<=k && k<=(int)a);
    			break;
    		case IIIB:
    			return (0<=k);
    			break;
    		case IV:
    			return (0<=k);
    			break;
    		case noType:
    			break;
    	}
    	return false;
    }

    inline double phypergeometric(
    	int x,		// Number of marked items in sample
    	int a, 				// Sample size
    	int n,      	// Total number of marked items
    	int N       		// Total number of items
    )
    {
    	if (x<maxm(0,a-(N-n)) || x>minm(a,n))
    		return NA_REAL;

    		// interchange n and a to get the fewest terms to sum
    	if (a<n) {
    		int k=a;
    		a=n;
    		n=k;
    	}

    	if (x == n) {
    		return 1.0;
    	}

    		// Switch tails if necessesary to minimize number of terms to sum
    	int xmin=maxm(0,n+a-N);
    	bool lowerTail=true;
    	if (x-xmin>n-x) {
    		x=n-x-1;
    		a=N-a;
    		xmin=maxm(0,n+a-N);
    		lowerTail=false;
    	}

    	int na_N=n+a-N;
    	double logP=loggamma((double)(a+1))+loggamma((double)(N-a+1))+loggamma((double)(n+1))+
    		loggamma((double)(N-n+1))-loggamma((double)(N+1))-loggamma((double)(a-xmin+1))-
    		loggamma((double)(n-xmin+1))-loggamma(xmin-na_N+1);

    	if (xmin!=0) {
    		logP-=loggamma((double)(xmin+1));
    	}

    		// Use normal approximation if can't do it
    	if (! R_FINITE(logP)){
    		double p=PeizerHypergeometric(x,a,n,N);
    		return lowerTail?p:1.0-p;
    	}

    	double term=1.0;
    	double sum=1.0;
    		// These are the terms of F[-a,-n;N-n-a+1;a], where F is the Gaussian
    		//  hypergeometric function -- i.e. coefficients of x^i in the expansion.
    	for (int k=xmin;k<x;k++) {
    		term*=((double)(a-k)*(double)(n-k))/((double)(k+1)*(double)(k+1-na_N));
    		sum+=term;
    	}
    		// Use normal aapproximation if can't do it
    	if (! R_FINITE(sum)){
    		double p=PeizerHypergeometric(x,a,n,N);
    		return lowerTail?p:1.0-p;
    	}

    	logP+=log(sum);
    	if (logP<-MAXEXP) {
    		return lowerTail?0.0:1.0;
    	}
    	else {
    		return lowerTail?exp(logP):1.0-exp(logP);
    	}
    }

    inline double pgenhypergeometric(
    	int x,
    	double a,
    	double n,
    	double N,
    	hyperType variety
    )
    {
    	double logP=0;
    	double b=0;
    	double temp=0;

    	switch (variety) {
    		case IAii:
    			temp=a;
    			a=n;
    			n=temp;
    			variety =IAi;
    		case IAi:
    			if (x == (int)n) {
    				return 1.0;
    			}
    			b=N-a;
    			logP=loggamma(b+1.0)+loggamma(N-n+1.0)-loggamma(b-n+1.0)-loggamma(N+1.0);
    			break;
    		case IIIA:
    			temp=a;
    			a=n;
    			n=temp;
    			variety =IIA;
    		case IIA:
    			if (x == (int)n) {
    				return 1.0;
    			}
    			b=N-a;
    			logP=loggamma(-b+n)+loggamma(-N)-loggamma(-b)-loggamma(-N+n);
    			break;
    		case IIIB:
    			temp=a;
    			a=n;
    			n=temp;
    			variety=IIB;
    		case IIB:
    			b=N-a;
    			// Can't use this because n is not an integer
    			//logP=loggamma(-b+n)+loggamma(-N)-loggamma(-b)-loggamma(-N+n);
    			logP=log(1.0/GaussianHypergometricFcn(-n,-a,b-n+1.0,1.0));
    			break;
    		case IB:
    			b=N-a;
    			logP=loggamma(b+1.0)+loggamma(N-n+1.0)-loggamma(b-n+1.0)-loggamma(N+1.0);
    			break;
    		case IV:
    			b=N-a;
    			logP=loggamma(b+1.0)+loggamma(N-n+1.0)-loggamma(b-n+1.0)-loggamma(N+1.0);
    			break;
    		default:
    			break;
    	}

    	double sum=1.0;
    	double Tr=1.0;
    	double bn=b-n;

    	for (int i=0;i<x;i++) {
    		double r=(double)i;
    		double rp=(double)(i+1);
    		Tr*=((r-a)*(r-n))/(rp*(bn+rp));
    		sum+=Tr;
    	}

    	if (! R_FINITE(sum)){
    		return NA_REAL;
    	}

    	logP += (double) log(sum);

    	if (logP < -MAXEXP) {
    		return 0.0;
    	}
    	else if (logP > 0) {
    		return 1.0;
    	}
    	else {
    		return exp(logP);
    	}


    }



    	// Finds the type of hypergeometric
    inline hyperType typeHyper(
    	double a, 				// Sample size
    	double m,      	// Total number of marked items
    	double N       		// Total number of items
    )
    {

    	hyperType variety;



    	if (0.0<a && 0.0<N && 0.0<m &&  isint(a) && isint(N) && isint(m)) {
    		variety=classic;
    	}

    	else
    	if (0.0<a && 0.0<N && 0.0<m && isint(m) && m-1.0<a && a<N-(m-1.0)) {
    		variety=IAi;
    	}
    	else
    	if (0.0<a && 0.0<N && 0.0<m && isint(a) && a-1.0<m && m<N-(a-1.0)) {
    		variety=IAii;
    	}
    	else
    	if (0.0<a && 0.0<N && 0.0<m &&  ! isint(a) && ! isint(m) && a+m-1.0<N &&
    		 		floor(a) == floor(m)) {
    		variety=IB;		// Specified 1.0<N to avoid problems with small parameters
    	}
    	else
    	if (a<0.0 && N<m+a-1.0 && 0.0<m && isint(m))	{ //Kemp&Kemp use b<0 && b!=-1, Ben Bolker mod
    		variety=IIA;
    	}
    	else
    	if (a<0.0 && -1.0<N && N<m+a-1.0 && 0.0<m && ! isint(m) &&
    				floor(m) == floor(m+a-1.0-N)) {
    		variety=IIB;
    	}
    	else
    	if (0.0<a && N<m-1.0 && m<0.0 && isint(a)) {
    		variety=IIIA;
    	}
    	else
    	if (0.0<a && -1.0<N && N<a+m-1.0 && m<0.0 && ! isint(a) &&
    				floor(a) == floor(a+m-1.0-N)) {
    		variety=IIIB;
    	}
    	else
    	if (a<0.0 && -1.0<N && m<0.0) {
    		variety=IV;
    	}
    	else {
    		variety=noType;
    	}

    	return variety;
    }

    /* End code from SuppDists */
  } // namespace <anonymous>

  template<typename T>
  concept DoubleOrInt = std::is_same_v<T, int> || std::is_same_v<T, double>;

  template<DoubleOrInt T_type>
  double pbnbinom(T_type q, double bnb_k, double bnb_alpha, double bnb_beta, bool lower_tail, bool lg){

    double p = NA_REAL;

    if constexpr (std::is_same_v<T_type, double>)
    {
      /*
      TODO: continuous approximation
      If the pre-treatment total count is high enough, then the beta conjugate will be close to gamma
      (as bnb_beta will be large). Also, ignoring Poisson variation means that we have a gamma distribution,
      with a gamma conjugate prior (with the same hyperparameters as the beta, i.e. alpha is the product of
      N and known shape parameter, and beta is the sum of observed data).
      If the post-treatment total count is also high enough, then we can ignore Poisson variation.
      We then have a compound gamma distribution, which is a special case of the generalised beta prime
      (https://en.wikipedia.org/wiki/Beta_prime_distribution#Generalization).
      The gbeta distribution (https://cran.r-project.org/web/packages/gbeta/ - https://github.com/stla/gbeta)
      implements this. Or we could just ignore this and use a gamma approximation to a compound gamma-Poisson
      i.e. calculate mean/variance of the BNB (https://en.wikipedia.org/wiki/Beta_negative_binomial_distribution)
      and derive parameters for a gamma accordingly...?
      Limitation:  we need alpha > 2 (product of k and N) for the BNB to have finite variance - this is not
      guaranteed purely from a high sum of pre-treatment data, as it could be a small sample but high mean
      */

    }else if constexpr (std::is_same_v<T_type, int>)
    {
    	// We redefine upper tail as inclusive:
    	if(q==0L){
    		return(1.0);
    	}else{
    	  q--;
    	}

    	// Convert from beta-NB to generalised hypergeometric:
    	double ghg_a = -bnb_beta;
    	double ghg_k = -bnb_k;
    	double ghg_N = bnb_alpha - 1.0;

    	hyperType variety = typeHyper(ghg_a, ghg_k, ghg_N);
    	if (! checkHyperArgument(q, ghg_a, ghg_k, ghg_N, variety)){
    	  p = NA_REAL;
    	}else if (variety==classic){
    	  p = phypergeometric(q, static_cast<int>(ghg_a), static_cast<int>(ghg_k), static_cast<int>(ghg_N));
    	}else{
    	  p = pgenhypergeometric(q, ghg_a, ghg_k, ghg_N, variety);
    	}
    	if(!lower_tail){
    	  p = 1.0 - p;
    	}
      if(lg){
        p = std::log(p);
      }

    }else
    {
      Rcpp::stop("This should not be possible given the concept/constraint used...!");
    }

    return(p);
}


} // namespace bayescount

#endif // BAYESCOUNT__PBNBINOM_H
