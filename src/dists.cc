/*
This file was adapted from SuppDists version 1.1-9
- Some code not necessary for implementation of the generalised hypergeometric has been removed
- Signature for first argument in hypergeometric functions has been changed from int to long long
Original file is copyright of Bob Wheeler, licensed GPL>=2
https://CRAN.R-project.org/package=SuppDists
Last updated by MJD February 2017 as part of the bayescount package
*/


// needs to come first, BDR 2013-10-30, and for gcc 6 we need NO_C_HEADERS
#include <new>
#include <cmath>
#include <cstdlib>



// for Solaris
using namespace std;
#if defined(__SUNPRO_CC) && (__cplusplus < 201103L)
// This is not C++98 and Solaris CC's headers do not include it
extern "C" int snprintf(char *str, size_t size, const char *format, ...);
#endif

//#include <float.h>

// No longer needed in R 3.3.0 and later:
// #define NO_C_HEADERS

#include <R.h>
#include <Rmath.h>

// snprintf is in fact not C++, but R.h used to include this
#include <stdio.h> // for snprintf
#include <cstring> // for memset

#include "dists.h"
//#include "datatabs.h"

// SuppDists by Robert E. Wheeler, March 2001

static const double LOG10=2.3025850929940456840179915;
static const double MAXEXP=LOG10*DBL_MAX_10_EXP;	// Maximum argument for exp()
static const double LOGSQRT2PI=0.9189385332046727417803296;
//static const double NA=-1e-12;


/*
	random normal deviates
*/	

void rgauss(
	double* normArray,
	int n,
	double mean,
	double sd
)
{ 
	int i;

	GetRNGstate();
	for (i=0;i<n;i++)
		normArray[i]=rnorm(mean,sd);
	PutRNGstate();


}

// Random chi squre
void	rdchisq(
	double *tArray,
	int n,
	int df
)
{
	int i;
	GetRNGstate();
	for (i=0;i<n;i++)
		tArray[i]=rchisq((double)df);
	PutRNGstate();
}

/*
	Derivitve of chisquared density
*/
double fpchisq(
	double x,
	int df
)
{
	double nu=(double)df/2.0;
	return dchisq(x,df,false)*((nu-1.0)/x-0.5);
}

/******************************************************************************
  	Inverse Gaussian -- Wald
*/



// UTILITIES *******************************************************************


/* 
	Natural logarithm of the gamma function.
	     CACM 291 due to M.C. Pike and I.D. Hill 
	     Accurate to at least 10 places. 
*/
double loggamma(double x)	 {
	const double  T1=1.0/12.0;
	const double  T3=1.0/360.0;
	const double  T5=1.0/1260.0;
	const double  T7=1.0/1680.0;
	const double  T9=1.0/1188.0;

	double f;

	if (x equals 1.0 || x equals 2.0) {
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


/* normdrv -- derivatives of normal */

/* normal density */
double phi0(double x)
{
   double constant=.398942280401433;

   return (constant*exp(-0.5*x*x));
}

/* third derivative z=phi0(x)*/
double phi3(double x,double z)
{
   return (z*x*(3.0-x*x));
}

/* fifth derivative z=phi0(x)*/
double phi5(double x,double z)
{
   double s;

   s=x*x;
   return (-z*x*((s-10.0)*s+15.0));
}

/* seventh derivative z=phi0(x)*/
double phi7(double x,double z)
{
   double s;

   s=x*x;
   return(z*x*(((-s+21.0)*s-105.0)*s+105.0));
}

/* Permute ***************************************************************
|  Randomly pemutes the n integers in a[] using the Fike
|  algorithm.  See Fike, "A permutation generation method"  The Computer
|  Journal, 18-1, Feb 75, 21-22.
*/

void Permute(
	int* a,
	int n
)
{
   int i;
   int j;
   int temp;

	GetRNGstate();
   for (i=1;i<n;i++) {
		j=(int)((double)(1+i)*unif_rand());  
      temp=a[j];
      a[j]=a[i];
      a[i]=temp;
   }
   PutRNGstate();
}







/**************************************************************************************** 
	Hypergeometric distribution 
	Given a total of N items, n of which are marked, select a sample of size S, and
	 find the probability that 
	 
	fhypergeometric:  that there are x marked items in the sample -- frequency
	phypergeometric:  that there are x or less in the sample -- cdf.
*/



 	// Normal approximation to the hypergeometric distribution function due to Peizer
	// See Ling, R.F. and Pratt, J.W. (1984) The accuracy of Peizer approximations
	//  to the hypergeometric distribution, with comparisons to some other 
	//  approximations. JASA 79-385. 49-60.
double PeizerHypergeometric(
	long long x,	  // Number of marked items in sample
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
	
	return pnorm(z,0,1,true,false);

}

char *hyperNames[]= {
	(char *)"classic",
	(char *)"IAi",
	(char *)"IAii",
	(char *)"IB",
	(char *)"IIA",
	(char *)"IIB",
	(char *)"IIIA",
	(char *)"IIIB",
	(char *)"IV",
	(char *)"no type"
};

	// Returns true if the double is an int
bool isint(
	double x
)
{
	return x equals floor(x);
}      



/* These are additional implementations using long long rather than int for bayescount */

bool checkHyperArgument(
	long long k,
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

double phypergeometric(
	long long x,		// Number of marked items in sample
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

	if (x equals n) {
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
	for (long long k=xmin;k<x;k++) { 
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

double qhypergeometric(
	long long x,		// Number of marked items in sample
	int a, 				// Sample size
	int n,      	// Total number of marked items
	int N       		// Total number of items
)
{
	return 1.0-phypergeometric(x,a,n,N);
}

 double pgenhypergeometric(
	long long x,	 
	double a,	 
	double n,	 
	double N,
	hyperType variety	
)
{
	double logP=0;
	double b=0;
	double temp=0;
	double P=0;

	switch (variety) {
		case IAii:
			temp=a;
			a=n;
			n=temp;
			variety =IAi;
		case IAi:
			if (x equals (int)n) {
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
			if (x equals (int)n) {
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
	
	for (long long i=0;i<x;i++) {
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

 double qgenhypergeometric(
	long long x,	 
	double a,	 
	double n,	 
	double N,
	hyperType variety
)
{
	return 1.0-pgenhypergeometric(x,a,n,N,variety);
}

/* End additional */



	// Finds the type of hypergeometric

hyperType typeHyper(
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
		 		floor(a) equals floor(m)) {
		variety=IB;		// Specified 1.0<N to avoid problems with small parameters
	}
	else
	if (a<0.0 && N<m+a-1.0 && 0.0<m && isint(m))	{ //Kemp&Kemp use b<0 && b!=-1, Ben Bolker mod
		variety=IIA;
	}
	else
	if (a<0.0 && -1.0<N && N<m+a-1.0 && 0.0<m && ! isint(m) && 
				floor(m) equals floor(m+a-1.0-N)) {
		variety=IIB;
	}
	else
	if (0.0<a && N<m-1.0 && m<0.0 && isint(a)) {
		variety=IIIA;
	}
	else
	if (0.0<a && -1.0<N && N<a+m-1.0 && m<0.0 && ! isint(a) && 
				floor(a) equals floor(a+m-1.0-N)) {
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

bool checkHyperArgument(
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

	
	// Reports the type and range to the user
DISTS_API void tghyperR(
	double *ap, 				// Sample size
	double *mp,      	// Total number of marked items
	double *Np,       		// Total number of items
	char **aString
)
{
	double a=*ap;
	double m=*mp;
	double N=*Np;

	hyperType variety=typeHyper(a,m,N);


	switch (variety) {
		case classic:
			snprintf(*aString,127,"type = %s -- %d <= x <= %d",hyperNames[(int)classic],maxm(0,(int)(a+m-N)),minm((int)a,(int)m));
			break;

		case IAi:
			snprintf(*aString,127,"type = %s -- 0 <= x <= %d",hyperNames[(int)IAi],(int)m);
			break;
		case IAii:
			snprintf(*aString,127,"type = %s -- 0 <= x <= %d",hyperNames[(int)IAii],(int)a);
			break;
		case IB:		// Specified 1.0<N to avoid problems with small parameters
			snprintf(*aString,127,"type = %s -- x = 0,1,2,...",hyperNames[(int)IB]);
			break;
		case IIA:
			snprintf(*aString,127,"type = %s -- 0 <= x <= %d",hyperNames[(int)IIA],(int)m);
			break;
		case IIB:
			snprintf(*aString,127,"type = %s -- x = 0,1,2,...",hyperNames[(int)IIB]);
			break;
		case IIIA:
			snprintf(*aString,127,"type = %s -- 0 <= x <= %d",hyperNames[(int)IIIA],(int)a);
			break;
		case IIIB:
			snprintf(*aString,127,"type = %s -- x = 0,1,2,...",hyperNames[(int)IIIB]);
			break;
		case IV:
			snprintf(*aString,127,"type = %s -- x = 0,1,2,...",hyperNames[(int)IV]);
			break;
		case noType:
			snprintf(*aString,127,"type = %s",hyperNames[(int)noType]);
	}


}

DISTS_API void dghyperR(
	int *kp,
	double *ap,
	double *np,
	double *Np,
	int *Mp,
	double *valuep
)
{
	hyperType variety;
	int M=*Mp;
	int i;

	for (i=0;i<M;i++) {
		variety=typeHyper(ap[i],np[i],Np[i]);
		if (variety==classic)
			valuep[i]=fhypergeometric(kp[i],(int)ap[i],(int)np[i],(int)Np[i]);
		else if (variety!=noType)
			valuep[i]=fgenhypergeometric(kp[i],ap[i],np[i],Np[i],variety);
		else valuep[i]=NA_REAL;
	}

		
}


double fhypergeometric(
	int x,		// Number of marked items in sample
	int S, 				// Sample size
	int n,      	// Total number of marked items
	int N       		// Total number of items
)
{
	double logP=loggamma((double)(n+1))+loggamma((double)(N-n+1))+loggamma((double)(S+1))+
		loggamma((double)(N-S+1));
	logP-=loggamma((double)(x+1))+loggamma((double)(n-x+1))+loggamma((double)(S-x+1))+
		loggamma((double)(N-S-n+x+1))+loggamma((double)(N+1));
	if (logP<-MAXEXP) {
		return 0.0;
	}
	else {
		return exp(logP);
	}

}


/*
	Essentially, one has a table:

		x  y | n
		?  ? | N-n
		----------
		a N-a| N
	
	(1) x ranges from a=max(0,a+n-N)) to b=min(n,a). If n<a, then the total range
			is N-A, for A=a+n-N or n for A=0.
	(2) One can interchange n and a if a is smaller, to minimize this range.
	(3) If a<N-a, then A=0, otherwise, one can interchange x and y and a and N-a to
		get A=0.
		Proof: if a<N-a, then A=a+n-N<N-a+n-N=n-a, and n<a hence A=0 by interchange.
*/


DISTS_API void pghyperR(
	int *kp,
	double *ap,
	double *np,
	double *Np,
	int *Mp,
	double *valuep
)
{
	hyperType variety;
	int M=*Mp;
	int i;

	for (i=0;i<M;i++) {
		variety=typeHyper(ap[i],np[i],Np[i]);
		if (! checkHyperArgument(kp[i],ap[i],np[i],Np[i],variety))
			valuep[i]=NA_REAL;
		else if (variety==classic)
			valuep[i]=phypergeometric(kp[i],(int)ap[i],(int)np[i],(int)Np[i]);
		else
			valuep[i]=pgenhypergeometric(kp[i],ap[i],np[i],Np[i],variety);
	}

		
}

double phypergeometric(
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

	if (x equals n) {
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


DISTS_API void ughyperR(
	int *kp,
	double *ap,
	double *np,
	double *Np,
	int *Mp,
	double *valuep
)
{
	hyperType variety;
	int M=*Mp;
	int i;

	for (i=0;i<M;i++) {
		variety=typeHyper(ap[i],np[i],Np[i]);
		if (! checkHyperArgument(kp[i],ap[i],np[i],Np[i],variety))
			valuep[i]=NA_REAL;
		else if (variety==classic)
			valuep[i]=qhypergeometric(kp[i],(int)ap[i],(int)np[i],(int)Np[i]);
		else
			valuep[i]=qgenhypergeometric(kp[i],ap[i],np[i],Np[i],variety);
	}

		
}


double qhypergeometric(
	int x,		// Number of marked items in sample
	int a, 				// Sample size
	int n,      	// Total number of marked items
	int N       		// Total number of items
)
{
	return 1.0-phypergeometric(x,a,n,N);
}

DISTS_API void qghyperR(
	double *pp,
	double *ap,
	double *np,
	double *Np,
	int *Mp,
	double *valuep
)
{
	hyperType variety;
	int M=*Mp;
	int i;

	for (i=0;i<M;i++) {
		variety=typeHyper(ap[i],np[i],Np[i]);
		if (variety==classic)
			valuep[i]=xhypergeometric(pp[i],(int)ap[i],(int)np[i],(int)Np[i]);
		else if (variety!=noType)
			valuep[i]=xgenhypergeometric(pp[i],ap[i],np[i],Np[i],variety);
		else valuep[i]=NA_REAL;
	}

		
}

	// returns smallest x such that  p<=Pr(X<=x|a,n,N) 
int  xhypergeometric(
	double p,		// cumulative probability
	int a, 				// Sample size
	int n,      	// Total number of marked items
	int N       		// Total number of items
)
{
	double T=qchisq(1.0-p,1,true,false);
	double z=(T*(p*(1.0-p)*(double)(a*(N-a))))/(double)(N-1);
	int x=(int)floor(0.5+p*(double)a+z*z);

	int minX=maxm(0,n+a-N);	
	int maxX=minm(n,a);
	x=maxm(x,minX);
	x=minm(x,maxX);

	if (0>p || p>1.)
		error("\nProbability must be in the 0 to 1 range");

   	bool larger=(p<=phypergeometric(x,a,n,N));
   	while (larger) {
      	if (x equals minX) {
      		return x;
		}
      	larger=(p<=phypergeometric(--x,a,n,N));
      	if (! larger) {
      		return ++x;
		}
   	}
    while (! larger){
      	larger=(p<=phypergeometric(++x,a,n,N));
      	if (larger) {
      		return x;
		}
	}

	return 0;
}


DISTS_API void rghyperR(
	double *ap,
	double *np,
	double *Np,
	int *Mp,
	int *Kp,
	double *valuep
)
{
	hyperType variety;

	int M=*Mp;
	int K=*Kp;
	int D;
	int j;
	int k;
	int loc;
	int cloc;
	double *tArray;

	if (K==1) {
		variety=typeHyper(*ap,*np,*Np);
		if (variety==classic)
			rhypergeometric(valuep,M,(int)*ap,(int)*np,(int)*Np);
		else if (variety!=noType)
			rgenhypergeometric(valuep,M,*ap,*np,*Np,variety);
		else 
			error("\nParameters are for no recognized type");
	}
	else { // Allow for random values for each element of nu and lambda
		D=(M/K)+((M%K)?1:0);
		tArray=(double *)S_alloc((long)D,sizeof(double));
		loc=0;
		for (j=0;j<K;j++) {
			variety=typeHyper(ap[j],np[j],Np[j]);
			if (variety==classic)
				rhypergeometric(tArray,D,(int)ap[j],(int)np[j],(int)Np[j]);
			else if (variety!=noType)
				rgenhypergeometric(tArray,D,ap[j],np[j],Np[j],variety);
			else 
				error("\nParameters are for no recognized type");
			for (k=0;k<D;k++) {
				cloc=loc+k*K;
				if (cloc<M)
					valuep[cloc]=tArray[k];
				else break;
			}
			loc++;
		}
	}

}



/*
	Random samples from hypergeometric
*/
void rhypergeometric(
	double* randArray,
	int n,	  // number of samples
	int a, 				// Sample size
	int m,      	// Total number of marked items
	int N       		// Total number of items
)
{
	GetRNGstate();

	for (int i=0;i<n;i++){	
		randArray[i]=(double)xhypergeometric(unif_rand(),a,m,N);
	}
	PutRNGstate();
}

DISTS_API void sghyperR(
	double *ap,
	double *mp,
	double *Np,
	int *Mp,
	double *meanp,
	double *medianp,
	double *modep,
	double *variancep,
	double *thirdp,
	double *fourthp
)
{
	int M=*Mp;
	int i;
	hyperType variety;

	for (i=0;i<M;i++) {
		variety=typeHyper(ap[i],mp[i],Np[i]);
		sghyper(ap[i],mp[i],Np[i],meanp+i,medianp+i,modep+i,variancep+i,thirdp+i,fourthp+i,variety);
	}
}

void sghyper(
	double a,
	double m,
	double N,
	double *mean,
	double *median,
	double *mode,
	double *variance,
	double *third,
	double *fourth,
	hyperType variety

)
{
	bool paramSet=false;
	double n=0;
	double A=0;
	double B=0;
	double T=0;

	double m1=0;
	double m2=0;
	double m3=0;
	double m4=0;


	switch (variety) {
		case IIIB:
			A=minm(m,a);
			n=maxm(m,a);
			B=N-A;
			paramSet=true;
		case IIB:
			if (! paramSet) {
				A=minm(m,a);
				n=maxm(m,a);
				B=N-A;
			}
 			T=A+B;
			*mode=(int)n+1;
			*median=(double)xgenhypergeometric(0.5,a,m,N,variety);
			*mean=NA_REAL;
			*variance=NA_REAL;
			*third=NA_REAL;
			*fourth=NA_REAL;
			break;
		case classic:
			n=minm(m,a);
			A=maxm(m,a);
			B=N-maxm(m,a);
			*median=(double)xhypergeometric(0.5,(int)a,(int)m,(int)N);
			paramSet=true;
		case IAi:
			if (! paramSet) {
				n=minm(m,a);
				A=maxm(m,a);
				B=N-A;
				*median=(double)xgenhypergeometric(0.5,a,m,N,variety);
				paramSet=true;
			}
		case IAii:
 			if (! paramSet) {
				n=minm(m,a);
				A=maxm(m,a);
				B=N-A;
				*median=(double)xgenhypergeometric(0.5,a,m,N,variety);
				paramSet=true;
			}
		case IIA:
 			if (! paramSet) {
				A=minm(m,a);
				n=maxm(m,a);
				B=N-A;
				*median=(double)xgenhypergeometric(0.5,a,m,N,variety);
				paramSet=true;
			}
		case IIIA:
 			if (! paramSet) {
				A=minm(m,a);
				n=maxm(m,a);
				B=N-A;
				*median=(double)xgenhypergeometric(0.5,a,m,N,variety);
				paramSet=true;
			}

 			T=A+B;

			if (n<=1.0) {
				*mean=0.0;
			}
			else {
				*mean=m1=(n*A)/T;
			}
			*mode=floor(((n+1.0)*(A+1.0))/(T+2.0));
			if (n<=2.0) {
				*variance=0.0;
			}
			else {
				m2=m1*(B*(T-n))/(T*(T-1.0));
				*variance=m2;
			}
			if (n<=3.0) {
				*third=0.0;
			}
			else {
				m3=m2*((B-A)*(T-2.0*n))/(T*(T-2.0));
				*third=m3;
			}
			if (n<=4.0) {
				*fourth=0.0;
			}
			else {
				m4=(m2/((T-2.0)*(T-3.0)))*(T*(T+1.0-6.0*n)+3.0*A*B*(n-2.0)+
						6.0*n*n+(3.0*A*B*n*(6.0-n))/(T)-(18.0*A*B*n*n)/(T*T));
				*fourth=m4;
			}
			break;
		case IB:
			if (! paramSet) {
				paramSet=true;
			}
		case IV:
			if (! paramSet) {
				paramSet=true;
			}

			*median=(double)xgenhypergeometric(0.5,a,m,N,variety);

			A=minm(m,a);
			n=maxm(m,a);
			B=N-A;
  			T=A+B;

			if (T<=0.0) {
				*mean=NA_REAL;
			}
			else {
				*mean=m1=(n*A)/T;
			}
			*mode=floor(((n+1.0)*(A+1.0))/(T+2.0));
			if (T<=1.0) {
				*variance=NA_REAL;
			}
			else {
				m2=m1*(B*(T-n))/(T*(T-1.0));
				*variance=m2;
			}
			if (T<=3.0) {
				*third=NA_REAL;
			}
			else {
				m3=m2*((B-A)*(T-2.0*n))/(T*(T-2.0));
				*third=m3;
			}
			if (T<=4.0) {
				*fourth=NA_REAL;
			}
			else {
				m4=(m2/((T-2.0)*(T-3.0)))*(T*(T+1.0-6.0*n)+3.0*A*B*(n-2.0)+
						6.0*n*n+(3.0*A*B*n*(6.0-n))/(T)-(18.0*A*B*n*n)/(T*T));
				*fourth=m4;
			}
			break;

		default:
			break;
	}

}


// Generalized hypergeometric

 double pgenhypergeometric(
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
	double P=0;

	switch (variety) {
		case IAii:
			temp=a;
			a=n;
			n=temp;
			variety =IAi;
		case IAi:
			if (x equals (int)n) {
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
			if (x equals (int)n) {
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
			P=1.0/GaussianHypergometricFcn(-n,-a,b-n+1.0,1.0);
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

	Rprintf("%f\n", logP);
	
 	double sum=1.0;
	double Tr=1.0;
	double bn=b-n;

	for (int i=0;i<x;i++) {
		double r=(double)i;
		double rp=(double)(i+1);
		Tr*=((r-a)*(r-n))/(rp*(bn+rp));
		sum+=Tr;
	}

	Rprintf("%f\n", sum);
	
	if (variety equals IIB) {
		P*=sum;
		return minm(P,1.0);	 // Occasional numerical error
	}
	else {
		logP+=log(sum);
		
		if (logP<-MAXEXP) {
			return 0.0;
		}
		else {
			return exp(logP);
		}
	}


}

 double qgenhypergeometric(
	int x,	 
	double a,	 
	double n,	 
	double N,
	hyperType variety
)
{
	return 1.0-pgenhypergeometric(x,a,n,N,variety);
}

 double fgenhypergeometric(
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
	double P=0;

	switch (variety) {
		case IAii:
			temp=a;
			a=n;
			n=temp;
			variety =IAi;
		case IAi:
			b=N-a;
			logP=loggamma(a+1.0)+loggamma(b+1.0)+loggamma(n+1.0)+loggamma(N-n+1.0);
			logP-=loggamma(x+1.0)+loggamma(a-x+1.0)+loggamma(n-x+1.0)+loggamma(b-n+x+1.0)+loggamma(N+1.0);
			break;
		case IIIA:
			temp=a;
			a=n;
			n=temp;
			variety =IIA;
		case IIA:
			b=N-a;
			logP=loggamma(-a+x)+loggamma(-b+n-x)+loggamma(n+1.0)+loggamma(-N);
			logP-=loggamma(x+1.0)+loggamma(-a)+loggamma(n-x+1.0)+loggamma(-b)+loggamma(-N+n);
			break;
		case IIIB:
			temp=a;
			a=n;
			n=temp;
			variety=IIB;
		case IIB:
			b=N-a;
				// (1) -1<b<0.  
				// (2) n is not itegral.
				// (3) must sum terms, cannot use loggamma  
			P=1.0/GaussianHypergometricFcn(-n,-a,b-n+1.0,1.0);
			{
				double Tr=1.0;
				double bn=b-n;

				for (int i=0;i<x;i++) {
					double r=(double)i;
					double rp=(double)(i+1);
					Tr*=((r-a)*(r-n))/(rp*(bn+rp));
				}
				P*=Tr;
			}

			break;
		case IB:
			b=N-a;
				// Assuming b>0 always
			logP=loggamma(a+1.0)+loggamma(b+1.0)+loggamma(n+1.0)+loggamma(N-n+1.0);
			logP-=loggamma(x+1.0)+loggamma(a-x+1.0)+loggamma(n-x+1.0)+loggamma(b-n+x+1.0)+loggamma(N+1.0);
			break;
		case IV:
			b=N-a;
			logP=loggamma(-a+x)+loggamma(b+1.0)+loggamma(-n+x)+loggamma(N-n+1.0);
			logP-=loggamma(x+1.0)+loggamma(-a)+loggamma(b-n+x+1.0)+loggamma(-n)+loggamma(N+1.0);
			break;
		default:
			break;
	}


	if (variety equals IIB) {
		return P;
	}
	else {
		if (logP<-MAXEXP) {
			return 0.0;
		}
		else {
			return exp(logP);
		}
	}
}


	// returns smallest x such that  p<=Pr(X<=x|a,n,N) 
int  xgenhypergeometric(
	double p,		// cumulative probability
	double a,	 
	double m,	 
	double N,
	hyperType variety
)
{
	double b=N-a;
	double n=m;

	double m1=(n*a)/N;
	double m2=(m1*(b*(a+b-n)))/(N*(N-1.0));

	if (0>p || p>1)
		error("\nProbability must be in the 0 to 1 range");

	int x=(int)(0.5+m1+sqrt(m2)*qnorm(p,0,1,true,false));
	x=maxm(0,x);


   	bool larger=(p<=pgenhypergeometric(x,a,m,N,variety));
   	while (larger) {
      	if (x equals 0) {
      		return x;
		}
      	larger=(p<=pgenhypergeometric(--x,a,m,N,variety));
      	if (! larger) {
      		return ++x;
		}
   	}
    while (! larger){
      	larger=(p<=pgenhypergeometric(++x,a,m,N,variety));
      	if (larger) {
      		return x;
		}
	}

	return 0;
}


/*
	Random samples from beta negative binomial
*/
 void rgenhypergeometric(
	double* randArray,
	int K,	  
	double a,	 
	double n,	 
	double N,
	hyperType variety
)
{
	GetRNGstate();

	for (int i=0;i<K;i++){		
		randArray[i]=(double)xgenhypergeometric(unif_rand(),a,n,N,variety);
	}

	PutRNGstate();
}




/*
	The Gaussian hypergeometric function, usually denoted as
	 F[a,b;c;x]
*/
 double GaussianHypergometricFcn(
	double a,
	double b,
	double c,
	double x
)
{
	int const MAXITERATES=100;

	if (c<0.0 && floor(c) equals c) 
		return NA_REAL;

	double sum=0.0;
	double term=1.0;
	int j=1;
	double dj;
	double djm1;
	repeat
		dj=(double)j;
		djm1=dj-1.0;
		sum+=term;
		term*=((a+djm1)*(b+djm1))/(c+djm1)*(x/dj);
		j++;
	until(sum+term equals sum || j>MAXITERATES);

	return sum;
}

