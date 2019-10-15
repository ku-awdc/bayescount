# News and version history of the bayescount package

After a relatively long period of inactivity, this package is currently undergoing a complete redesign (and overhaul of the code base) towards a CRAN release of version 1.  Contributions and feedback are welcome, so if you find a bug, or if you have any comments/criticisms/questions about this package, then please feel free to email the package developer.

# Open Issues

## Bugs

* The approximation-of-last-resort in fecrt.cc::bnbpval needs to be fixed

## Features to be implemented

* Migration of older functions from version 0.9.x is not yet complete, so some older code may not yet work with this version of the package.  In particular, updates to the MCMC-based parameter estimation are planned.
* Vignettes and help files are currently incomplete
* Automated package tests are currently lacking
* Additional tweaks to the shiny interfaces are planned
* Code for the generalised hypergeometric distribution will hopefully be removed from this package once the direct interface is available within the SuppDists package


# Version history

## Version 1.0.0-10 - October 2019

### New features

* Migration of the package to its new home at https://github.com/ku-awdc/bayescount
* Restructing of the package to conform to a tidy structure
* Incorporation of underlying C/C++ code for the statistical methods outlined in "A hypothesis testing framework for the ratio of means of two negative binomial distributions: classifying the efficacy of anthelmintic treatment against intestinal parasites"
* Addition of shiny interfaces referenced from http://fecrt.com/framework

### Bug fixes

None


## Version 0.9.99 - April 2015

### New features

* A more robust fecrt function with more options

* Functions to be removed from version 1 are now deprecated

### Bug fixes

* Issue fixed preventing the citation file being parsed correctly by CRAN

* New contact details and bug fixes in manuals



## Version 0.9.9

### New features

* Help file concerning the likelihoods produced by bayescount and bayescount.single updated for clarity.

* power analysis

* fecrt function now takes confidence parameter - used for all CI intervals and boot/mcmc classifications (from probability)

* fecrt function and all gamma models in bayescount now use a T(10^-200,) on the gamma distributions to prevent crashes caused by a failure to calculate log densities associated with very small gamma values


### Bug fixes

* Bug fixed which could cause a crash with bayescount when using the scale mean option with data consisting of all zeros

* Bug fixed which would prevent bayescount from running with data that could have been multiplied by a constant

* Citation information updated for compatibility with R version 2.10

* Numerous other minor bug fixes


## Version 0.9.0

### New features

* Functions for running JAGS developed into 'runjags' package, which is now a requirement for bayescount.

* Running of models to convergence and calculation of necessary run time is now handled using the autorun.jags function.

* NAMESPACE file used now

* Model specification for (ZI)GP and LP models slightly altered.  The new model specifications produce superior results.

* bayescount.single now calculates 95% credible intervals using highest posterior density estimates rather than quantiles

* Likelihood calculations now also return the maximum likelihood estimate

* Repeat counts can now be analysed (although the likelihood cannot be calculated with irregular numbers of counts or with any repeat counts for the Weibull models)

* Coefficient of variation can now be specified as a parameter to lnormal.params and normal.params, and functions now return lists rather than a matrix.

* FECRT function now implemented.  See ?fecrt

* Several other new features


### Bug fixes

* Superficial changes to the feedback provided by the likelihood function.

* Superficial changes to wording of feedback provided by bayescount and bayescount.single.

* Several other bug fixes.


## Older versions

* Function to call JAGS for any given model/data/inits with any number of chains created November 2007.

* Modified to use JAGS 0.9.99 and lognormal / Weibul / indpendant Poisson models November 2007.  Also modified to extend rather than re-run for extra convergence iterations.

* Function to analyse faecal egg count data with either a (zero-inflated) Gamma Poisson, Lognormal Poisson, Weibul Poisson, independant Poisson or single Poisson model then check converegence with 2 chains created October 2007.

* Source script to handle/call function to analyse faecal egg count data created 3rd July 2007.


