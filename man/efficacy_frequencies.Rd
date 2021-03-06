% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/efficacy_frequencies.R
\name{efficacy_frequencies}
\alias{efficacy_frequencies}
\title{Obtain a Monte Carlo estimation of the expected frequency of typologies associated with the given parameter values}
\usage{
efficacy_frequencies(
  r = NA,
  paired = TRUE,
  T_I = 0.99,
  T_A = 0.95,
  N = c(20, 20),
  R = c(1, 1),
  S = c(1, 1),
  mean = 20,
  k = c(1, 0.7),
  cor = 0.5,
  iterations = 10^3,
  alpha = 0.025,
  bnb_priors = c(0, 0),
  use_delta = NA,
  beta_iters = 10^4,
  use_ml = TRUE,
  binomial_priors = c(1, 1),
  binomial_cl_adj = 0.2
)
}
\arguments{
\item{r}{a vector representing a set of true population reductions to simulate.  If NA then values of r over an appropriate range of interest (based on T_I and T_A) are automatically generated.}

\item{paired}{logical flag for a paired or unpaired study design assumption on which to simulate (and analyse) the data}

\item{T_I}{the threshold for inferioirty (target efficacy of the intervention)}

\item{T_A}{the threshold for non-inferioirty (target efficacy of the intervention minus a non-inferioirty margin delta)}

\item{N}{a length-2 integer giving the number of control and treatment individuals for an unpaired study, or a length-1 integer giving the number of individuals for a paired study}

\item{R}{a length-2 integer giving the number of replicates used for control/pre-treatment and treatment/post-treatment individuals}

\item{S}{a length-2 numeric variable giving the counting sensitivity used for control/pre-treatment and treatment/post-treatment individuals}

\item{mean}{the control/pre-treatment mean count}

\item{k}{a length-2 strictly positive continuous number giving the over-dispersion of the control/pre-treatment and treatment/post-treatment data}

\item{cor}{a length-1 number between 0 and 1 giving the correlation between pre- and post-treatment data (ignored for an unpaired study)}

\item{iterations}{the number of iterations to use for the Monte Carlo approximation}

\item{alpha}{the significance level to use for the classification, where 0.025 corresponds to 95\\% CI and 0.05 corresponds to 90\\% CI}

\item{bnb_priors}{priors to use for the BNB method}

\item{use_delta}{logical flag to use the delta method approximation for the BNB method (NA means to use it unless it fails)}

\item{beta_iters}{number of iterations to use for the Monte Carlo approximation of the beta distribution transformation for the BNB method (when use_delta==FALSE)}

\item{use_ml}{a logical flag controlling if maximum likelihood is used (rather than less accurate but faster methods) to estimate the over-dispersion parameter k}

\item{binomial_priors}{the priors to use for the Binomial method}

\item{binomial_cl_adj}{the adjustment in confidence level to use for the Binomial method relative to the others (the default corresponds to e.g. 99\\% rather than 95\\% CI) - where 0 > adj >= 1}
}
\value{
A data frame with columns reflecting the method, population reduction r, typology, classification, frequency, and proportion within method & r.  The parameter values used for the simulation are given as an attribute list.
}
\description{
this calls underlying C++ code to repeatedly simulate data from a (paired) set of negative binomial distributions with given parameters, analyse the resultant data with each of four/five statistical methods, and return the expected frequencies of each typology encountered using each method
}
\examples{
# Monte Carlo estimates of typology frequencies for all five methods
# at the critical values corresponding to T_I and T_A:
( critical_values <- efficacy_frequencies(r=c(0.95, 0.99), paired = TRUE, T_I = 0.99, T_A = 0.95) )
# The same for an unpaired analysis (excluding the Binimial method):
( critical_values <- efficacy_frequencies(r=c(0.95, 0.99), paired = FALSE, T_I = 0.99, T_A = 0.95) )
}
\seealso{
\code{\link{bayescount}}
}
