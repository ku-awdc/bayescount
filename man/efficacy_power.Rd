% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/efficacy_power.R
\name{efficacy_power}
\alias{efficacy_power}
\title{Obtain a Monte Carlo estimation of the expected frequency of typologies associated with the given parameter values}
\usage{
efficacy_power(
  N = seq(5, 50, by = 5),
  method = "BNB",
  paired = TRUE,
  T_I = 0.99,
  T_A = 0.95,
  ...
)
}
\arguments{
\item{N}{the sample sizes to examine for type 1 and 2 error rates, either as an integer vector where it is presumed that N_1 = N_2, or a 2-column matrix of paired N_1 and N_2 values to examine at each combination given by row}

\item{method}{the statistical method(s) for which to return type 1 and type 2 error rates}

\item{paired}{logical flag for a paired or unpaired study design assumption on which to simulate (and analyse) the data}

\item{T_I}{the threshold for inferioirty (target efficacy of the intervention, and also the value of efficacy tested for the non-inferiority test)}

\item{T_A}{the threshold for non-inferioirty (target efficacy of the intervention minus a non-inferioirty margin delta, and also the value of efficacy tested for the inferiority test)}

\item{...}{other parameters to be passed to the underlying \code{\link{efficacy_frequencies}} and \code{\link[parallel]{mclapply}} functions}
}
\value{
A data frame with columns reflecting the method, N_1, N_2, hypothesis test, error type, and error rate.
}
\description{
this calls underlying C++ code to repeatedly simulate data from a (paired) set of negative binomial distributions with given parameters, analyse the resultant data with each of four/five statistical methods, and return the expected frequencies of each typology encountered using each method
}
\examples{
rates <- efficacy_power(N = seq(5,50,by=5), paired = TRUE, T_I = 0.99, T_A = 0.95)
if(requireNamespace("ggplot2")){
  ggplot(rates, aes(x=N_1, y=Rate, ymin=LCI, ymax=UCI, col=Test)) +
    geom_errorbar() + geom_line(lty='dashed') + geom_point() +
    geom_hline(aes(yintercept=h), data.frame(Type=c("Power","Type1Error"), h=c(0.8,0.025)), lty='dashed') +
    facet_wrap(~Type, scales='free_y') + xlab("N")
}
}
\seealso{
\code{\link{bayescount}}
}
