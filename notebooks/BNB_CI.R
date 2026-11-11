# Dependencies:
install.packages(c("dplyr", "rlang", "tidyr", "readr", "stringr", "purrr", "readxl", "R6", "shiny", "shinythemes", "rhandsontable", "ggplot2", "rmarkdown", "Rcpp", "RcppDist", "RcppArmadillo"))

# Install development version (requires compilers, and you will get many compilation warnings - but it succeeds at least on macOS):
remotes::install_github("ku-awdc/bayescount")
stopifnot(packageVersion("bayescount")>="1.2.1-3")

# For example:
N <- 20
(data1 <- rnbinom(N, mu=10, size=1))
(data2 <- rnbinom(N, mu=1*rbinom(1,1,0.5), size=1))

# To activate the secret (and largely untested) option:
Sys.setenv(BAYESCOUNT_SECRET_CI=TRUE)

# Run analysis:
efficacy_analysis(data1, data2, paired=FALSE)[c(1,4,5),1:3]

# The BNB method is based on a beta negative binomial distribution (https://arxiv.org/abs/1910.06667)
# The Gamma method is Bruno's delta method
# The WAAVP method is Lyndal-Murphy et al 2014
