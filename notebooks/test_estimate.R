## Testing estimate class

pre <- rnbinom(10, 1, mu=20)
post <- rnbinom(10, 1, mu=5)

paired <- TRUE
k_type <- "fix"
mean_ratio <- 1.0
H0_1 <- 0.9
H0_2 <- 0.95
lci <- 0.025
uci <- 0.025
conjugate_priors <- c(0.0,0.0)
delta <- "sometimes"
beta_iters <- 1e3
approx <- "sometimes"
dobson_priors <- c(1.0,1.0)
true_effk_pre <- 1.0
true_effk_post <- 1.0

bayescount:::Rcpp_estimate_fecrt(pre, post, TRUE, k_type, mean_ratio, H0_1, H0_2, lci, uci, conjugate_priors, delta, beta_iters, approx, dobson_priors, true_effk_pre, true_effk_post)


stop()
est <- bayescount:::Rcpp_estimator_pair_unfix$new(100, 1.0, 1.0)
est$push_data(pre, post)
est$estimate()[1:4]

1 - mean(post)/mean(pre)
var(pre)
var(post)
cov(pre,post)

