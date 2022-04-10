## Testing estimate class

pre <- rnbinom(10, 1, mu=20)
post <- rnbinom(10, 1, mu=5)

est <- bayescount:::Rcpp_estimator_pair_unfix$new(100, 1.0, 1.0)
est$push_data(pre, post)
est$estimate()[1:4]

1 - mean(post)/mean(pre)
var(pre)
var(post)
cov(pre,post)

