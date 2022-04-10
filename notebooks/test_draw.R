## Testing multiariate draws

n <- 1e5
mus <- c(10,10)
cvs <- c(10,1)
lcor <- 0

tt1 <- bayescount:::Rcpp_draw_lambda(n, mus, cvs, lcor, "lnormpois")
tt2 <- bayescount:::Rcpp_draw_count(n, mus, cvs, lcor, "lnormpois")
colMeans(tt1); colMeans(tt2)

tt1 <- bayescount:::Rcpp_draw_lambda(n, mus, cvs, lcor, "mvlnormpois")
tt2 <- bayescount:::Rcpp_draw_count(n, mus, cvs, lcor, "mvlnormpois")
colMeans(tt1); colMeans(tt2)

tt1 <- bayescount:::Rcpp_draw_lambda(n, mus, cvs, lcor, "nbinom")
tt2 <- bayescount:::Rcpp_draw_count(n, mus, cvs, lcor, "nbinom")
colMeans(tt1); colMeans(tt2)

tt1 <- bayescount:::Rcpp_draw_lambda(n, mus, cvs, lcor, "mvnbinom")
tt2 <- bayescount:::Rcpp_draw_count(n, mus, cvs, lcor, "mvnbinom")
colMeans(tt1); colMeans(tt2)
apply(tt1,2,sd)/colMeans(tt1); apply(tt2,2,sd)/colMeans(tt2)
qqplot(rgamma(1e5, 1/cvs[1]^2, mus[1]/cvs[1]^2), tt1[,1], pch=20); abline(0,1)
qqplot(log(rgamma(1e5, 1/cvs[1]^2, mus[1]/cvs[1]^2)), log(tt1[,1]), pch=20); abline(0,1)


tt1 <- bayescount:::Rcpp_draw_lambda(n, mus, c(0,0), lcor, "poisson")
tt2 <- bayescount:::Rcpp_draw_count(n, mus, c(0,0), lcor, "poisson")
colMeans(tt1); colMeans(tt2)



# Means are over-high!
lc <- seq(-1,1,by=0.05)
ec <- sapply(lc, function(x){
	tt <- bayescount:::Rcpp_draw_lambda(1e5, c(10,10), c(10,1), x, "mvlnormpois")
	cor(tt[,1], tt[,2])
})

plot(lc, ec, type='l'); abline(0,1)


tt <- bayescount:::Rcpp_draw_lambda(1e5, c(10,10), c(1,1), 0.5, "mvlnormpois")
colMeans(tt)
apply(tt,2,sd)/colMeans(tt)
cor(tt[,1], tt[,2])
cor(log(tt[,1]), log(tt[,2]))

mean <- 10
tau <- 1
(lsd <- sqrt(log(tau^2 +1)))^2
log(mean)
(lmu <- log(mean) - ((lsd^2) / 2))

mean(rlnorm(1e5, lmu, lsd))


tt <- bayescount:::Rcpp_draw_count(1e5, c(10,10), c(1,1), 0.5, "mvlnormpois")
head(tt)
colMeans(tt)
apply(tt,2,sd)/colMeans(tt)
cor(tt[,1], tt[,2])
