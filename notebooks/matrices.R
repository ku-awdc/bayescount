library('tidyverse')

### WAAVP guidelines matrixes

Ns=N <- seq(10, 30, by=2)
margin <- c(0.005, 0.006, 0.008, 0.009, 0.01, 0.011, 0.012, 0.015, 0.02, 0.03, 0.04, 0.05)
target <- 0.99
mu <- 100.0
k_pre <- 1.0
k_post <- 1.0
k_c <- 2.0
iters <- 1000L
conjugate_priors <- c(1,1)
delta <- 1L
beta_iters <- 1000L
approx <- 1L
tail <- 0.05
useml <- TRUE

pmat <- bayescount:::RCPP_power_matrix_paired(N, margin, target, mu, k_pre, k_post, k_c, iters, conjugate_priors, delta, beta_iters, approx, tail, useml)

smat <- pmat %>%
	group_by(N, Margin, Hypothesis, Reduction) %>%
	summarise(power = mean(RejectH0), meanred = mean(ObsReduction), varred = var(ObsReduction))

smat

ggplot(smat, aes(x=N, y=power, col=factor(Margin))) +
	geom_line() +
	facet_wrap(~ Hypothesis) +
	scale_y_continuous(trans='logit')

ggplot(smat, aes(x=Margin, y=power, col=factor(N))) +
	geom_line() +
	facet_wrap(~ Hypothesis) +
	scale_x_continuous(trans='log10')

ggplot(smat %>% filter(power < 0.99), aes(x=Margin, y=power, col=factor(N))) +
	geom_line() +
	facet_wrap(~ Hypothesis) +
	scale_x_continuous(trans='log10') +
	scale_y_continuous(trans='logit')

cmat <- smat %>%
	group_by(N, Margin) %>%
	summarise(power = min(power))
cmat %>%
	spread(Margin, power) %>%
	print(n=Inf)

library('mgcv')
model <- gam(plogis(power) ~ s(N) + s(Margin), data=cmat)
plot(model, pages=1)
cmat$pred <- predict(model)
plot(qlogis(cmat$pred), cmat$power)


model <- lm(plogis(power) ~ poly(N, 5) + poly(Margin, 10), data=cmat)
cmat$pred <- predict(model)
plot(qlogis(cmat$pred), cmat$power); abline(0,1)

