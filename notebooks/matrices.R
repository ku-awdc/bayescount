library('tidyverse')

### WAAVP guidelines matrixes


numbers <- read_csv("/Users/matthewdenwood/Documents/Research/Projects/WAAVP/Power paper/analyses/median_numbers.csv")

numbers %>%
	mutate(kr = k2/k1) %>%
	summarise(kr = mean(kr, na.rm=TRUE), correlation = mean(correlation, na.rm=TRUE))

# Fill in numbers for Ascaris and Trichuris:
numbers <- numbers %>%
	mutate(k2 = case_when(
		!is.na(k2) ~ k2,
		TRUE ~ k1 * 0.9
	)) %>%
	mutate(correlation = case_when(
		!is.na(correlation) ~ correlation,
		TRUE ~ 0.3
	)) %>%
	rowid_to_column() %>%
	split(.$rowid) %>%
	map_df( ~ .x %>% mutate(kc = bayescount:::adjust_k(k1, k2, correlation)$correlated_k))
numbers

combos <- tribble(~Host, ~Parasite, ~Drug, ~Target,
									"Calf", "Nematodes", "All", 99,
									"Cattle", "Nematodes", "All", 99,
									"Equine", "Nematodes", "AVM", 99.9,
									"Equine", "Nematodes", "Pyrantel", 98,
									"Pig", "Ascaris", "IVM", 97,
									"Pig", "Ascaris", "BZ", 98,
									"Pig", "Oesophagostomum", "IVM", 95,
									"Pig", "Oesophagostomum", "BZ", 98,
									"Pig", "Trichuris", "All", 60,
									"Sheep", "Nematodes", "All", 99
									) %>%
	expand_grid(mu = c(10.0, 50.0)) %>%
	left_join(numbers, by = c("Host", "Parasite")) %>%
	select(-rowid)

combos <- combos %>%
	filter(k1 > 0.5)

Ns <- seq(10, 30, by=1)
margin <- seq(0.01, 0.1, by=0.005)
margin <- seq(0.01, 0.1, by=0.01)
iters <- 5000L
conjugate_priors <- c(1,1)
delta <- 1L
beta_iters <- 1000L
approx <- 1L
tail <- 0.05
useml <- FALSE

theme_set(theme_light())
pdf("all_graphs.pdf")
for(r in seq_len(nrow(combos))){

	pmat <- bayescount:::RCPP_power_matrix_paired(Ns, margin, combos$Target[r]/100.0, combos$mu[r], combos$k1[r], combos$k2[r], combos$kc[r], iters, conjugate_priors, delta, beta_iters, approx, tail, useml)


smat <- pmat %>%
	group_by(N, Margin, Hypothesis, Reduction) %>%
	summarise(power = mean(RejectH0), meanred = mean(ObsReduction), varred = var(ObsReduction))
cmat <- smat %>%
	group_by(N, Margin) %>%
	summarise(power = min(power))
cmat %>%
	spread(Margin, power) %>%
	print(n=Inf)


pt <- ggplot(smat %>% mutate(Margin = factor(Margin, levels=rev(margin))), aes(x=N, y=power, col=Margin)) +
	#	geom_point() +
	geom_smooth(se=FALSE) +
	facet_wrap(~ Hypothesis) +
	scale_y_continuous(breaks=seq(0,1,by=0.05)) +
	coord_cartesian(ylim=c(0.5, 1)) +
	geom_hline(yintercept=0.8) +
	ggtitle(combos %>% slice(r) %>% {function(x) paste(x[c(1:3, 5)], collapse=" - ")}())

print(pt)
print(r)

}
dev.off()


ggplot(cmat, aes(x=factor(N), y=factor(Margin), fill=power)) +
	geom_tile() +
	scale_fill_stepsn(n.breaks=10, colors=c(rep('red',7), 'orange', 'green', 'blue'))

smat

ggplot(smat, aes(x=Margin, y=power, col=factor(N))) +
	geom_smooth() +
	facet_wrap(~ Hypothesis)
+
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

