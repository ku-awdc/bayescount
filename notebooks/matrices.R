library('tidyverse')
library('bayescount')

theme_set(theme_light())

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
	rowid_to_column()
numbers

combos <- tribble(~Host, ~Parasite, ~Drug, ~Target, ~Research, ~Clinical,
									"Calf", "Nematodes", "", 99, 0.04, 0.09,
									"Equine", "Nematodes", "AVM", 99.9, 0.019, 0.049,
									"Equine", "Nematodes", "Pyrantel", 98, 0.06, 0.08,
									"Pig", "Ascaris", "IVM", 97, 0.07, 0.15,
									"Pig", "Ascaris", "BZ", 98, 0.08, 0.18,
									"Pig", "Oesophagostomum", "IVM", 95, 0.07, 0.12,
									"Pig", "Oesophagostomum", "BZ", 98, 0.03, 0.08,
									"Sheep", "Nematodes", "All", 99, 0.02, 0.04
) %>%
	left_join(numbers, by = c("Host", "Parasite")) %>%
	select(-rowid)

theme_set(theme_light())
pdf("all_graphs.pdf", width=10, height=8)
for(r in 1:nrow(combos)){

	if(r %in% 4:5) next

	epgrange <- c(combos$muMin[r], combos$muQ3[r])
	margins <- c(combos$Research[r], combos$Clinical[r])
	cor <- combos$correlation[r]

	flist <- fecrt_sample_size(EPGrange=epgrange, target=combos$Target[r]/100, margin=margins, k1=combos$k1[r], k2=combos$k2[r], cor=cor, iters=5000L)
 # function(Nrange = c(5, 100), EPGrange = c(50, 1000), EDT = c(1, 2.5, 5, 10, 25, 50), target = 0.99, margin = c(0.04, 0.09), k1 = 1, k2 = 0.6, cor = 0.3, Nlength=50, mulength=20, EPGlength=100)


	fmu <- flist$fmu %>%
		mutate(N = pmax(N, 5))

	smu <- fmu %>%
		mutate(N = ceiling(N)) %>%
		mutate(EDT = factor(EDT, levels=rev(sort(unique(EDT))))) %>%
		mutate(Protocol = factor(Margin, levels=sort(unique(Margin)), labels=str_c(c("Research", "Clinical"), ": ", sort(unique(Margin)))))

	pt <- ggplot(smu, aes(x=EPG, y=N, col=EDT)) +
		geom_line() +
		# scale_y_continuous(trans='log10', breaks=c(5,7,10,14,20,30,50), minor_breaks=5:50) + coord_cartesian(ylim=c(5, 50)) +
		coord_cartesian(ylim=c(0, 50)) +
		# scale_x_continuous(trans='log10') +
		ggtitle(combos %>% slice(r) %>% {function(x) paste(x[1:3], collapse=" ")}()) +
		facet_wrap(~Protocol, ncol=1) +
		ylab("Sample Size") +
		xlab("Eggs per Gram (arithmetic mean)")

	print(pt)

}
dev.off()


###### OLDER STUFF

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
									"Equine", "Nematodes", "AVM", 99.9,
									"Equine", "Nematodes", "Pyrantel", 98,
									"Pig", "Ascaris", "IVM", 97,
									"Pig", "Ascaris", "BZ", 98,
									"Pig", "Oesophagostomum", "IVM", 95,
									"Pig", "Oesophagostomum", "BZ", 98,
									"Sheep", "Nematodes", "All", 99
									) %>%
	# expand_grid(EDT = c(1, 5, 10, 25, 50)) %>%
	expand_grid(mu = c(10, 50)) %>%
	left_join(numbers, by = c("Host", "Parasite")) %>%
	select(-rowid) %>%
	mutate(kpre = k1 / (1.0 - correlation)) %>%
	mutate(kpost = k2 / (1.0 - correlation))

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
useml <- 0L

theme_set(theme_light())
pdf("all_graphs.pdf")
for(r in seq_len(nrow(combos))){

	pmat <- bayescount:::RCPP_power_matrix_paired(Ns, margin, combos$mu[r], combos$Target[r]/100.0, combos$kpre[r], combos$kpost[r], combos$kc[r], iters, conjugate_priors, delta, beta_iters, approx, tail, useml)


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




## Multiple mu:
r <- 1

# Pick a range based on Q1 and Q3 of mean and possible edt of:
EDT <- c(1, 2.5, 5, 10, 25, 50)
#murange <- c(floor(combos$muQ1[r]/max(EDT)), ceiling(combos$muQ3[r]/min(EDT)))
epgrange <- c(50, 1000)
murange <- epgrange/rev(range(EDT))

Ns <- unique(round(10^seq(log10(5), log10(100), length=50)))
combos$mu[r]
mu <- exp(seq(log(murange[1]), log(murange[2]), length=20))
margin <- c(0.04, 0.09)
pmat <- bayescount:::RCPP_power_matrix_paired(Ns, margin, mu, combos$Target[r]/100.0, combos$kpre[r], combos$kpost[r], combos$kc[r], iters, conjugate_priors, delta, beta_iters, approx, tail, useml)


smat <- pmat %>%
	group_by(N, Margin, Mean, Hypothesis, Reduction) %>%
	summarise(power = mean(RejectH0), meanred = mean(ObsReduction), varred = var(ObsReduction))
cmat <- smat %>%
	group_by(N, Margin, Mean) %>%
	summarise(power = min(power))
cmat %>%
	spread(Margin, power) %>%
	print(n=Inf)


pt <- ggplot(smat %>% mutate(Mean = factor(Mean)), aes(x=N, y=power, col=Mean)) +
	#	geom_point() +
	geom_line() +
	facet_grid(Hypothesis ~ Margin) +
	scale_y_continuous(breaks=seq(0,1,by=0.05)) +
	scale_x_continuous(trans='log10') +
	coord_cartesian(ylim=c(0.5, 1)) +
	geom_hline(yintercept=0.8) +
	ggtitle(combos %>% slice(r) %>% {function(x) paste(x[c(1:3, 5)], collapse=" - ")}())

print(pt)

nmu <- cmat %>%
	ungroup() %>%
	mutate(MarginMean = str_c(Margin, "_", Mean)) %>%
	split(.$MarginMean) %>%
	map_df( ~
			 	with(.x, smooth.spline(power, log10(N))) %>%
			 		predict(0.8) %>%
			 		as_tibble() %>%
			 		select(power=x, N=y) %>%
			 		bind_cols(.x %>% select(Margin, Mean) %>% slice(1))
		) %>%
	mutate(N = 10^N) %>%
	arrange(Margin, Mean)

## Fake data set:
fakedata <- expand_grid(EPG=seq(epgrange[1], epgrange[2], by=10), EDT = EDT) %>%
	mutate(mu = EPG/EDT)


fmu <- nmu %>%
	ungroup() %>%
	split(.$Margin) %>%
	map_df( ~ with(.x, smooth.spline(log10(Mean), log10(N))) %>%
						predict(log10(fakedata$mu)) %>%
						as_tibble() %>%
						select(N=y) %>%
						bind_cols(fakedata)
	, .id="Margin") %>%
	mutate(N = 10^N) %>%
	filter(mu >= murange[1], mu <= murange[2], N >= min(Ns), N <= max(Ns))


ggplot(nmu, aes(x=Mean, y=N)) +
	geom_point() +
	geom_line() +
	#	geom_smooth(se=FALSE) +
	scale_y_continuous(trans='log10') +
	scale_x_continuous(trans='log10', breaks=c(10,20,30,50,100,200,300,500,1000)) +
	ggtitle(combos %>% slice(r) %>% {function(x) paste(x[c(1:3, 5)], collapse=" - ")}()) +
	facet_wrap(~Margin)

ggplot(fmu, aes(x=EPG, y=N, col=factor(EDT))) +
	geom_line() +
	scale_y_continuous(trans='log10') +
#	scale_x_continuous(trans='log10') +
	ggtitle(combos %>% slice(r) %>% {function(x) paste(x[c(1:3, 5)], collapse=" - ")}()) +
	facet_wrap(~Margin) +
	coord_cartesian(ylim=c(5, 50))



ss <- with(nmu, smooth.spline(log10(Mean), log10(N)))
plot(ss)
pp <- predict(ss, seq(0.5,2.5,0.1))
plot(pp$x, pp$y, type='l')
points(ss)
plot(10^pp$x, pp$y, type='l')
plot(10^pp$x, 10^pp$x*pp$y, type='l')


ggplot(fakedata, aes(x=EPG, y=10^N, col=factor(EDT))) +
	geom_line() +
#	scale_x_continuous(trans='log10') +
#	scale_y_continuous(trans='log10') +
	coord_cartesian(ylim=c(5, 50))

stop()


ggplot(cmat, aes(x=factor(N), y=factor(Margin), fill=power)) +
	geom_tile() +
	scale_fill_stepsn(n.breaks=10, colors=c(rep('red',7), 'orange', 'green', 'blue'))


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

