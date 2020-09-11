#' Title
#'
#' @param Nrange
#' @param EPGrange
#' @param EDT
#' @param target
#' @param margin
#' @param kpre
#' @param kpost
#' @param cor
#'
#' @return
#'
#' @export
fecrt_sample_size <- function(Nrange = c(5, 100), EPGrange = c(50, 1000), EDT = c(1, 2.5, 5, 10, 25, 50), target = 0.99, margin = c(0.04, 0.09), k1 = 1, k2 = 0.6, cor = 0.3, Nlength=50, mulength=20, EPGlength=100, iters=1000L){

	Ns <- unique(round(10^seq(log10(Nrange[1]), log10(Nrange[2]), length=Nlength)))

	murange <- EPGrange/rev(range(EDT))
	mu <- round(10^(seq(log10(max(5,murange[1])), log10(min(1000,murange[2])), length=mulength)), 2)

	## TODO: Make into arguments
	conjugate_priors <- c(1,1)
	delta <- 2L  # 0=never, 1=unless_fails, 2=always
	beta_iters <- 1000L
	approx <- 0L  # 0=never, 1=if_necessary, 2=always
	tail <- 0.05
	useml <- 0L

	adjk <- adjust_k(k1, k2, cor)
	kc <- adjk$correlated_k
	## Note: k1 and k2 are total over-dispersions
	kpre <- k1
	kpost <- k2

	pmat <- RCPP_power_matrix_paired(Ns, margin, mu, target, kpre, kpost, kc, iters, conjugate_priors, delta, beta_iters, approx, tail, useml)

	if(any(pmat$Typology=="error")){
		cat("Note:", sum(pmat$Typology=="error"), "errors\n")
	}

	smat <- pmat %>%
		filter(Typology!="error") %>%
		group_by(N, Margin, Mean, Hypothesis, Reduction) %>%
		summarise(power = mean(RejectH0), meanred = mean(ObsReduction), varred = var(ObsReduction))
	cmat <- smat %>%
		group_by(N, Margin, Mean) %>%
		summarise(power = min(power))

	## Remove any interpolations outside the range of obs:
	check <- smat %>%
		group_by(Margin, Mean) %>%
		summarise(Nunder = sum(power < 0.8), Nover = sum(power > 0.8), .groups='drop')
	## NB: if the power is always > 0.8 then just set N to min(Nrange)
	##     otherwise if power is always < 0.8 then set N to unknown

	cj <- check %>%
		filter(Nunder>=1, Nover>=1)

	nmu <- cmat %>%
		right_join(cj, by=c('Margin', 'Mean')) %>%
		ungroup() %>%
		mutate(MarginMean = str_c(Margin, "_", Mean)) %>%
		split(.$MarginMean) %>%
		map_df( ~
							with(.x, smooth.spline(power, log10(N), tol=1e-6 * IQR(Ns))) %>%
							predict(0.8) %>%
							as_tibble() %>%
							select(power=x, N=y) %>%
							bind_cols(.x %>% select(Margin, Mean) %>% slice(1))
		) %>%
		mutate(N = 10^N) %>%
		arrange(Margin, Mean)

	## Fake data set:
	fakedata <- expand_grid(EPG=seq(EPGrange[1], EPGrange[2], length.out=EPGlength), EDT = EDT) %>%
		mutate(Mean = EPG/EDT)

	fmu <- nmu %>%
		ungroup() %>%
		split(.$Margin) %>%
		map_df( ~ with(.x, smooth.spline(log10(Mean), log10(N))) %>%
							predict(log10(fakedata$Mean)) %>%
							as_tibble() %>%
							select(N=y) %>%
							mutate(N = 10^N) %>%
							bind_cols(fakedata)# %>%
							#filter(Mean >= min(.x$Mean), Mean <= max(.x$Mean))
						, .id="Margin")# %>%
		#filter(Mean >= murange[1], Mean <= murange[2], N >= min(Ns), N <= max(Ns))
		#filter(Mean >= min(mu), Mean <= max(mu))

	return(list(smat = smat, nmu = nmu, fmu = fmu))

	smu <- fmu %>%
		mutate(N = ceiling(N))


	ggplot(fmu, aes(x=EPG, y=N, col=factor(EDT))) +
		geom_line() +
		scale_y_continuous(trans='log10') +
		scale_x_continuous(trans='log10') +
		ggtitle(combos %>% slice(r) %>% {function(x) paste(x[c(1:3, 5)], collapse=" - ")}()) +
		facet_wrap(~Margin) +
		coord_cartesian(ylim=c(5, 50))

	ggplot(fmu, aes(x=Mean, y=N)) +
		geom_line() +
		scale_y_continuous(trans='log10') +
		#	scale_x_continuous(trans='log10') +
		ggtitle(combos %>% slice(r) %>% {function(x) paste(x[c(1:3, 5)], collapse=" - ")}()) +
		facet_wrap(~Margin) +
		coord_cartesian(ylim=c(5, 50), xlim=c(0, 10))


	ggplot(nmu, aes(x=Mean, y=N)) +
		geom_point() +
		geom_line() +
		#	geom_smooth(se=FALSE) +
		scale_y_continuous(trans='log10') +
		scale_x_continuous(trans='log10', breaks=c(10,20,30,50,100,200,300,500,1000)) +
		ggtitle(combos %>% slice(r) %>% {function(x) paste(x[c(1:3, 5)], collapse=" - ")}()) +
		facet_wrap(~Margin)



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


}
