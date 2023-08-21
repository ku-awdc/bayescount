function(input, output, session) {

	res <- fecrt_sample_size(EPGrange=50*20, EDT=50, mulength=1, margin=0.04, iters=10000L)
	# EPGrange=epgrange, mulength=length(epgrange), EDT=1, target=using$Target[1]/100, margin=unique(using$Delta/100), k1=using$k1[1], k2=using$k2[1], cor=using$cor[1], iters=iters
	ggplot(res$smat, aes(x=N, y=power, col=Hypothesis)) +
		geom_line() +
		geom_point() +
		scale_x_continuous(trans="log10")

}
