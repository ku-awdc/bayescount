# TODO:  efficacy_simulate function to simulate paired/unpaired data


estimate_k <- function(data_1, data_2, paired){
	# TODO: use findtheta in C++ code:
	#	summarise(N=n(), mu=mean(PreTx), ko1 = findtheta(PreTx), ko2 = findtheta(PostTx), cor=cor(PreTx, PostTx)) %>%
	#	mutate(ko2 = ifelse(is.na(cor), ko1, ko2), cor = ifelse(is.na(cor), 0, cor))
	
	# Then use adjust_k if paired
}

adjust_k <- function(total_k1, total_k2, cor){
	
	checksingleprob(cor)
	checksingleposdouble(total_k1)
	checksingleposdouble(total_k2)
	
	kc <- sqrt(total_k1*total_k2) / cor
	kp1 <- ((kc+1)*total_k1) / (kc-total_k1)
	kp2 <- ((kc+1)*total_k2) / (kc-total_k2)
	
	return(list(total_k = c(total_k1, total_k2), correlated_k = kc, uncorrelated_k = c(kp1, kp2)))
}
