function(input, output, session) {
	
	rv <- reactiveValues(calculating=FALSE)
	
	observeEvent(input$simulate, {
		
		# TODO: work out how to set a message saying 'calculating'
		rv$calculating <- TRUE
		
		paired <- input$paired
		if(paired){
			Ns <- c(input$N, input$N)
			cor <- input$cor
		}else{
			Ns <- c(input$N_1, input$N_2)
			cor <- 0
		}
		mean <- input$mean
		ks <- c(input$k_1, input$k_2)
		T_I <- input$T_I
		T_A <- input$T_A
		
		simres <- efficacy_frequencies(r = c(T_I,T_A), paired = paired, T_I = T_I, T_A = T_A, N = Ns, R = c(1, 1), S = c(1, 1), mean = mean, k = ks, cor = cor,
			iterations = citers, tail = 0.025, bnb_priors = c(0, 0), use_delta = NA, beta_iters = 10^4, use_ml = TRUE, binomial_priors = c(1, 1), binomial_cl_adj = 0.2) %>%
			filter(Method==input$method) %>%
			group_by(Method, Efficacy, Classification) %>%
			summarise(Frequency = sum(Frequency), Proportion = sum(Proportion)) %>%
			ungroup()
		simres <- simres %>%
			full_join(expand.grid(Efficacy = c(T_I, T_A), Method=unique(simres$Method), Classification = factor(c("Reduced", "Inconclusive", "Borderline", "Adequate","Method_Failure"),levels=c("Reduced", "Inconclusive", "Borderline", "Adequate","Method_Failure")), stringsAsFactors=FALSE), by = c("Method", "Efficacy", "Classification"))
		simres$Frequency[is.na(simres$Frequency)] <- 0
		simres$Proportion[is.na(simres$Proportion)] <- 0		
			
		simres_a <- simres %>% filter(Efficacy == T_A)
		simres_i <- simres %>% filter(Efficacy == T_I)
		
		mkperc <- function(x) paste0(round(x*100, 1), '%')
		
		output$tbl_a <- renderTable(simres_a %>% mutate(Probability = mkperc(Proportion)) %>% select(-Method, -Efficacy, -Frequency, -Proportion))
		output$tbl_i <- renderTable(simres_i %>% mutate(Probability = mkperc(Proportion)) %>% select(-Method, -Efficacy, -Frequency, -Proportion))

		output$tta <- renderText(paste0("<h4>Inferiority test with efficacy = ", T_A*100, "%:</p>"))
		output$tti <- renderText(paste0("<h4>Non-inferioirty test with efficacy = ", T_I*100, "%:</p>"))
		
		power_a <- simres_a %>% filter(Classification %in% c("Reduced","Borderline")) %>% summarise(prob=sum(Proportion)) %>% pull(prob)
		error_a <- simres_a %>% filter(Classification %in% c("Adequate")) %>% summarise(prob=sum(Proportion)) %>% pull(prob)		
		output$bta <- renderText(paste0("<h5>Power = ", mkperc(power_a), "</h5><h5>Type-1 error = ", mkperc(error_a), "</h5>"))

		power_i <- simres_i %>% filter(Classification %in% c("Adequate","Borderline")) %>% summarise(prob=sum(Proportion)) %>% pull(prob)
		error_i <- simres_i %>% filter(Classification %in% c("Reduced")) %>% summarise(prob=sum(Proportion)) %>% pull(prob)		
		output$bti <- renderText(paste0("<h5>Power = ", mkperc(power_i), "</h5><h5>Type-1 error = ", mkperc(error_i), "</h5>"))
		
		if(power_a > 0.8 && power_i > 0.8){
			output$interp <- renderText("<h4 style='text-align:center; color:green; '>The sample size is large enough to achieve 80% power for both tests</span>")
		}else{
			output$interp <- renderText("<h4 style='text-align:center; color:red; '>The sample size is not large enough to achieve 80% power for both tests</span>")
		}
		
		rv$calculating <- FALSE
		
	})
	
	 output$footer <- renderText(paste0("<p align='center'>Sample size calculation tool for ERR / FECRT data by Matthew Denwood<br>", footeraddtext, "<a href='http://www.fecrt.com/framework/', target='_blank'>Click here for more information (opens in a new window)</a></p>"))
}
