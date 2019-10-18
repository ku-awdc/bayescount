function(input, output, session) {

	observe({
		
		sums <- c(input$Sum_1, input$Sum_2)
		paired <- input$paired
		if(paired){
			Ns <- c(input$N, input$N)
			cor <- input$cor
		}else{
			Ns <- c(input$N_1, input$N_2)
			cor <- 0
		}
		ks <- c(input$k_1, input$k_2)
		T_I <- input$T_I
		T_A <- input$T_A
		
		results <- efficacy_typologies(sum=sums, N=Ns, k = ks, cor = cor, paired = paired, T_I = T_I, T_A = T_A, R = c(1, 1), S = c(1, 1),
			method = "all", tail = 0.025, bnb_priors = c(0, 0), use_delta = NA, beta_iters = 10^4, binomial_priors = c(1, 1), binomial_cl_adj = 0.2)
		
		
		results$LCI <- ifelse(is.na(results$LCI), "", paste0(round(100*results$LCI,2), "%"))
		results$UCI <- ifelse(is.na(results$UCI), "", paste0(round(100*results$UCI,2), "%"))
		pp <- function(x) ifelse(x < 0.001, rep("<0.001",length(x)), format(c(round(x,3),0.001))[1:length(x)])
		results$pI <- ifelse(is.na(results$pI), "", pp(results$pI))
		results$pA <- ifelse(is.na(results$pA), "", pp(results$pA))
		
		results$Classification_Denwood <- results$Classification
		results$Classification <- NULL
		results$Classification_Coles <- results$Typology
		levels(results$Classification_Coles) <- list("Resistance Present" = c("1ab", "1c", "2a", "2b"), "Resistance Suspected" = c("2c", "3", "4a"), "Susceptible" = c("4bc"), "Error" = c("100%red", "<0%red", "error"))
		
		output$efficacy <- renderText(paste0("<h4 align='left'>Resulting efficacy is ", round(100*(1-(sums[2]/Ns[2])/(sums[1]/Ns[1])),2), "%, with classifications:</p>"))
		output$tbl <- renderTable(results) 
		
	})
	
    output$footer <- renderText(paste0("<p align='center'>Typologies exploration tool for ERR / FECRT data by Matthew Denwood<br>", footeraddtext, "<a href='http://www.fecrt.com/framework/', target='_blank'>Click here for more information (opens in a new window)</a></p>"))
	
}
