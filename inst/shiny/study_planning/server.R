function(input, output, session) {

	# res <- fecrt_sample_size(EPGrange=50*20, EDT=50, mulength=1, margin=0.04, iters=10000L)
	# # EPGrange=epgrange, mulength=length(epgrange), EDT=1, target=using$Target[1]/100, margin=unique(using$Delta/100), k1=using$k1[1], k2=using$k2[1], cor=using$cor[1], iters=iters
	# ggplot(res$smat, aes(x=N, y=power, col=Hypothesis)) +
	# 	geom_line() +
	# 	geom_point() +
	# 	scale_x_continuous(trans="log10")

	rv <- reactiveValues(lowerThreshold = 95, results_table = data.frame(N=numeric(0), power=numeric(0), Hypothesis=character(0)), results_text = "", parameter_feedback="", plot=ggplot(), status=0)

	observe({
		input$targetEfficacy
		input$delta

		rv$lowerThreshold <- input$targetEfficacy - input$delta
	})

	observeEvent(input$calculate, {

		rv$parameter_feedback <- "Calculation unavailable due to a software update - please try again in a few hours"

		if(FALSE){
		## Validate inputs:
		if(input$species=="INVALID"){
			rv$parameter_feedback <- "Please select a host/parasite species and click Calculate again"
			validate(rv$parameter_feedback)
		}
		if(input$multfact==0){
			rv$parameter_feedback <- "Multiplication factor must be greater than zero"
			validate(rv$parameter_feedback)
		}
		guideline <- input$species
		k1 <- switch(guideline,
			"cattle" = 1.1,
			"sheep" = 1.4,
			"goats" = 3.8,
			"cyath" = 1.3,
			"foals" = 1.1,
			"pig" = 0.8,
			NA_real_
		)
		k2 <- switch(guideline,
			"cattle" = 0.4,
			"sheep" = 0.8,
			"goats" = 2.4,
			"cyath" = 0.5,
			"foals" = 1.0,
			"pig" = 1.2,
			NA_real_
		)
		kc <- switch(guideline,
			"cattle" = 0.1,
			"sheep" = 0.5,
			"goats" = 0.5,
			"cyath" = 0.5,
			"foals" = 0.2,
			"pig" = 0.4,
			NA_real_
		)
		stopifnot(!is.na(k1), !is.na(k2), !is.na(kc))

		rv$parameter_feedback <- "Calculation unavailable due to a software update - please try again in a few hours"
		validate(rv$parameter_feedback)


		withProgress(message = "Calculating...", value= 0, {

			setProgress(1/4)
			res <- fecrt_sample_size(Nrange = c(5, input$maxN), Nlength=100L, EPGrange=input$meanEPG, EDT=input$multfact, mulength=1, target=input$targetEfficacy/100, margin=input$delta/100, iters=10000L, k1=k1, k2=k2, cor=kc)

			setProgress(2/4)
			rv$results_table <- res
			res$smat |>
				filter(power>=0.8) |>
				group_by(Hypothesis) |>
				arrange(N) |>
				slice(1) ->
				reqss

			if(nrow(reqss)<2){
				rv$parameter_feedback <- "A suitable sample size could not be found - try increasing the maximum sample size"
				validate(rv$parameter_feedback)
			}
			stopifnot(nrow(reqss)==2)

			stopifnot(all(c("Efficacy","Resistance") %in% reqss$Hypothesis))
			rv$results_text <- str_c(
				"Required sample size is ", max(reqss$N), "<br>",
				"(", reqss |> filter(Hypothesis=="Resistance") |> pull(N), " for the test for resistance and ", reqss |> filter(Hypothesis=="Efficacy") |> pull(N), " for the test for Susceptibility)"
				)

			setProgress(3/4)
			rv$results_plot <- 	ggplot(res$smat, aes(x=N, y=power, col=Hypothesis)) +
				geom_hline(yintercept=0.8, lty="dashed") +
				geom_line() +
				geom_point() +
				theme(legend.pos='bottom')# + scale_x_continuous(trans="log10")

			setProgress(4/4)

		})

		rv$status <- 1
		}
	})

	fluidPage({
		output$lowerThreshold <- renderText(str_c("<h5>Lower efficacy target = ", rv$lowerThreshold, "%</h5>"))
	})

	fluidPage({
		output$parameter_feedback <- renderText(rv$parameter_feedback)
		output$results_text <- renderText(rv$results_text)
		output$results_plot <- renderPlot(print(rv$results_plot))
	})

	output$status <- renderText(rv$status)
	outputOptions(output, 'status', suspendWhenHidden=FALSE)

}
