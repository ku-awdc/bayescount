# Should not be called in this code, but bayescount needs to have been
# installed this way before deployment:
# devtools::install_github('mdenwood/bayescount')

library('shiny')
library('bayescount')
library('dplyr')
library('tibble')
library('tidyr')
library('ggplot2')

options(stringsAsFactors=FALSE)

# Permanent settings:
citers <- 10^3

# Read settings from a website so they can be dynamically altered:
load(url('https://ku-awdc.github.io/rsc/shinyresults.Rdata'))

pubdate <- pardate

# Not currently using:
logslidvals <- c(-2,2,0.025)


getsimres <- function(simpars){

	tsp <- as.data.frame(simpars) %>%
		mutate(ti=target, ta=(target-delta), ko1=k1, ko2=k2, kd=k2/k1) %>%
		mutate(kc = sqrt(ko1*ko2) / cor, ku1 = ((kc+1)*ko1)/(kc-ko1), ku2 = ((kc+1)*ko2)/(kc-ko2))

	# Get lower and upper as function of ti and ta:
	upper <- 1 - max(0, tsp$ta-0.05)
	lower <- 1 - min(1, tsp$ti+0.05)

	# Make by a sensible number to give around 25 values depending on upper-lower:
	by <- seq(0,upper-lower,length.out=25)[2]
	choices <- c(0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05)
	by <- choices[which.min(abs(by-choices))]

	rvs <- seq(lower,upper,by=by)
	if(0 %in% rvs)
		rvs <- unique(c(0, choices[1:6], rvs))
	rvs <- unique(c(rvs, (1-tsp$ta)+c(-by/2,0,by/2), (1-tsp$ti)+c(-by/2,0,by/2)))

	args <- as.list((tsp %>% select(N, mu, k1, k2, kc, ta, ti)))
	args <- c(args, list(rvs=rvs, iters=citers, approx=1, priors=c(0,0), useml=TRUE, plot=FALSE))

	res <- do.call('powersim_paired', args)
	
    res$Classification <- factor(res$Classification)
    levels(res$Classification) <- classifications[["Typologies"]]
    stopifnot(all(!is.na(res$Classification)))
  
    res$Method[res$Method=='Levecke'] <- 'Gamma'
    res$Method[res$Method=='MLE'] <- 'Asymptotic'
  
    typres <- res %>%
      mutate(Efficacy = 1-Reduction) %>%
      group_by(Method, Efficacy, Classification) %>%
      summarise(Number = n()) %>%
      # Need to join with full dataset:
      full_join( crossing(Efficacy=1-args$rvs, Method=unique(res$Method), Classification=factor(levels(res$Classification), levels=levels(res$Classification))) , by = c("Method", "Efficacy", "Classification")) %>%
      # Any blank numbers are zeros:
      mutate(Number = ifelse(is.na(Number), 0, Number)) %>%
      group_by(Method, Efficacy) %>%
      mutate(Proportion = Number / sum(Number)) %>%
      ungroup()
 	
	return(list(simpars=simpars, simres=typres))
}

inits <- simparameters %>% slice(1)
initsimpars <- list(N=inits$N, mu=inits$mu, target=inits$Target/100, delta=inits$Delta/100, pthresh=0.025, k1=inits$k1, k2=inits$k2, cor=inits$cor)


function(input, output, session) {
	
	rv <- reactiveValues(simpars=initsimpars, storedres=simresults[[1]], calculating=FALSE)
	
	observe({
		preset <- input$preset
		
		m <- which(preset==presets)
		stopifnot(length(m)==1)
		sp <- simparameters %>% slice(m)
		updateNumericInput(session, "N", value=sp$N)
		updateNumericInput(session, "mu", value=sp$mu)
		updateNumericInput(session, "target", value=sp$Target)
		updateNumericInput(session, "delta", value=sp$Delta)
		updateNumericInput(session, "pthresh", value=0.025)
		#updateSliderInput(session, "k1", value=log10(sp$k1))
		#updateSliderInput(session, "k2", value=log10(sp$k2))
		#updateSliderInput(session, "cor", value=sp$cor)
		updateNumericInput(session, "k1", value=sp$k1)
		updateNumericInput(session, "k2", value=sp$k2)
		updateNumericInput(session, "cor", value=sp$cor)
		
	})

	observeEvent(input$simulate, {
		
		# TODO: work out how to set a message saying 'calculating'
		rv$calculating <- TRUE
		
		#simpars <- list(StudyType=input$type, N=input$N, mu=input$mu, target=input$target/100, delta=input$delta/100, pthresh=input$pthresh, kpreset=input$kpreset, k1=10^input$k1, k2=10^input$k2, cor=input$cor)
		simpars <- list(N=input$N, mu=input$mu, target=input$target/100, delta=input$delta/100, pthresh=input$pthresh, k1=input$k1, k2=input$k2, cor=input$cor)
		
		# Find an existing set if possible:
		matchpars <- simparameters %>% filter(N==input$N, mu==input$mu, Target==input$target, Delta==input$delta, k1==input$k1, k2==input$k2, cor==input$cor)
		
		stopifnot(nrow(matchpars)<=1)
		if(input$pthresh==0.025 && nrow(matchpars)==1){
			rv$storedres <- simresults[[matchpars$Set]]
		}else{
			newres <- getsimres(simpars)
			rv$storedres <- newres$simres
		}
		
		rv$simpars <- simpars
		
		rv$calculating <- FALSE
		
	})
	
	
	fluidPage(
		output$classifications <- renderText(frameworks[[input$statframe]])
	)
 	
	observe({
		statframe <- input$statframe
		
		wd <- 900
		he <- 600
		
		if(input$statframe=="Denwood"){
			output$typologies <- renderImage({
				list(src="typologies_denwood.png", contentType = 'image/png', width = wd, height = he, alt = "Typologies image")
			}, deleteFile=FALSE)
		}else if(input$statframe=="Typologies"){
			output$typologies <- renderImage({
				list(src="typologies_raw.png", contentType = 'image/png', width = wd, height = he, alt = "Typologies image")
			}, deleteFile=FALSE)
		}else if(input$statframe=="Kaplan"){
			output$typologies <- renderImage({
				list(src="typologies_kaplan.png", contentType = 'image/png', width = wd, height = he, alt = "Typologies image")
			}, deleteFile=FALSE)
		}else if(input$statframe=="Coles"){
			output$typologies <- renderImage({
				list(src="typologies_coles.png", contentType = 'image/png', width = wd, height = he, alt = "Typologies image")
			}, deleteFile=FALSE)
		}else{
			stop('Unrecognised statistical framework')
		}
		
	})
	
	
	output$noninc_plot <- renderPlot({
		
		typres <- rv$storedres
		levels(typres$Classification) <- classifications[[input$statframe]]
		
		methodsusing <- input$statmethod
		
		typres <- typres %>%
			filter(Method %in% methodsusing) %>%
		    group_by(Method, Efficacy, Classification) %>%
		    summarise(Number = sum(Number), Proportion = sum(Proportion))
			

		theme_set(theme_light())
		
		ti <- rv$simpars$target
		ta <- ti - rv$simpars$delta
		incres <- typres %>%
	    	filter(Classification %in% inconclusives[[input$statframe]] ) %>%
		    group_by(Method, Efficacy) %>%
		    summarise(Number = sum(Number), Proportion = sum(Proportion)) %>%
			mutate(region = ifelse(Efficacy <= ta, 'one', ifelse(Efficacy < ti, 'two', 'three')))
	
		ggplot(incres, aes(x=Efficacy, y=1-Proportion, group=region)) +
			geom_rect(xmin=ta, xmax=ti, ymin=0, ymax=1, fill='grey80', col='grey80', alpha=0.5) +
		#	geom_line(stat='smooth', method='loess', formula=y~x, se=FALSE, alpha=0.5, lwd=2) +
		#	geom_line(stat='smooth', method='lm', formula=y~poly(x,2), se=FALSE, alpha=0.5, lwd=2) +
			geom_line(alpha=0.5, lwd=2) +
			geom_point() +
			geom_vline(xintercept=ti, lty='dashed') +
			geom_vline(xintercept=ta, lty='dotted') +
			geom_hline(yintercept=0.8, lty='solid') +
			facet_grid(rows="Method") +
			theme(legend.pos='none') +
			scale_y_continuous(breaks=seq(0,1,by=0.2), limits=c(0,1))
	})
	
	output$typologies_plot <- renderPlot({

		typres <- rv$storedres
		levels(typres$Classification) <- classifications[[input$statframe]]
		
		methodsusing <- input$statmethod
		
		typres <- typres %>%
			filter(Method %in% methodsusing) %>%
		    group_by(Method, Efficacy, Classification) %>%
		    summarise(Number = sum(Number), Proportion = sum(Proportion))
			

		theme_set(theme_light())
		
		ti <- rv$simpars$target
		ta <- ti - rv$simpars$delta
		ggplot(typres, aes(x=Efficacy, y=Proportion, col=Classification)) +
			geom_rect(xmin=ta, xmax=ti, ymin=0, ymax=1, fill='grey80', col='grey80', alpha=0.5) +
		#	geom_line(stat='smooth', method='lm', formula=y~poly(x,10), se=FALSE, alpha=0.5, lwd=2) +
		#	geom_line(stat='smooth', method='loess', formula=y~x, se=FALSE, alpha=0.5, lwd=2) +
			geom_line(alpha=0.5, lwd=2) +
			geom_point() +
			geom_vline(xintercept=ti, lty='dashed') +
			geom_vline(xintercept=ta, lty='dotted') +
			geom_hline(yintercept=0.025, lty='solid') +
			facet_grid(rows="Method") +
			theme(legend.pos='bottom')
		
	})
	
	output$parameters <- renderText(paste0('<p>Showing results for the ', if(input$customise) '(customised) ', 'guidelines for ', input$preset, '</p><p>[Parameters: N = ', rv$simpars$N, ', mu = ', rv$simpars$mu, ', k1 = ', rv$simpars$k1, ', k2 = ', rv$simpars$k2, ', correlation = ', rv$simpars$cor, ', target = ', rv$simpars$target, ', delta = ', rv$simpars$delta, ', alpha = ', rv$simpars$pthresh, ']</p>Probabilities calculated using the ', input$statmethod, ' method and ', if(input$statframe=='Typologies') 'raw typologies' else paste0(input$statframe, ' et al framework'), '</p>'))
	
	output$footer <- renderText(paste0('<p align="center">Power calculation tool for FECRT data by Matthew Denwood (last updated ', pubdate, ')<br><a href="http://www.fecrt.com/", target="_blank">Learn more (opens in a new window)</a></p>'))
}
