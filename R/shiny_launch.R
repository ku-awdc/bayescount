shiny_launch <- function(...){
	
	# Might have different shiny apps within here at some point?
	
	ad <- system.file("shiny_typologies", package = "bayescount")
	if(ad == ""){
		stop("Package directory not found - is the bayescount package installed?", call. = FALSE)
	}
	  
	shiny::runApp(appDir=ad, ...)
	
}
