#' @name launch_shiny
#' @title Launch an interactive Shiny app contained within the bayescount package
#'
#' @description
#' This function provides a convinient way in which to launch one of the interactive shiny apps contained within this package, as detailed at http://www.fecrt.com/framework
#' 
#' @param appname the name of an app to launch (see the details section)
#'
#' @details
#' The apps that are currently available are:
#' - **data_analysis** : analyse an existing dataset
#' - **framework** : examine the implications of different interpretative frameworks for ERR / FECR
#' * **study_planning** : sample size calculations for a prospective study
#' 
#' @seealso \code{\link{bayescount}}
#' @export

launch_shiny <- function(appname, ...){
	
	ad <- system.file("shiny", package = "bayescount")
	if(ad == ""){
		stop("Package directory not found - is the bayescount package installed?", call. = FALSE)
	}
	apa <- list.files(ad)
	  
	if(missing(appname) || !is.character(appname) || length(appname)!=1 || is.na(appname) || ! appname %in% apa){
		stop("Please specify one of the following appnames to run:\n\t", paste(apa, collpase=' '), call. = FALSE)
	}
	
	if(appname == "typologies"){
		pn <- c("dplyr", "tibble", "tidyr", "ggplot2")

	}else if(appname == "data_analysis"){
		pn <- "rhandsontable"
		parasitology <- FALSE

	}else if(appname == "framework"){
		stop("not yet implemented")

	}else if(appname == "study_planning"){
		stop("not yet implemented")

	}else{
		warning("The requirements for the specified shiny app are not known")
	}
	
	if(!all(sapply(c('shiny', 'shinythemes', pn), requireNamespace, quietly=TRUE))){
		stop("Additional packages must be installed for the specified shiny app - you should run:", iptext(c('shiny','shinythemes'), pn), call. = FALSE)
	}
		
	shiny::runApp(appDir=file.path(ad, appname), ...)
	
}

iptext <- function(...) paste0('\n\tinstall.packages(c(', do.call('paste', c(lapply(list(...), function(x) paste0('"', paste(x, collapse='", "'), '"')), list(sep=', '))), '))\n')
