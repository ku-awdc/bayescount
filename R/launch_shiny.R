#' @name launch_shiny
#' @title Launch an interactive Shiny app contained within the bayescount package
#'
#' @description
#' This function provides a convinient way in which to launch one of the interactive shiny apps contained within this package, as detailed at http://www.fecrt.com/framework
#'
#' @param appname the name of an app to launch (see the details section)
#' @param ... further arguments to be passed to \code{\link[shiny]{runApp}}
#' @details
#' The apps that are currently available are:
#' - **data_analysis** : analyse an existing dataset
#' - **typologies** : examine the typology resulting from a hypothetical dataset and implications of different interpretative frameworks for ERR / FECR
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

	if(missing(appname) || !is.character(appname) || length(appname)!=1 || is.na(appname)){
		stop("Please specify one of the following appnames to run:\n\t", paste(apa, collpase=' '), call. = FALSE)
	}
	appname <- apa[pmatch(appname, apa)]
	if(is.na(appname)){
		stop("Appname not recognised: please specify one of the following appnames to run:\n\t", paste(apa, collpase=' '), call. = FALSE)
	}
	
	# TODO: fix waavp_committee app to work with new bayescount
	if(appname=="waavp_committee")
		stop("That app is currently broken, sorry!")
	
	# TODO: replace this with the app-specific package requirements from top of global.R inside requested app:
	pn <- c('shiny', 'shinythemes', 'rhandsontable')
	# Also write a test to make sure they are all listed as imports/suggests by bayescount

	if(!all(sapply(pn, requireNamespace, quietly=TRUE))){
		stop("Additional packages must be installed for the specified shiny app - you should run:", iptext(pn), call. = FALSE)
	}

	shiny::runApp(appDir=file.path(ad, appname), ...)

}

iptext <- function(...) paste0('\n\tinstall.packages(c(', do.call('paste', c(lapply(list(...), function(x) paste0('"', paste(x, collapse='", "'), '"')), list(sep=', '))), '))\n')
