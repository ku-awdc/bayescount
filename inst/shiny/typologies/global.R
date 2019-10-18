## The required packages must be supplied like this on the second line of the file as it is used to check the packages are installed before launching:
packages <- c("bayescount","shiny","shinythemes")

# Needs to be explicitly in here as library etc calls for deployment to shinyapps.io:
library("bayescount")
library("shiny")
library("shinythemes")
# TODO: handle this more nicely

# Load the packages:
if(!all(sapply(packages, require, character.only=TRUE))) stop("One or more required package missing")

# Other global settings and options:
options(stringsAsFactors=FALSE)

headscript <- ""
footeraddtext <- ""
