## The required packages must be supplied like this on the second line of the file as it is used to check the packages are installed before launching:
packages <- c("shiny","shinythemes","rhandsontable")

# Load the packages:
if(!all(sapply(packages, require, character.only=TRUE))) stop("One or more required package missing")

# Other global settings and options:
options(stringsAsFactors=FALSE)

### For testing:
testing <- FALSE
# testing <- TRUE
###

parasitology <- FALSE

headscript <- ""
footeraddtext <- ""
