## The required packages must be supplied like this on the second line of the file as it is used to check the packages are installed before launching:
packages <- c("bayescount","shiny","shinythemes","rhandsontable")

# Needs to be explicitly in here as library etc calls for deployment to shinyapps.io:
library("bayescount")
library("shiny")
library("shinythemes")
library("rhandsontable")
library("stringr")
library("ggplot2")
library("dplyr")
library("purrr")
theme_set(theme_light())
# TODO: handle this more nicely

# Load the packages:
if(!all(sapply(packages, require, character.only=TRUE))) stop("One or more required package missing")

# Other global settings and options:
options(stringsAsFactors=FALSE)

### For testing:
testing <- FALSE
# testing <- TRUE
###

parasitology <- FALSE

headscript <- "Hello"
footeraddtext <- "Footer"

library("markdown")
colwidth <- 6L

species_choices <- list(
	`*SELECT*` = "INVALID",
	Ruminants = list(`Cattle - Nematodes` = "cattle", `Sheep - Nematodes` = "sheep", `Goats - Nematodes` = "goats"),
	Equine = list(`Adults - Cyathostomins` = "cyath", `Foals - Parascaris` = "foals"),
	`Swine` = list(`Oesophagostomum` = "pig")
)

host_choices <- list(
	`*SELECT*` = "INVALID",
	Ruminants = list(`Cattle` = "cattle", `Sheep` = "sheep", `Goats` = "goats"),
	Equine = list(`Horses (adults)` = "horse_adult", `Horses (foals)` = "Horses (foals)", `Donkeys (adults)` = "donk_adult", `Donkeys (foals)` = "donk_foal"),
	`Swine` = "pigs",
	`Other` = "other"
)

repnull <- function(x) rep(NA_integer_, if(is.null(x)) 0L else x)
changed <- function(x,input,settings) any(sapply(x, function(y) !is.null(input[[y]]) && !identical(input[[y]],settings[[y]])))
