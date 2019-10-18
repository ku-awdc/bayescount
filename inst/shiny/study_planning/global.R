## The required packages must be supplied like this on the second line of the file as it is used to check the packages are installed before launching:
packages <- c("bayescount","shiny","shinythemes","dplyr","ggplot2")

# Load the packages:
if(!all(sapply(packages, require, character.only=TRUE))) stop("One or more required package missing")

# Other global settings and options:
options(stringsAsFactors=FALSE)

# Permanent settings:
citers <- 10^4



# Default min, max and step for log sliders:
logslidvals <- c(-2,2,0.025)

# logifySlider javascript function
JS.logify <-
"
// function to logify a sliderInput
function logifySlider (sliderId, sci = false) {
  if (sci) {
    // scientific style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return ('10<sup>'+num+'</sup>'); }
    })
  } else {
    // regular number style
    $('#'+sliderId).data('ionRangeSlider').update({
      'prettify': function (num) { return (Math.round(Math.pow(10, num)*100)/100); }
    })
  }
}"

# call logifySlider for each relevant sliderInput
JS.onload <-
"
// execute upon document loading
$(document).ready(function() {
  // wait a few ms to allow other scripts to execute
  setTimeout(function() {
    // include call for each slider
    logifySlider('k1', sci = false)
    logifySlider('k2', sci = false)
  }, 5)})
"



headscript <- ""
footeraddtext <- ""
