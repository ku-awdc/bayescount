library('shiny')
library('shinythemes')

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

# Read settings from a website so they can be dynamically altered:
load(url('https://ku-awdc.github.io/rsc/shinyresults.Rdata'))


# Default min, max and step for log sliders:
logslidvals <- c(-2,2,0.025)

fluidPage(
	# For log scale sliders:
    tags$head(tags$script(HTML(JS.logify))),
    tags$head(tags$script(HTML(JS.onload))),

	theme = shinytheme("cerulean"),
	
	hr(),
#	h3("Select Parameter Values", style="text-align:center; "),
#	hr(),
	
	tabsetPanel(
		tabPanel("Parameters",
			fluidRow(
				hr(),
				column(12,
					# div(style = "height:200px;background-color: yellow;", "Topleft"),
					h4("Choose Study Design Parameters", style="text-align:center; "),
					hr(),
					selectInput("preset", "Select pre-set study design parameters:", presets, selected=presets[1], width='100%'),
					checkboxInput("customise", "Customise parameters?", value=FALSE),

					hr(),
					actionButton("simulate", "Click to run the simulation with these parameters [this takes up to a few seconds to process]", width="100%"),
					hr()
			
				)
			),

			conditionalPanel(
		 		condition = "input.customise == true",
			
				fluidRow(
					column(4, 
			 

						 selectInput("type", "Study type", c("Simple Paired"), selected="Simple Paired", width='100%'),
						 # Also allow Simple Unpaired and more complex designs with replicates, pooling, unequal EDT etc?

						 numericInput('N', 'Sample size', value=10, min=5, step=1, width='100%'),
						 numericInput('mu', 'Mean pre-treatment count (EDT * mean EPG)', value=20, min=1, step=0.1, width='100%')
			 
					),
					column(4,
						numericInput('target', 'Target efficay (%)', min=0, max=100, value=95, step=0.1, width='100%'),
						numericInput('delta', 'Non-inferiority margin (% points)', min=0, max=25, value=5, step=0.1, width='100%'),
						numericInput('pthresh', 'Threshold for significance (p)', min=0, max=0.1, value=0.025, step=0.025, width='100%')
					),
					column(4,

						#sliderInput('k1', 'Pre-treatment k', min=logslidvals[1], max=logslidvals[2], value=0, step=logslidvals[3], width='100%'),
						#sliderInput('k2', 'Post-treatment k', min=logslidvals[1], max=logslidvals[2], value=-0.30, step=logslidvals[3], width='100%'),
						#sliderInput('cor', 'Correlation', min=0.01, max=0.99, value=0.2, step=0.01, width='100%')
						
						numericInput('k1', 'Pre-treatment k', min=0.1, max=10, value=1, step=0.1, width='100%'),
						numericInput('k2', 'Post-treatment k', min=0.1, max=10, value=1, step=0.1, width='100%'),
						numericInput('cor', 'Correlation', min=0.01, max=0.99, value=0.2, step=0.01, width='100%')
						
					)
				)
			),
			
			column(12,
				hr(),
				htmlOutput('footer'),
				hr()
			)
			
		),
		
		tabPanel("Probabilities",
			fluidRow(
				hr(),
				column(6,
					h4("Choose Statistical Method", style="text-align:center; "),
					selectInput("statmethod", "Statistical Method", c("BNB", "WAAVP", "Gamma", "Dobson", "Asymptotic"), selected="BNB", width='100%')
				),
				column(6,
					h4("Choose Statistical Framework", style="text-align:center; "),
					selectInput("statframe", "Statistical Framework", names(frameworks), selected="Denwood", width='100%')
				)
			),
			fluidRow(
				hr(),
				column(12,
					h4("Simulation Parameters"),
					htmlOutput('parameters'),
					hr(),
					h4("Probability of a non-inconclusive result"),
					plotOutput("noninc_plot"),
					hr(),
					h4("Probability of each classification type"),
					plotOutput("typologies_plot"),
					p("For an explanation of these classifications click on the 'Typologies' tab above.", style="text-align:center; "),
					hr()
				)
			)
		),
		
		tabPanel("Typologies",
			fluidRow(
				hr(),
				column(12,
					h4("To see typologies under a different framework, change the 'Statistical Framework' selection under the 'Probabilites' tab", style="text-align:center; "),
					hr(),
					imageOutput("typologies")
				)
			)
		)
	
	)
)

