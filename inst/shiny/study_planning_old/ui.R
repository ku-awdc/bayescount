fluidPage(

    tags$head(headscript),

	theme = shinytheme("cerulean"),

	hr(),
	h4("Select parameter values", style="text-align:center; "),
	
	fluidRow(
		column(12,
			checkboxInput('paired', 'Paired study?', value=TRUE)
		)
	),
		
	conditionalPanel(
		condition = "input.paired == false",
		fluidRow(
			column(6,
				numericInput('N_1', 'Control sample size', min=1, value=10, step=1, width='100%')
			),
			column(6,
				numericInput('N_2', 'Ttreatment sample size', min=1, value=10, step=1, width='100%')
			)
		)
	),

	conditionalPanel(
		condition = "input.paired == true",
		fluidRow(
			column(6,
				numericInput('N', 'Sample size', min=1, value=10, step=1, width='100%')
			),
			column(6,
				numericInput('cor', 'Correlation', min=0, max=1, value=0.25, step=0.01, width='100%')
			)
		)
	),
	
	fluidRow(
		column(12,
			numericInput('mean', 'Control/Pre-treatment mean count', min=1, value=50, step=0.1, width='100%')
		)
	),
	
	fluidRow(
		column(6,
			numericInput('k_1', 'Control/Pre-treatment k', min=0.01, value=1, step=0.01, width='100%'),
			numericInput('T_I', 'Threshold for inferioirty', min=0, max=1, value=0.99, step=0.001, width='100%')
		),
		column(6,
			numericInput('k_2', 'Treatment/Post-treatment k', min=0.01, value=0.75, step=0.01, width='100%'),
			numericInput('T_A', 'Threshold for non-inferioirty', min=0, max=1, value=0.95, step=0.001, width='100%')
		)
	),

	hr(),
	fluidRow(
		column(12,
			selectInput('method', 'Select analysis method', choices=c("BNB", "WAAVP", "Gamma", "Asymptotic", "Binomial"), selected="BNB", width='100%')
		)
	),
	hr(),
	h4("Run the Monte Carlo simulation", style="text-align:center; "),
	fluidRow(
		column(6, align="center", offset = 3,
			actionButton("simulate", "Click to run the simulation with these parameters [this takes up to a few seconds to process]", icon = icon("paper-plane"), width='100%')
		)
	),
	hr(),
		
	h4("Results", style="text-align:center; "),
	hr(),
	
	htmlOutput('results'),
	fluidRow(
		column(6,
			htmlOutput('tta'),
			tableOutput('tbl_a'),
			htmlOutput('bta')
		),
		column(6,
			htmlOutput('tti'),
			tableOutput('tbl_i'),
			htmlOutput('bti')
		)
	),
	
	hr(),
	fluidRow(
		column(12,
			htmlOutput('interp')
		)
	),
	
	hr(),
	htmlOutput('footer'),
	hr()
)
