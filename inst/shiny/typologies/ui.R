fluidPage(

    tags$head(headscript),

	theme = shinytheme("cerulean"),

	hr(),
	h4("Select parameter values:", style="text-align:left; "),
	
	fluidRow(
		column(12,
			checkboxInput('paired', 'Paired study?', value=FALSE)
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
		column(6,
			numericInput('Sum_1', 'Control/Pre-treatment sum count', min=1, value=200, step=1, width='100%'),
			numericInput('k_1', 'Control/Pre-treatment k', min=0.01, value=1, step=0.01, width='100%'),
			numericInput('T_I', 'Threshold for inferioirty', min=0, max=1, value=0.99, step=0.001, width='100%')
		),
		column(6,
			numericInput('Sum_2', 'Treatment/Post-treatment sum count', min=0, value=5, step=1, width='100%'),
			numericInput('k_2', 'Treatment/Post-treatment k', min=0.01, value=1, step=0.01, width='100%'),
			numericInput('T_A', 'Threshold for non-inferioirty', min=0, max=1, value=0.95, step=0.001, width='100%')
		)
	),
	
	hr(),
	htmlOutput('efficacy'),
	tableOutput('tbl'),
	
	hr(),
	htmlOutput('footer'),
	hr()
)
