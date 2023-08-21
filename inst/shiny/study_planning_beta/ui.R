fluidPage(
	theme = shinytheme("readable"),

	tags$head(includeHTML("head.html")),
#	tags$header(hr()),

	navbarPage("fecrt.com",

		## INSTRUCTIONS
		tabPanel("Instructions",
			includeMarkdown("instructions.md")
		),
		## /INSTRUCTIONS


		## PARAMETERS
		tabPanel("Parameters",

			fluidRow(
				column(colwidth,
					selectInput("species",
						"Select host/parasite species:",
						choices = species_choices,
						selected = character(0),
						selectize = FALSE
					)
				),
				column(colwidth,
					numericInput("targetEfficacy",
						"Select target efficacy (%):",
						value = 99,
						min = -100,
						max = 100,
						step = NA
					)
				)
			),

			fluidRow(
				column(colwidth,
					numericInput("delta",
						"Select non-inferiority margin (%):",
						value = 4,
						min = 1,
						max = 100,
						step = NA
					)
				),
				column(colwidth,
					htmlOutput("lowerThreshold")
				)
			),

			fluidRow(
				column(colwidth,
					numericInput("meanEPG",
						"Select expected pre-treatment mean EPG:",
						value = 1000,
						min = 0,
						max = NA,
						step = NA
					)
				),
				column(colwidth,
					numericInput("multfact",
						"Select multiplication factor:",
						value = 50,
						min = 0,
						max = NA,
						step = NA
					)
				),

				fluidRow(
					column(colwidth*2,
						numericInput("maxN",
							"Select maximum sample size:",
							value = 100,
							min = 10,
							max = NA,
							step = NA
						)
					)
				)
			)
		),
		## /PARAMETERS


		## CALCULATION
		tabPanel("Calculation",

			fluidRow(
				column(colwidth*2,
					actionButton("calculate",
						"Click to Calculate",
						icon = NULL
					)
				)
			),

			conditionalPanel("output.status == 0",
				fluidRow(
					column(colwidth*2,
						htmlOutput("parameter_feedback", width="100%")
					)
				)
			),

			conditionalPanel("output.status == 1",
				fluidRow(
					column(colwidth*2,
						htmlOutput("results_text", width="100%")
					)
				),
				fluidRow(
					column(colwidth*2,
						plotOutput("results_plot", width="100%")
					)
				)
			),

		)
	),

# 	tags$footer(HTML("
#                     <!-- Footer -->
#                            <footer class='page-footer font-large indigo'>
#                            <!-- Copyright -->
#                            <div class='footer-copyright text-center py-3'>© 2018 Copyright:
#                            <a href='https://mdbootstrap.com/education/bootstrap/'> MDBootstrap.com</a>
#                            </div>
#                            <!-- Copyright -->
#
#                            </footer>
#                            <!-- Footer -->")
# 		)

	tags$footer(hr())
)
