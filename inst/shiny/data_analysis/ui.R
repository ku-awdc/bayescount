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


		## DATA
		tabPanel("Data",

			## Top row:
			fluidRow(
				column(colwidth,
					radioButtons("entryType",
						"Select data entry method:",
						choices = c(`File upload` = "file", `Direct entry` = "direct", `Demonstration` = "demo"),
						selected = character(0),
						inline = FALSE
					)
				),
				conditionalPanel(
					"input.entryType == 'direct'",
					column(colwidth,

						## Select input for study design:
						selectInput("design", "Study design:",
							choices = c(`Paired` = "paired", `Treatment / Control` = "unpaired"),
							selected = "paired",
							selectize = FALSE
						)
					)
				),
				conditionalPanel(
					"input.entryType == 'file'",
					column(colwidth,
						fileInput("dataFile", "Choose data file(s):",
							multiple = TRUE,
							accept = c(".csv", ".xls", ".xlsx")
						)
					)
				)
			),

			## Conditional row for feedback on file upload (needed for reactive):
			conditionalPanel(
				"input.entryType == 'file'",
				fluidRow(
					column(colwidth*2,
						htmlOutput("upload_feedback", width="100%")
					)
				)
			),

			## Conditional row for direct entry paired N:
			conditionalPanel(
				"input.entryType == 'direct' && input.design == 'paired'",

				fluidRow(
					column(colwidth*2,
						numericInput("directN_paired",
							"Enter number of animals:",
							value = 20,
							min = 1,
							step = 1
						)
					)
				)
			),

			## Conditional row for direct entry unpaired N:
			conditionalPanel(
				"input.entryType == 'direct' && input.design == 'unpaired'",

				# Non-equal N currently disabled as broken for BNB:
				# fluidRow(
				# 	column(colwidth,
				# 		numericInput("directN_ctl",
				# 			"Enter control animals:",
				# 			value = 20,
				# 			min = 1,
				# 			step = 1
				# 		)
				# 	),
				# 	column(colwidth,
				# 		numericInput("directN_txt",
				# 			"Enter treatment animals:",
				# 			value = 20,
				# 			min = 1,
				# 			step = 1
				# 		)
				# 	)
				# )
				fluidRow(
					column(colwidth*2,
						numericInput("directN_ctl",
							"Enter number of animals:",
							value = 20,
							min = 1,
							step = 1
						)
					)
				)
			),

			## Conditional row for initialise:
			conditionalPanel(
				"input.entryType == 'direct' && output.status == 0",
				fluidRow(
					column(colwidth*2,
						actionButton("initialise_data",
							"Initialise data entry"
						)
					)
				)
			),

			## Conditional row for data paired:
			conditionalPanel(
				"input.entryType == 'direct' && input.design == 'paired' && output.status >= 1",
				fluidRow(
					column(colwidth*2,
						rHandsontableOutput("data_paired")
					)
				)
			),

			## Conditional row for data unpaired:
			conditionalPanel(
				"input.entryType == 'direct' && input.design == 'unpaired' && output.status >= 1",
				fluidRow(
					column(colwidth,
						rHandsontableOutput("data_ctl")
					),
					column(colwidth,
						rHandsontableOutput("data_txt")
					)
				)
			),

			## Conditional row for data demo:
			conditionalPanel(
				"input.entryType == 'demo'",
				fluidRow(
					column(colwidth*2,
						h5("Demonstration data (paired study; N=20):")
					)
				),
				fluidRow(
					column(colwidth*2,
						rHandsontableOutput("data_demo")
					)
				)
			)
		),
		## /DATA


		## PARAMETERS
		tabPanel("Parameters",

			## Row for pre-selection and host:
			fluidRow(
				column(colwidth,
					radioButtons("parameterType",
						"Select scenario:",
						#choices = c(`Guidelines: clinical` = "clinical", `Guidelines: research` = "research", `Custom` = "custom"),
						choices = c(`Guidelines: clinical` = "clinical", `Guidelines: research` = "research"),
						selected = "waavp", #character(0),
						inline = FALSE
					)
				),

				conditionalPanel("input.parameterType == 'clinical' || input.parameterType == 'research'",
					column(colwidth,
						selectInput("waavpSetup", "Species and anthelmintic:",
							choices = waavp_choices,
							selected = NULL,
							selectize = FALSE
						)
					)
				),
				conditionalPanel("input.parameterType == 'custom'",
					column(colwidth,
						selectInput("hostSpecies", "Host species:",
							choices = host_choices,
							selected = NULL,
							selectize = FALSE
						)
					)
				)
			),

			## Optional rows for custom parameters:
			conditionalPanel("input.parameterType == 'custom'",
				fluidRow(
					## TODO


					## Custom settings for custom option:
					conditionalPanel("input.parameterType == 'custom'",

						## Select input for host (always displayed):
						selectInput("hostSpecies", "Host species:",
							choices = list(
								`*SELECT*` = "INVALID",
								Ruminants = list(`Cattle` = "cattle", `Sheep` = "sheep", `Goats` = "goats"),
								Equine = list(`Horses (adults)` = "horse_adult", `Horses (foals)` = "Horses (foals)", `Donkeys (adults)` = "donk_adult", `Donkeys (foals)` = "donk_foal"),
								`Swine` = "pigs",
								`Other` = "other"
							),
							selected = NULL,
							selectize = FALSE
						),
						## Text input for host (if hostSpecies->other)
						conditionalPanel(
							condition = "input.hostSpecies == 'other'",
							textInput("host_other", "Enter host species:")
						),

						## Select input for parasite (if ! hostSpecies->other)
						conditionalPanel("input.hostSpecies != 'other'",
							## Note: choices are reactive
							selectInput("parasiteSpecies", "Parasite species:",
								choices = NULL,
								selectize = FALSE
							)
						),
						## Text input for parasite (if hostSpecies->other or parasiteSpecies->other)
						conditionalPanel("input.hostSpecies == 'other' || input.parasiteSpecies == 'other'",
							textInput("parasite_other", "Enter parasite species:")
						),

						## Expected variabilities x3 - NOTE: reactive
						numericInput("k1", "", "", min=0, max=20, step=0.1),
						numericInput("k2", "", "", min=0, max=20, step=0.1),
						numericInput("kc", "", "", min=0, max=20, step=0.1),


						## Select input for anthelmintic - note: choices are reactive
						selectInput("anthelmintic", "Anthelmintic class:",
							choices = NULL,
							selectize = FALSE
						),
						## Text input for anthelmintic (if anthelmintic->other)
						conditionalPanel("input.anthelmintic == 'other'",
							textInput("anthelmintic_other", "Enter anthelmintic class:")
						),

						## Text input for parasite (if hostSpecies->other or parasiteSpecies->other)
						conditionalPanel("input.hostSpecies == 'other' || input.parasiteSpecies == 'other'",
							## Note: choices are reactive
							textInput("parasite_other", "Enter parasite species:")
						),


						## Target, Lower

						## Button input for pre-selection:
						## alpha level
						numericInput("alphaLevel",
							"Desired type 1 error:",
							value = 0.05,
							min = 0,
							max = 0.5,
							step = 0.005
						)

					)

				)
			),

			## Common rows for multiplication factor:
			conditionalPanel("input.parameterType == 'clinical' || input.parameterType == 'research' || input.parameterType == 'custom'",
				conditionalPanel("input.design == 'paired'",
					fluidRow(
						column(colwidth,
							div(strong("Multiplication factor (pre-treatment):")),
							numericInput("mf_pre",
								"",
								NULL,
								min=0, max=100, step=1
							)
						),
						column(colwidth,
							checkboxInput("mfp_fixed",
								"Post- same as pre-treatment",
								value=TRUE
							),
							conditionalPanel("input.mfp_fixed == false",
								numericInput("mf_post",
									"",
									NULL,
									min=0, max=100, step=1
								)
							)
						)
					)
				),

				conditionalPanel("input.design == 'unpaired'",
					fluidRow(
						column(colwidth,
							div(strong("Multiplication factor (controls):")),
							numericInput("mf_ctl",
								"",
								NULL,
								min=0, max=100, step=1
							)
						),
						column(colwidth,
							checkboxInput("mfu_fixed",
								"Treatment same as control",
								value=TRUE
							),
							conditionalPanel("input.mfu_fixed == false",
								numericInput("mf_txt",
									"",
									NULL,
									min=0, max=100, step=1
								)
							)
						)
					)
				),

				## Row for optional info:
				fluidRow(
					column(colwidth,
						textInput("region",
							"Country/region (optional):",
							value = ""
						),
					),
					column(colwidth,
						textInput("identifier",
							"Study identifier (optional):",
							value = ""
						)
					)
				)
			)
		),
		## /PARAMETERS


		## RESULTS
		tabPanel("Results",

			## First some feedback text to say what is missing:
			conditionalPanel("output.status < 3",
				fluidRow(
					column(colwidth*2,
						htmlOutput('status_feedback', width="100%")
					)
				)
			),

			## Then a button to calculate:
			conditionalPanel("output.status == 3",
				fluidRow(
					column(colwidth*2,
						actionButton("calculate",
							"Click to Calculate",
							icon = NULL
						)
					)
				)
			),

			## Then results and download tools:
			conditionalPanel("output.status == 4",
				fluidRow(
					column(colwidth*2,
						htmlOutput("result_summary", width="100%"),
						hr(),
						selectInput("downloadType",
							"Report download format:",
							choices = c(`PDF` = "pdf", `Microsoft Word` = "word"),
							selected="pdf"
						),
						downloadButton("downloadReport", "Download report")
					)
				)
			),

			## TODO: icon - https://fontawesome.com/icons or https://getbootstrap.com/docs/3.3/components/#glyphicons

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
