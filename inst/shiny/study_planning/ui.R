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

			h5("This tool is currently down for maintenance - please check back on 2023-08-21")

		),
		## /PARAMETERS


		## CALCULATION
		tabPanel("Calculation",

		h5("This tool is currently down for maintenance - please check back on 2023-08-21")

		),
		## /CALCULATION


		## RESULTS
		tabPanel("Results",
		h5("This tool is currently down for maintenance - please check back on 2023-08-21")
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
