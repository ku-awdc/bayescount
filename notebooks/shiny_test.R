library(shiny)

ui <- fluidPage(

	actionButton("calculate",
		"Click to Calculate",
		icon = NULL
	),

	conditionalPanel("output.ready == 'yes'",
		selectInput("downloadType",
			"Select report file type for download:",
			choices = c(`PDF` = "pdf", `Microsoft Word` = "word"),
			selected="pdf"
		)
	)
)

server <- function(input, output) {

	rv <- reactiveValues(ready="no")
	observeEvent(input$calculate, {
		rv$ready <- "yes"
	})

	output$ready <- renderText(rv$ready)

	outputOptions(output, 'ready', suspendWhenHidden=FALSE)
}

shinyApp(ui, server)
