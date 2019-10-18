# Nav list for tabs:  https://shiny.rstudio.com/articles/layout-guide.html

fluidPage(

    tags$head(headscript),

	theme = shinytheme("cerulean"),

	hr(),
#	h3("Select Parameter Values", style="text-align:center"),
#	hr(),
	
	tabsetPanel(
		tabPanel("Parameters",
			conditionalPanel(
				condition = "output.showreset == 1",
				hr(),
				h4("Select Study Design Parameters", style="text-align:center; "),
				hr(),
				if(parasitology){
				fluidRow(
					column(6,				
						selectInput("type", "Study Type", c("Unpaired","Paired"), selected="Unpaired", width='100%')
					),
					column(6,
						selectInput("scale", "Enter Data as", c("Raw Counts", "Eggs Per Gram"), selected="Raw Counts", width='100%')
					)
				)
				}else{
					fluidRow(
						column(12,				
							selectInput("type", "Study Type", c("Unpaired","Paired"), selected="Unpaired", width='100%')
						)
					)
				}
			),
			
			conditionalPanel(
				condition = "output.showreset == 1 && input.type == 'Paired'",
				
				fluidRow(
					column(6,
						numericInput('Npre', 'Pre-treatment sample size', min=1, value=10, step=1, width='100%'),
						numericInput('Rpre', 'Pre-treatment replicates', min=1, value=1, step=1, width='100%'),
						if(parasitology) numericInput('EDTpre', 'Pre-treatment group Egg Detection Threshold', value=25, min=0, width='100%')
					),
					column(6,
						numericInput('Npost', 'Post-treatment sample size', min=1, value=10, step=1, width='100%'),
						numericInput('Rpost', 'Post-treatment replicates', min=1, value=1, step=1, width='100%'),
						if(parasitology) numericInput('EDTpost', 'Post-treatment group Egg Detection Threshold', value=25, min=0, width='100%')
					)
				)
			),
			conditionalPanel(
				condition = "output.showreset == 1 && input.type == 'Unpaired'",
				
				fluidRow(
					column(6,
						numericInput('Ncont', 'Control group sample size', min=1, value=10, step=1, width='100%'),
						numericInput('Rcont', 'Control group replicates', min=1, value=1, step=1, width='100%'),
						if(parasitology) numericInput('EDTcont', 'Control group Egg Detection Threshold', value=25, min=0, width='100%')
					),
					column(6,
						numericInput('Ntx', 'Treatment group sample size', min=1, value=10, step=1, width='100%'),
						numericInput('Rtx', 'Treatment group replicates', min=1, value=1, step=1, width='100%'),
						if(parasitology) numericInput('EDTtx', 'Treatment group Egg Detection Threshold', value=25, min=0, width='100%')
					)
				)
			),
			
			conditionalPanel(
				condition = "output.showreset == 1",
				hr(),
				h4("Initialise the Data Inputs", style="text-align:center; "),			
				hr(),
				fluidRow(
					column(6, align="center", offset = 3,
						actionButton("reset", "Click to Initialise Data Inputs", icon = icon("paper-plane"), width='100%')
					)
				),
				br(),
				htmlOutput("initerrors", style="text-align:center; color:red; ", width="100%")
			),
				
			conditionalPanel(
				condition = "output.showreset == 0",
				br(),
				hr(),
				h4("Now enter your data into the 'Data Input' tab above", style="text-align:center; ")
			)	
		),
		tabPanel("Data Input",
			hr(),
			h4("Copy and paste your data from Excel, then click on the 'Results' tab above", style="text-align:center; "),			
			htmlOutput('scalelabel', style="text-align:center; ", width="100%"),
			hr(),
			fluidRow( 
				column(6,
					htmlOutput('prelabel', style="text-align:left; ", width="100%"),
					rHandsontableOutput("predata")
				),
				column(6,
					htmlOutput('postlabel', style="text-align:left; ", width="100%"),
					rHandsontableOutput("postdata")
				)
			)
		),
		tabPanel("Results", 
			hr(),
			h4("Select Test Parameters", style="text-align:center; "),			
			hr(),
	        sliderInput('target', 'Target efficay (%)', min=0, max=100, value=99, step=0.1, width='100%'),
	        sliderInput('nim', 'Non-inferiority margin (% points)', min=0, max=25, value=4, step=0.1, width='100%'),
			numericInput('pthresh', 'Threshold for significance (p)', min=0, max=0.2, value=0.025, step=0.025, width='100%'),
			hr(),
			fluidRow(
				column(6, align="center", offset = 3,
					actionButton("calculate", "Click to Calculate", icon = icon("paper-plane"), width='100%')
				)
			),
			htmlOutput("dataerrors", style="text-align:center; color:red; ", width="100%"),
			conditionalPanel(
				condition = "output.showresults == 1",
				hr(),
				h4("Results", style="text-align:center; "),			
				hr(),
				htmlOutput('summaries', width="100%")
			)
		)
	),

	hr(),
	htmlOutput('footer'),
	hr()
)
