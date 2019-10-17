# TODO:  add pooling and individuals per pool (pre and post?) questions and adjust reported k, and change prob_priors

initl <- "<br><h4>Error: Data inputs have not been initialised</h4>"
blankdf <- data.frame(Data=numeric(0))

function(input, output, session) {

	rv <- reactiveValues(predata = blankdf, postdata = blankdf, prebackup = blankdf, postbackup = blankdf, summaries="Select study design and enter data before calculating results!", initerrors="", dataerrors="", showreset=1, showresults=0, datainit=0, prelabel=initl, postlabel=initl, scalelabel="", prestr="control group / pre-treatment", poststr="treatment group / post-treatment", edt1=1, edt2=1)

	observeEvent(input$reset, {
		# Required to reset if dims don't change:
		rv$predata <- NULL
		rv$postdata <- NULL
		
		if(input$type == "Unpaired"){
			row1 <- input$Ncont
			row2 <- input$Ntx
			col1 <- input$Rcont
			col2 <- input$Rtx
			edt1 <- input$EDTcont
			edt2 <- input$EDTtx
		}else{
			row1 <- input$Npre
			row2 <- input$Npost
			col1 <- input$Rpre
			col2 <- input$Rpost
			edt1 <- input$EDTpre
			edt2 <- input$EDTpost
		}
		
		if(!parasitology){
			edt1 <- 1
			edt2 <- 1
		}
		
		prestr <- ifelse(input$type == "Unpaired", "control group", "pre-treatment")
		poststr <- ifelse(input$type == "Unpaired", "treatment group", "post-treatment")
		
		errors <- character(0)
		if(col1 != col2){
			errors <- c(errors, paste0("It is currently a requirement that the number of pre-treatment replicates equals the number of post-treatment replicates"))
		}
		if(row1 < 5 || round(row1)!=row1){
			errors <- c(errors, paste0("The ", prestr, " sample size must be a whole number >= 5"))
		}
		if(row2 < 5 || round(row2)!=row2){
			errors <- c(errors, paste0("The ", poststr, " sample size must be a whole number >= 5"))
		}
		if(col1 < 1 || round(col1)!=col1){
			errors <- c(errors, paste0("Zero, negative or non-integer ", prestr, " replicates"))
		}
		if(col2 < 1 || round(col2)!=col2){
			errors <- c(errors, paste0("Zero, negative or non-integer ", poststr, " replicates"))
		}
		if(edt1 <= 0){
			errors <- c(errors, paste0("Zero or negative ", prestr, " egg detection threshold"))
		}
		if(edt2 <= 0){
			errors <- c(errors, paste0("Zero or negative ", poststr, " egg detection threshold"))
		}
		if(length(errors)==0){
			rv$initerrors <- ""
		}else if(length(errors)==1){
			rv$initerrors <- paste0("<br>Error:  ", errors)
			return(1)
		}else{
			rv$initerrors <- paste0("<br>Errors:  ", paste(errors, collapse=", "))
			return(1)
		}
		rv$edt1 <- edt1
		rv$edt2 <- edt2
		# Don't save row and col here as it could be changed by the user
		
		if(testing){
			newdf <- lapply(1:col1, function(x) rnbinom(row1, 1, mu=15))
		}else{
			newdf <- lapply(1:col1, function(x) rep("", row1))
		}
		if(col1==1){
#			names(newdf) <- ifelse(input$type == "Unpaired", "Control", "PreTx")
			if(parasitology){
				names(newdf) <- ifelse(input$scale=="Raw Counts", "FEC", "EPG")
			}else{
				names(newdf) <- "Count"
			}
		}else{
#			names(newdf) <- paste0(ifelse(input$type == "Unpaired", "Control_Rep", "PreTx_Rep"), 1:col1)
			names(newdf) <- paste0("Rep_", 1:col1)
		}
		rv$predata <- as.data.frame(newdf)
		rv$prebackup <- rv$predata
		
		if(testing){
			tp <- sample(1:3, 1)
			if(tp==1) tvals <- 0
			if(tp==2) tvals <- 0:5
			if(tp==3) tvals <- 5:15				
			newdf <- lapply(1:col2, function(x) sample(tvals, row2, TRUE))
		}else{
			newdf <- lapply(1:col2, function(x) rep("", row2))
		}
		if(col2==1){
#			names(newdf) <- ifelse(input$type == "Unpaired", "Treatment", "PostTx")
			if(parasitology){
				names(newdf) <- ifelse(input$scale=="Raw Counts", "FEC", "EPG")
			}else{
				names(newdf) <- "Count"
			}
		}else{
#			names(newdf) <- paste0(ifelse(input$type == "Unpaired", "Treatment_Rep", "PostTx_Rep"), 1:col2)
			names(newdf) <- paste0("Rep_", 1:col2)
		}
		rv$postdata <- as.data.frame(newdf)
		rv$postbackup <- rv$postdata
		
		if(parasitology){
			scalelabel <- ifelse(input$scale=="Raw Counts", "(Enter data as raw egg counts", "(Enter data as eggs per gram")
		}else{
			scalelabel <- "(Enter count data"
		}
		if(col1 > 1 || col2 > 1){
			scalelabel <- paste0(scalelabel, ", with individuals in rows and replicates in columns)")
		}else{
			scalelabel <- paste0(scalelabel, ")")
		}
		rv$scalelabel <- scalelabel
		
		units <- "" #ifelse(edt1==1, "(raw counts)", "(eggs per gram)")
		rv$prelabel <- paste0(ifelse(input$type == "Unpaired", "<h4>Control Data ", "<h4>Pre-treatment Data "), units, "</h4>")
		units <- "" #ifelse(edt2==1, "(raw counts)", "(eggs per gram)")
		rv$postlabel <- paste0(ifelse(input$type == "Unpaired", "<h4>Treatment Data ", "<h4>Post-treatment Data "), units, "</h4>")
		rv$prestr <- prestr
		rv$poststr <- poststr

		rv$summaries <- ""
		
		# The reset buttton and nrow selectors can be hidden by setting to 0:
		if(!testing)
			rv$showreset <- 0
		
		rv$datainit <- 1
	})
	
	observeEvent(input$calculate, {
		
		rv$showresults <- 0
		rv$dataerrors <- ""
		
		if(rv$datainit==0){
			rv$dataerrors <- paste0("Error: The data inputs have not been initialised")
			return(1)
		}
		if(is.null(input$predata)){
			rv$dataerrors <- paste0("Error: The ", rv$prestr, " data has not been entered")
			return(1)
		}
		if(is.null(input$postdata)){
			rv$dataerrors <- paste0("Error: The ", rv$poststr, " data has not been entered")
			return(1)
		}
		
		te <- try(predata <- hot_to_r(input$predata))
		if(inherits(te,'try-error')){
			rv$dataerrors <- paste0("Error: Failed to read the ", rv$prestr, " data - this can happen when manually resizing the table - try entering the data again")
			rv$predata <- NULL
			rv$predata <- rv$prebackup
			return(1)
		}
		te <- try(postdata <- hot_to_r(input$postdata))
		if(inherits(te,'try-error')){
			rv$dataerrors <- paste0("Error: Failed to read the ", rv$poststr, " data - this can happen when manually resizing the table - try entering the data again")
			rv$postdata <- NULL
			rv$postdata <- rv$postbackup
			return(1)
		}
		
		if(nrow(predata)==0 || ncol(predata)==0){
			rv$dataerrors <- paste0("Error: Failed to initialise the ", rv$prestr, " data")
			return(1)
		}
		if(nrow(postdata)==0 || ncol(postdata)==0){
			rv$dataerrors <- paste0("Error: Failed to initialise the ", rv$poststr, " data")
			return(1)
		}
		
		if(any(is.na(predata)) || any(predata=="")){
			rv$dataerrors <- paste0("Error: Blank cells detected in the ", rv$prestr, " data")
			return(1)
		}
		if(any(is.na(postdata)) || any(postdata=="")){
			rv$dataerrors <- paste0("Error: Blank cells detected in the ", rv$poststr, " data")
			return(1)
		}
		
		predata <- as.matrix(as.data.frame(lapply(predata, as.numeric)))
		postdata <- as.matrix(as.data.frame(lapply(postdata, as.numeric)))

		if(any(is.na(predata))){
			rv$dataerrors <- paste0("Error: Non-numeric cells detected in the ", rv$prestr, " data")
			return(1)
		}
		if(any(is.na(postdata))){
			rv$dataerrors <- paste0("Error: Non-numeric cells detected in the ", rv$poststr, " data")
			return(1)
		}
		
		if(!parasitology || input$scale=="Raw Counts"){
			
			if(any(predata%%1 != 0)){
				rv$dataerrors <- paste0("Error: Non-integer cells detected in the ", rv$prestr, " data")
				return(1)
			}
			if(any(postdata%%1 != 0)){
				rv$dataerrors <- paste0("Error: Non-integer cells detected in the ", rv$poststr, " data")
				return(1)
			}
			
		}else{
			predata <- predata / rv$edt1
			postdata <- postdata / rv$edt2
		
			if(any(predata%%1 != 0)){
				rv$dataerrors <- paste0("Error: Non-integer cells detected in the ", rv$prestr, " data (after accounting for EDT)")
				return(1)
			}
			if(any(postdata%%1 != 0)){
				rv$dataerrors <- paste0("Error: Non-integer cells detected in the ", rv$poststr, " data (after accounting for EDT)")
				return(1)
			}
			
		}
		
		if(input$pthresh <= 0 || input$pthresh >= 1){
			rv$dataerrors <- paste0("Error:  The threshold for significance must be between 0-1")
			return(1)
		}
		if(! input$target < 100 ){
			rv$dataerrors <- paste0("Error:  The Target Effiacy value must be less than 100%")
			return(1)
		}
		
		# Now do analyses:
		
		premean <- mean(predata) * rv$edt1
		postmean <- mean(postdata) * rv$edt2
		
		Npre <- nrow(predata)
		Rpre <- ncol(predata)
		Npost <- nrow(postdata)
		Rpost <- ncol(postdata)

		paired <- input$type=="Paired"
		tail <- input$pthresh
		
		# TODO: make sure mean_rato works OK
		mean_ratio <- 1
		
		results <- efficacy_analysis(predata, postdata, input$target/100, (input$target-input$nim)/100, paired, tail)
		
		obsred <- round( (1 - (mean(postdata)/mean(predata))) * 100 , 1)
		
		if(parasitology){
			outstring <- paste0("<strong>Summary statistics:</strong><br> &nbsp; The ", rv$prestr, " mean is ", round(premean, 1), "EPG (sample size = ", Npre, if(Rpre>1) paste0("x", Rpre), ")<br> &nbsp; The ", rv$poststr, " mean is ", round(postmean, 1), "EPG (sample size = ", Npost, if(Rpost>1) paste0("x", Rpost), ")<br>")
		}else{
			outstring <- paste0("<strong>Summary statistics:</strong><br> &nbsp; The ", rv$prestr, " mean is ", round(premean, 1), " (sample size = ", Npre, if(Rpre>1) paste0("x", Rpre), ")<br> &nbsp; The ", rv$poststr, " mean is ", round(postmean, 1), " (sample size = ", Npost, if(Rpost>1) paste0("x", Rpost), ")<br>")
		}
		outstring <- paste0(outstring, " &nbsp; The threshold for inferiority (T_I) is: &nbsp; ", input$target, "%<br> &nbsp; The threshold for non-inferiority (T_A) is: &nbsp; ", input$target-input$nim, "%<br><br>")
		
		colouredclass <- function(x) paste0("<span style='color:", switch(x, "Reduced"="red", "Inconclusive"="grey", "Marginal"="orange", "Adequate"="blue", "black"), ";'>", x, if(x!="Inconclusive" && !grepl("Error",x)) " Efficacy", "</span>")
		
		pci <- function(x, y, adj=1) paste0(100 * (1-(input$pthresh*adj*2)), "% CI: &nbsp; ", round(x*100,1), "-", round(y*100,1), "%")
		pp <- function(x) if(x < 0.001) "<0.001" else format(c(round(x,3),0.001))[1]
		
		resrow <- results[results$Method=="BNB",]
		outstring <- paste0(outstring, "<strong>BNB method:</strong> &nbsp; &nbsp; p_I: ", pp(resrow$pI), "; &nbsp; &nbsp; p_A:  ", pp(resrow$pA), "<br>&nbsp; &nbsp; ", colouredclass(resrow$Classification), "<br><br>")

		resrow <- results[results$Method=="WAAVP",]
		outstring <- paste0(outstring, "<strong>WAAVP method:</strong> &nbsp; &nbsp; ", pci(resrow$LCI, resrow$UCI), "<br>&nbsp; &nbsp; ", colouredclass(resrow$Classification), "<br><br>")

		resrow <- results[results$Method=="Gamma",]
		outstring <- paste0(outstring, "<strong>Gamma method:</strong> &nbsp; &nbsp; ", pci(resrow$LCI, resrow$UCI), "<br>&nbsp; &nbsp; ", colouredclass(resrow$Classification), "<br><br>")

		resrow <- results[results$Method=="Asymptotic",]
		outstring <- paste0(outstring, "<strong>Asymptotic method:</strong> &nbsp; &nbsp; ", pci(resrow$LCI, resrow$UCI), "<br>&nbsp; &nbsp; ", colouredclass(resrow$Classification), "<br><br>")
		
		if(paired){
			resrow <- results[results$Method=="Binomial",]
			outstring <- paste0(outstring, "<strong>Binomial method:</strong> &nbsp; &nbsp; ", pci(resrow$LCI, resrow$UCI, 0.2), "<br>&nbsp; &nbsp; ", colouredclass(resrow$Classification), "<br><br>")
		}
		
		if(FALSE){
		
			# Track if we get a bad pvalue:
			anerr <- FALSE
			res <- list(p_1 = results$pA[results$Method=="BNB"], p_2 = results$pI[results$Method=="BNB"])
		
		
		
			# If the observed reduction is above the target then don't report the inf test:
			if(obsred >= input$target){
				inf <- FALSE
				if(obsred == input$target){
					outstring <- paste0(outstring, "<strong>Inferiority test:</strong> The observed FECR of ", obsred, "% is equal to the specified target efficacy<br>")
				}else{
					outstring <- paste0(outstring, "<strong>Inferiority test:</strong> The observed FECR of ", obsred, "% is greater than the specified target efficacy<br>")
				}
			}else{
				pval <- res$p_2
				if(is.na(pval) || pval==Inf || pval < -0.001 || pval > 1){
					anerr <- TRUE
					outstring <- paste0(outstring, "<strong>Inferiority test:</strong> The non-inferiority test result could not be calculated<br>")
				}else{
					pval <- round(pval, 3)			
					pstr <- ifelse(pval < 0.001, "<0.001", paste0("=",pval))
					if(pval <= input$pthresh){
						outstring <- paste0(outstring, "<strong>Inferiority test:</strong> The observed FECR of ", obsred, "% is <span style='color:red;'>significantly inferior</span> to the target of ", input$target, "% (p", pstr, ")<br>")
						inf <- TRUE
					}else{
						outstring <- paste0(outstring, "<strong>Inferiority test:</strong> The observed FECR of ", obsred, "% is not significantly inferior to the target of ", input$target, "% (p", pstr, ")<br>")
						inf <- FALSE
					}
				}
			}
		
			# If the observed reduction is below the margin then don't report the non-inf test:
			if(obsred < (input$target-input$nim)){
				ninf <- FALSE
				outstring <- paste0(outstring, "<strong>Non-inferiority test:</strong> The observed FECR of ", obsred, "% is below the specified non-inferiority margin of the target efficacy<br>")
			}else{
				pval <- res$p_1
				if(is.na(pval) || pval==Inf || pval < -0.001 || pval > 1){
					anerr <- TRUE
					outstring <- paste0(outstring, "<strong>Non-inferiority test:</strong> The non-inferiority test result could not be calculated<br>")
				}else{
					pval <- round(pval, 3)			
					pstr <- ifelse(pval < 0.001, "<0.001", paste0("=",pval))
					if(pval <= input$pthresh){
						outstring <- paste0(outstring, "<strong>Non-inferiority test:</strong> The observed FECR of ", obsred, "% is <span style='color:blue;'>significantly non-inferior</span> to the target of ", input$target, "% with given margin (p", pstr, ")<br>")
						ninf <- TRUE
					}else{
						outstring <- paste0(outstring, "<strong>Non-inferiority test:</strong> The observed FECR of ", obsred, "% is not significantly non-inferior to the target of ", input$target, "% with given margin (p", pstr, ")<br>")
						ninf <- FALSE
					}
				}
			}
		
			if(anerr){
				class <- "The classification could not be determined"
			}else{
				if(inf && !ninf){
					class <- "<span style='color:red;'>Reduced Efficacy</span>"
				}else if(!inf && !ninf){
					class <- "<span style='color:grey;'>Inconclusive</span>"
				}else if(inf && ninf){
					class <- "<span style='color:orange;'>Marginal Efficacy</span>"
				}else if(!inf && ninf){
					class <- "<span style='color:blue;'>Adequate Efficacy</span>"
				}else{
					class <- "ERROR DETERMINING CLASS"
				}
			}
			outstring <- paste0(outstring, "<br><strong>Classification:</strong> ", class)
		
		}
		
		
		
		rv$showresults <- 1
		rv$summaries <- outstring
	})
	
	fluidPage(
		output$predata <- renderRHandsontable({
			rhandsontable(rv$predata, rowHeaders=NULL, useTypes = FALSE, stretchH = "none")
		}),
		output$postdata <- renderRHandsontable({
			rhandsontable(rv$postdata, rowHeaders=NULL, useTypes = FALSE, stretchH = "none")
		}),
		output$summaries <- renderText({
			rv$summaries
		}),
		output$showreset <- renderText({
			rv$showreset
		}),
		output$showresults <- renderText({
			rv$showresults
		}),
		output$initerrors <- renderText({
			rv$initerrors
		}),
		output$dataerrors <- renderText({
			paste0("<br>", rv$dataerrors)
		}),
		output$prelabel <- renderText({
			rv$prelabel
		}),
		output$postlabel <- renderText({
			rv$postlabel
		}),
		output$scalelabel <- renderText({
			rv$scalelabel
		})
	)
	
    output$footer <- renderText(paste0("<p align='center'>Data analysis tool for ERR / FECRT data by Matthew Denwood<br>", footeraddtext, "<a href='http://www.fecrt.com/framework/', target='_blank'>Click here for more information (opens in a new window)</a></p>"))
	

	outputOptions(output, "showreset", suspendWhenHidden=FALSE)
	outputOptions(output, "showresults", suspendWhenHidden=FALSE)
	
	# This breaks stuff:
#	outputOptions(output, "predata", suspendWhenHidden=FALSE)
#	outputOptions(output, "postdata", suspendWhenHidden=FALSE)

}
