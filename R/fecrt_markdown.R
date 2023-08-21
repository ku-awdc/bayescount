## Internal function taking results and producing a character vector of markdown file
pvround <- function(x){
	x <- str_c("p = ", format(c(round(x, 3),0.001))[1:length(x)])
	x[x=="p = 0.000"] <- "p < 0.001"
	x
}

fecrt_markdown <- function(input, parameters, design){

	stopifnot(length(design)==1, design %in% c("paired","unpaired"))
	results <- input$full

	if(is.na(parameters$identifier)){
		title <- "FECRT analysis"
	}else{
		title <- str_c("Analysis of FECRT data from ", parameters$identifier)
	}

	author <- "fecrt.com data analysis tool"

	output <- str_c('---
title: "', title, '"
author: "', author, '"
date: "', as.character(Sys.Date()), '"
---
')

	if(design=="paired"){
		text1 <- "pre-treatment"
		text2 <- "post-treatment"
		is_paired <- TRUE
	}else if(design=="unpaired"){
		text1 <- "control"
		text2 <- "treatment"
		is_paired <- FALSE
	}else{
		stop("Unrecognised design")
	}

	stopifnot(input$n > 0)

	output <- str_c(output, '

## Efficacy classification', if(input$n>1) 's', '

', str_replace_all(input$headline, "<br>", "\n\n"), '
')

	for(i in seq_len(input$n)){

		if(input$n==1){
			output <- str_c(output, '

## Detailed results

')
		}else{
			output <- str_c(output, '

## Detailed results for ', results[[i]]$name, '

')
		}

		output <- str_c(output, '

### Summary statistics

', if(is_paired){
		str_c('Number of animals: ', results[[i]]$sumstats$n1)
	}else{
		str_c('Number of ', text1, ' animals: ', results[[i]]$sumstats$n1, '
Number of ', text2, ' animals: ', results[[i]]$sumstats$n2)
	}, '\n
Mean of ', text1, ' data:  ', results[[i]]$sumstats$m1, '\n
Mean of ', text2, ' data:  ', results[[i]]$sumstats$m2, '\n
Variance of ', text1, ' data:  ', results[[i]]$sumstats$v1, '\n
Variance of ', text2, ' data:  ', results[[i]]$sumstats$v2, '\n
Estimated over-dispersion (k) of ', text1, ' data:  ', results[[i]]$sumstats$ks[1], '\n
Estimated over-dispersion (k) of ', text2, ' data:  ', results[[i]]$sumstats$ks[2], '\n
', if(is_paired) str_c('Estimated within-animal correlation:  ', results[[i]]$sumstats$ks[3]), '\n', '

')

		allres <- results[[i]]$all |>
			mutate(CIstring = case_when(
				!is.na(LCI) & !is.na(UCI) ~ str_c('90% CI = ', round(LCI*100,1), '% - ', round(UCI*100,1), '%'),
				TRUE ~ "[90% CI uncalculable]"
			)) |>
			mutate(PVstring = case_when(
				!is.na(pI) & !is.na(pA) ~ str_c('Test for Resistance: ', pvround(pI), ';  Test for Susceptibility: ', pvround(pA)),
				TRUE ~ "[test p-values uncalculable]"
			))

		output <- str_c(output, '

### Results from the Delta method (Levecke et al.)

Classification:  ', allres |> filter(Method=="Gamma") |> pull(Classification), '\n
', allres |> filter(Method=="Gamma") |> pull(CIstring), '\n

Notes:

- This method is non-parametric, so is robust to distributional assumptions, and variances of the ', text1, ' and ', text2, ' data are estimated independently
- This method cannot be used when the ', text2, ' data are all zero
- This method may give misleading results with fewer than five observations, and when fewer than three ', text2, ' observations are non-zero, due to unstable variance estimates
- This is the preferred method when the sample size is greater than or equal to 5, and where at least three ', text2, ' observations are non-zero


### Results from the WAAVP method (Coles et al. and Pepper et al.)

Classification:  ', allres |> filter(Method=="WAAVP") |> pull(Classification), '\n
', allres |> filter(Method=="WAAVP") |> pull(CIstring), '\n

Notes:

- This method is non-parametric, so is robust to distributional assumptions, and variances of the ', text1, ' and ', text2, ' data are estimated independently
- This method cannot be used when the ', text2, ' data are all zero
- This method may give misleading results with fewer than five observations, and when fewer than three ', text2, ' observations are non-zero, due to unstable variance estimates


### Results from the BNB method (Denwood et al.) version A

Classification:  ', allres |> filter(Method=="BNB") |> pull(Classification), '\n
', allres |> filter(Method=="BNB") |> pull(PVstring), '\n

Notes:

- This method is parametric, and assumes that the data follow a negative binomial distribution:  the classification will be unavailable if the multiplication factor you entered does not match the data
- The over-dispersion is estimated independently for the ', text1, ' and ', text2, 'data
- This method may give misleading results with fewer than five observations, and when fewer than three ', text2, ' observations are non-zero, due to unstable estimates of over-dispersion


### Results from the BNB method (Denwood et al.) version B

Classification:  ', allres |> filter(Method=="BNB_FixK2") |> pull(Classification), '\n
', allres |> filter(Method=="BNB_FixK2") |> pull(PVstring), '\n

Notes:

- This method is parametric, and assumes that the data follow a negative binomial distribution:  the classification will be unavailable if the multiplication factor you entered does not match the data
- The over-dispersion is estimated for the ', text1, ' data, but the over-dispersion in the ', text2, 'data is assumed to be proportional to that of the ', text1, 'data
- This method may give misleading results with fewer than five observations due to unstable estimates of over-dispersion
- This is the preferred method when the sample size is greater than or equal to 5, and where fewer than three ', text2, ' observations are non-zero


### Results from the BNB method (Denwood et al.) version C

Classification:  ', allres |> filter(Method=="BNB_KnownKs") |> pull(Classification), '\n
', allres |> filter(Method=="BNB_KnownKs") |> pull(PVstring), '\n

Notes:

- This method is parametric, and assumes that the data follow a negative binomial distribution:  the classification will be unavailable if the multiplication factor you entered does not match the data
- The over-dispersion is not estimated from the data, but is assumed to follow published estimates for typical over-dispersion in the host/parasite species you have selected
- This method may give misleading results in some groups of animals where the population over-dispersion is in fact different to published estimates
- This is the only preferred method when the sample size is less than 5 (and is the only viable method with a sample size of 1)

')

	}

	return(output)
}
